#!/usr/bin/env bash
# Build a GRCh38 supplementary-annotation directory (--sa-dir) for fastVEP:
# gnomAD v4.1 exomes, REVEL v1.3, and SpliceAI + PhyloP distilled from gnomAD.
#
# This is docs/ACMG_SETUP.md's recipe as one resumable script. Every step
# skips itself when its output already exists, so rerunning after an
# interruption continues where it stopped.
#
#   scripts/build_sa_grch38.sh --targets regions.bed --data-dir DATA --out SA_DIR
#
# --targets  BED of the regions you will annotate (Ensembl "1" or UCSC "chr1"
#            names). gnomAD is pulled only for these regions.
# --data-dir Shared downloads (reference FASTA, GFF3, REVEL). Reused by
#            scripts/mc3_oncogenicity/run_mc3.sh.
# --out      The SA directory to write.
# --work     Scratch for gnomAD slices and extracts [<out>_work]. Tied to
#            the targets like --out is: never share it between target sets.
# --fastvep  fastvep binary [target/release/fastvep in this repo].
# --merge-gap  Merge target intervals closer than this many bp before
#            querying gnomAD [5000]. Fewer, larger remote reads; the extra
#            sequence only adds coverage.
#
# The resulting directory is only valid for variants inside --targets.
# Outside them gnomAD was never queried, so a missing gnomAD record there is
# not absence - and both PM2 (germline) and OP4 (somatic) read a missing
# record as absence. MANIFEST.json records the targets for that reason.
#
# Tools: curl, unzip, python3, samtools, bcftools, tabix, bgzip
# (scripts/mc3_oncogenicity/environment.yml provides them).

set -euo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TARGETS="" DATA="" OUT="" WORK="" FASTVEP="$REPO/target/release/fastvep" MERGE_GAP=5000

usage() { sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --targets) TARGETS=$2; shift 2 ;;
    --data-dir) DATA=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --work) WORK=$2; shift 2 ;;
    --fastvep) FASTVEP=$2; shift 2 ;;
    --merge-gap) MERGE_GAP=$2; shift 2 ;;
    -h|--help) usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done
[[ -n $TARGETS && -n $DATA && -n $OUT ]] || usage 1
[[ -s $TARGETS ]] || { echo "targets BED $TARGETS is missing or empty" >&2; exit 1; }
[[ -x $FASTVEP ]] || { echo "fastvep not found at $FASTVEP (cargo build --release -p fastvep-cli)" >&2; exit 1; }
for t in curl unzip python3 samtools bcftools tabix bgzip sha256sum; do
  command -v "$t" >/dev/null || { echo "missing tool: $t (conda activate fastvep-mc3)" >&2; exit 1; }
done
WORK=${WORK:-${OUT%/}_work}
mkdir -p "$DATA" "$OUT" "$WORK"
DATA=$(cd "$DATA" && pwd) OUT=$(cd "$OUT" && pwd) WORK=$(cd "$WORK" && pwd)
TARGETS=$(cd "$(dirname "$TARGETS")" && pwd)/$(basename "$TARGETS")

# Resuming skips any step whose output exists, and the gnomAD slices and
# databases are only right for the targets they were pulled for. So both
# directories are stamped with the targets they belong to, and a resume with
# different targets is refused rather than mixed.
TARGETS_SHA=$(sha256sum "$TARGETS" | cut -d' ' -f1)
for d in "$OUT" "$WORK"; do
  stamp=$d/.targets.sha256
  if [[ -s $stamp && $(cat "$stamp") != "$TARGETS_SHA" ]]; then
    echo "$d was built for different targets; use a new --out (and --work)." >&2
    exit 1
  fi
  echo "$TARGETS_SHA" > "$stamp"
done

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# Download once; a partial file never takes the final name.
fetch() {
  local url=$1 dest=$2
  if [[ -s $dest ]]; then return; fi
  log "download $url"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

built() { compgen -G "$OUT/$1.*" >/dev/null; }

ENSEMBL=https://ftp.ensembl.org/pub/release-115
FASTA_URL=$ENSEMBL/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GFF3_URL=$ENSEMBL/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
REVEL_URL=https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
# Plain http on purpose. Over https, libcurl negotiates HTTP/2 with Google
# Cloud Storage, and region queries through htslib fail within seconds with
# "Libcurl reported error 16 (Error in the HTTP2 framing layer)" (libcurl 8.21
# from the fastvep-mc3 env; reproducible). htslib has no setting to force
# HTTP/1.1, and plain http gets it. The data is public and every BGZF block
# carries a CRC32, so a corrupted block is rejected rather than read; what
# plain http does not protect against is deliberate tampering in transit.
GNOMAD_URL=http://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes/gnomad.exomes.v4.1.sites.chr%s.vcf.bgz
CHROMS=$(echo {1..22} X Y)

# ── 1. Reference sequence and gene models ──────────────────────────────────
FASTA=$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GFF3=$DATA/Homo_sapiens.GRCh38.115.gff3
for pair in "$FASTA_URL $FASTA" "$GFF3_URL $GFF3"; do
  set -- $pair
  if [[ ! -s $2 ]]; then
    fetch "$1" "$2.gz"
    log "gunzip $(basename "$2.gz")"
    gunzip -c "$2.gz" > "$2.part" && mv "$2.part" "$2"
  fi
done
[[ -s $FASTA.fai ]] || samtools faidx "$FASTA"

# ── 2. REVEL, split per chromosome to bound sa-build's memory ──────────────
REVEL_ZIP=$DATA/revel-v1.3_all_chromosomes.zip
fetch "$REVEL_URL" "$REVEL_ZIP"
REVEL_SPLIT=$DATA/revel_split   # independent of targets, so shared
if [[ ! -e $REVEL_SPLIT/.done ]]; then
  log "split REVEL per chromosome"
  rm -rf "$REVEL_SPLIT" && mkdir -p "$REVEL_SPLIT"
  unzip -p "$REVEL_ZIP" revel_with_transcript_ids \
    | awk -F, -v d="$REVEL_SPLIT" 'NR==1 {h=$0; next}
        { f = d "/revel_chr" $1 ".csv"
          if (!(f in seen)) { print h > f; seen[f] = 1 }
          print > f }'
  touch "$REVEL_SPLIT/.done"
fi
for f in "$REVEL_SPLIT"/revel_chr*.csv; do
  name=$(basename "$f" .csv)
  built "$name" && continue
  log "sa-build revel $name"
  "$FASTVEP" sa-build --source revel -i "$f" -o "$OUT/$name" --assembly GRCh38 --no-progress
done

# ── 3. gnomAD exomes for the target regions, then SpliceAI and PhyloP ──────
# One remote pull per chromosome into a local slice; gnomAD, SpliceAI and
# PhyloP are all built from that slice, so the bucket is read once.
SLICES=$WORK/gnomad_slices EXTRACTS=$WORK/extracts BEDS=$WORK/target_beds
mkdir -p "$SLICES" "$EXTRACTS" "$BEDS"
for c in $CHROMS; do
  bed=$BEDS/chr$c.bed
  if [[ ! -e $bed ]]; then
    # gnomAD names contigs chr-prefixed; accept either style in --targets.
    awk -v c="chr$c" -v n="$c" -v gap="$MERGE_GAP" 'BEGIN {OFS="\t"}
        ($1 == c || $1 == n) {print c, $2, $3}' "$TARGETS" \
      | sort -k2,2n \
      | awk -v gap="$MERGE_GAP" 'BEGIN {OFS="\t"}
          NR == 1 {c=$1; s=$2; e=$3; next}
          $2 <= e + gap {if ($3 > e) e = $3; next}
          {print c, s, e; s=$2; e=$3}
          END {if (NR) print c, s, e}' > "$bed.part"
    mv "$bed.part" "$bed"
  fi
  if [[ ! -s $bed ]]; then
    log "chr$c: no targets, skipped"
    continue
  fi

  slice=$SLICES/gnomad.exomes.v4.1.chr$c.targets.vcf.gz
  if [[ ! -s $slice.tbi ]]; then
    log "chr$c: pull gnomAD for $(wc -l < "$bed") regions"
    url=$(printf "$GNOMAD_URL" "$c")
    # htslib caches a remote index in the working directory; keep it in the work dir.
    (cd "$SLICES" && bcftools view --no-version -R "$bed" -Oz -o "$slice.part" "$url")
    mv "$slice.part" "$slice"
    tabix -f -p vcf "$slice"
  fi
  if [[ $(bcftools index -n "$slice") -eq 0 ]]; then
    log "chr$c: gnomAD has no records in the targets"
    continue
  fi

  built "gnomad_chr$c" || {
    log "sa-build gnomad chr$c"
    "$FASTVEP" sa-build --source gnomad -i "$slice" -o "$OUT/gnomad_chr$c" --assembly GRCh38 --no-progress
  }
  if ! built "spliceai_chr$c" || ! built "phylop_chr$c"; then
    python3 "$REPO/scripts/extract_gnomad_scores.py" --chrom "$c" --url "$slice" --out-dir "$EXTRACTS"
    built "spliceai_chr$c" || "$FASTVEP" sa-build --source spliceai -i "$EXTRACTS/spliceai_chr$c.vcf" \
        -o "$OUT/spliceai_chr$c" --assembly GRCh38 --no-progress
    built "phylop_chr$c" || "$FASTVEP" sa-build --source phylop -i "$EXTRACTS/phylop_chr$c.tsv" \
        -o "$OUT/phylop_chr$c" --assembly GRCh38 --no-progress
    rm -f "$EXTRACTS/spliceai_chr$c.vcf" "$EXTRACTS/phylop_chr$c.tsv"
  fi
done

# ── 4. Manifest ────────────────────────────────────────────────────────────
log "write MANIFEST.json"
python3 - "$OUT/MANIFEST.json" <<EOF
import hashlib, json, os, subprocess, sys, datetime

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def git(*args):
    try:
        return subprocess.check_output(["git", "-C", "$REPO", *args], text=True).strip()
    except Exception:
        return None

manifest = {
    "assembly": "GRCh38",
    "built": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
    "fastvep_commit": git("rev-parse", "HEAD"),
    "fastvep_dirty": bool(git("status", "--porcelain")),
    "valid_only_inside_targets": (
        "gnomAD was queried only inside the targets below. Outside them a "
        "missing gnomAD record is not absence, so OP4 and PM2 must not be "
        "read for variants there."
    ),
    "targets": {"path": "$TARGETS", "sha256": sha256("$TARGETS"), "merge_gap_bp": $MERGE_GAP},
    "sources": {
        "gnomad": {"url": "$GNOMAD_URL".replace("%s", "{chrom}"), "release": "4.1 exomes"},
        "spliceai": {"from": "gnomAD 4.1 INFO/spliceai_ds_max (scripts/extract_gnomad_scores.py)"},
        "phylop": {"from": "gnomAD 4.1 INFO/phylop, Zoonomia (scripts/extract_gnomad_scores.py)"},
        "revel": {"url": "$REVEL_URL", "sha256": sha256("$REVEL_ZIP")},
    },
    "databases": sorted(f for f in os.listdir("$OUT") if not f.startswith(".") and f != "MANIFEST.json"),
}
json.dump(manifest, open(sys.argv[1], "w"), indent=2)
EOF
log "done: $OUT"
