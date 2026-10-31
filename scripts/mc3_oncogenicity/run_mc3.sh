#!/usr/bin/env bash
# Run fastVEP's somatic oncogenicity classifier (--somatic) over TCGA MC3.
#
#   run_mc3.sh all --data-dir DATA --work WORK [--sa-dir SA] [--limit N]
#   run_mc3.sh prepare|annotate|summarize --data-dir DATA --work WORK [...]
#
# Stages
#   prepare    Download the MC3 public MAF (GRCh37), keep FILTER=PASS rows,
#              convert with vcf2maf's maf2vcf.pl, lift to GRCh38 with Picard
#              LiftoverVcf, re-check every REF against GRCh38, and write
#              WORK/mc3.targets.bed for the SA build.
#   annotate   fastvep annotate --somatic --pick, JSON, gzipped.
#   summarize  scripts/mc3_oncogenicity/summarize.py -> WORK/report.md and a
#              per-variant TSV.
#   all        prepare, build the SA directory with scripts/build_sa_grch38.sh
#              (unless --sa-dir is given), annotate, summarize.
#
# Options
#   --data-dir  Shared downloads; the same directory build_sa_grch38.sh uses.
#   --work      Outputs of this run.
#   --sa-dir    Existing GRCh38 SA directory. Default WORK/sa_grch38, built
#               by `all` for exactly this run's targets.
#   --limit N   Keep only the first N PASS rows (smoke test).
#   --fastvep   fastvep binary [target/release/fastvep in this repo].
#   --force     Annotate even if the SA directory's targets are not this run's.
#
# Tools: conda env from scripts/mc3_oncogenicity/environment.yml.
#
# Why each conversion step is the way it is:
#   - MC3 is GRCh37; fastVEP's annotation stack is GRCh38. Annotating across
#     builds does not fail, it looks up the wrong positions.
#   - maf2vcf.pl writes one VCF column per tumor/normal pair (~20,000 for
#     MC3). Only distinct variants matter here, so every row is given the same
#     placeholder pair; the per-variant tumor count is kept separately.
#   - The GRCh37 key (chr_pos_ref_alt as written in the MAF) goes into a
#     variant_id column, which maf2vcf.pl copies to the VCF ID. It survives
#     the liftover and reaches fastVEP's JSON `id`, so results join back to
#     MC3.

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
STAGE=${1:-}; shift || true
DATA="" WORK="" SA_DIR="" LIMIT=0 FASTVEP="$REPO/target/release/fastvep" FORCE=0

usage() { sed -n '2,36p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --data-dir) DATA=$2; shift 2 ;;
    --work) WORK=$2; shift 2 ;;
    --sa-dir) SA_DIR=$2; shift 2 ;;
    --limit) LIMIT=$2; shift 2 ;;
    --fastvep) FASTVEP=$2; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) usage 0 ;;
    *) echo "unknown argument: $1" >&2; usage 1 ;;
  esac
done
[[ $STAGE =~ ^(prepare|annotate|summarize|all)$ && -n $DATA && -n $WORK ]] || usage 1
mkdir -p "$DATA" "$WORK"
DATA=$(cd "$DATA" && pwd) WORK=$(cd "$WORK" && pwd)
SA_GIVEN=${SA_DIR:+1}
SA_DIR=${SA_DIR:-$WORK/sa_grch38}

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
need() { for t in "$@"; do command -v "$t" >/dev/null || { echo "missing tool: $t (conda activate fastvep-mc3)" >&2; exit 1; }; done; }
fetch() {
  local url=$1 dest=$2
  if [[ -s $dest ]]; then return; fi
  log "download $url"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}
gunzip_to() { [[ -s $2 ]] || { log "gunzip $(basename "$1")"; gunzip -c "$1" > "$2.part" && mv "$2.part" "$2"; }; }

MAF_URL=https://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc
GRCH37_URL=https://ftp.ensembl.org/pub/grch37/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
GRCH38_URL=https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
GFF3_URL=https://ftp.ensembl.org/pub/release-115/gff3/homo_sapiens/Homo_sapiens.GRCh38.115.gff3.gz
CHAIN_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

MAF=$DATA/mc3.v0.2.8.PUBLIC.maf.gz
GRCH37=$DATA/Homo_sapiens.GRCh37.dna.primary_assembly.fa
GRCH38=$DATA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GFF3=$DATA/Homo_sapiens.GRCh38.115.gff3
CHAIN=$DATA/hg19ToHg38.over.chain.gz
VCF38=$WORK/mc3.grch38.vcf
TARGETS=$WORK/mc3.targets.bed
COUNTS=$WORK/mc3.counts.tsv
FUNNEL=$WORK/funnel.json
JSON_OUT=$WORK/mc3.oncogenicity.json.gz

prepare() {
  need curl gunzip awk perl samtools bcftools picard maf2vcf.pl python3
  fetch "$MAF_URL" "$MAF"
  fetch "$CHAIN_URL" "$CHAIN"
  fetch "$GRCH37_URL" "$GRCH37.gz"; gunzip_to "$GRCH37.gz" "$GRCH37"
  fetch "$GRCH38_URL" "$GRCH38.gz"; gunzip_to "$GRCH38.gz" "$GRCH38"
  [[ -s $GRCH37.fai ]] || samtools faidx "$GRCH37"
  [[ -s $GRCH38.fai ]] || samtools faidx "$GRCH38"
  # Picard reads the target dictionary from <fasta minus .fa>.dict.
  [[ -s ${GRCH38%.fa}.dict ]] || samtools dict "$GRCH38" -o "${GRCH38%.fa}.dict"

  # ── PASS rows, one placeholder tumor/normal pair, GRCh37 key as variant_id
  local pass_maf=$WORK/mc3.pass.maf
  log "filter MAF to PASS (limit=${LIMIT})"
  gunzip -c "$MAF" | awk -F'\t' -v OFS='\t' -v limit="$LIMIT" \
      -v counts="$COUNTS" -v stats="$WORK/maf_stats.tsv" -f "$HERE/filter_maf.awk" > "$pass_maf.part"
  mv "$pass_maf.part" "$pass_maf"

  # ── MAF -> VCF (GRCh37) with vcf2maf's maf2vcf.pl ─────────────────────────
  local m2v=$WORK/maf2vcf vcf37=$WORK/mc3.grch37.vcf
  rm -rf "$m2v"
  log "maf2vcf.pl"
  maf2vcf.pl --input-maf "$pass_maf" --output-dir "$m2v" --output-vcf "$vcf37.raw" --ref-fasta "$GRCH37"
  # Sites only, and no MT: UCSC's hg19 chrM is not GRCh37's rCRS, so the
  # hg19 chain cannot lift it (its chains are dropped below too).
  bcftools view --no-version -G "$vcf37.raw" | awk -F'\t' '/^#/ || $1 != "MT"' > "$vcf37"
  rm -f "$vcf37.raw"

  # ── Liftover GRCh37 -> GRCh38 with Picard ─────────────────────────────────
  # UCSC's chain names contigs chr1..; Ensembl FASTAs and MC3 say 1... Picard
  # needs them to agree, so the chain is renamed rather than the FASTAs.
  local chain_ens=$WORK/hg19ToHg38.ensembl.chain
  zcat "$CHAIN" | awk '
    $1 == "chain" {
      keep = ($3 !~ /_/ && $8 !~ /_/ && $3 != "chrM" && $8 != "chrM")
      sub(/^chr/, "", $3); sub(/^chr/, "", $8)
    }
    keep' > "$chain_ens"
  local lifted=$WORK/mc3.grch38.lifted.vcf rejected=$WORK/mc3.liftover_rejected.vcf
  log "picard LiftoverVcf"
  # RECOVER_SWAPPED_REF_ALT stays off: if GRCh38 carries the tumor's ALT as
  # reference, swapping them would describe a different mutation.
  picard -Xmx8g LiftoverVcf -I "$vcf37" -O "$lifted" -CHAIN "$chain_ens" \
      -REJECT "$rejected" -R "$GRCH38" --WARN_ON_MISSING_CONTIG true \
      --WRITE_ORIGINAL_POSITION true --VALIDATION_STRINGENCY LENIENT \
      --CREATE_INDEX false

  # A second, independent REF check: abort on any mismatch rather than annotate it.
  log "sort and check REF against GRCh38"
  bcftools sort "$lifted" -Ou \
    | bcftools norm --no-version --check-ref e -f "$GRCH38" -Ov -o "$VCF38"

  awk -F'\t' -v OFS='\t' '!/^#/ { s = $2 - 51; if (s < 0) s = 0; print $1, s, $2 - 1 + length($4) + 50 }' "$VCF38" \
    | sort -k1,1 -k2,2n > "$TARGETS"

  # ── Funnel counts ─────────────────────────────────────────────────────────
  python3 - "$FUNNEL" "$WORK/maf_stats.tsv" "$m2v" "$vcf37" "$rejected" "$VCF38" <<'EOF'
import collections, glob, json, sys
out, stats, m2v, vcf37, rejected, vcf38 = sys.argv[1:]
f = {k: int(v) for k, v in (l.split("\t") for l in open(stats))}
skipped = glob.glob(f"{m2v}/*.skipped.tsv")
f["maf2vcf_skipped_ref_mismatch_rows"] = sum(1 for _ in open(skipped[0])) - 1 if skipped else 0
records = lambda p: sum(1 for l in open(p) if not l.startswith("#"))
f["grch37_vcf_records"] = records(vcf37)
reasons = collections.Counter(l.split("\t")[6] for l in open(rejected) if not l.startswith("#"))
f["liftover_rejected"] = dict(reasons)
f["grch38_vcf_records"] = records(vcf38)
json.dump(f, open(out, "w"), indent=2)
print(json.dumps(f, indent=2), file=sys.stderr)
EOF
  log "prepared $VCF38 and $TARGETS"
}

check_sa_targets() {
  local manifest=$SA_DIR/MANIFEST.json
  [[ -s $manifest ]] || { log "warning: $SA_DIR has no MANIFEST.json; cannot check its gnomAD targets"; return; }
  local want have
  want=$(sha256sum "$TARGETS" | cut -d' ' -f1)
  have=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["targets"]["sha256"])' "$manifest")
  if [[ $want != "$have" ]]; then
    echo "SA directory $SA_DIR was built for different targets than this run." >&2
    echo "Outside its targets a missing gnomAD record is not absence, so OP4 would be wrong there." >&2
    [[ $FORCE == 1 ]] || { echo "Rebuild it with these targets, or pass --force." >&2; exit 1; }
  fi
}

annotate() {
  [[ -x $FASTVEP ]] || { echo "fastvep not found at $FASTVEP (cargo build --release -p fastvep-cli)" >&2; exit 1; }
  [[ -s $VCF38 ]] || { echo "$VCF38 missing; run prepare first" >&2; exit 1; }
  [[ -d $SA_DIR ]] || { echo "SA directory $SA_DIR missing; run build_sa_grch38.sh or pass --sa-dir" >&2; exit 1; }
  [[ -s $GFF3 ]] || { fetch "$GFF3_URL" "$GFF3.gz"; gunzip_to "$GFF3.gz" "$GFF3"; }
  check_sa_targets
  log "fastvep annotate --somatic"
  # --fasta matters beyond HGVS: without it gnomAD lookups for indels are not
  # normalised, a miss is not trustworthy, and OP4 is reported not_evaluable.
  "$FASTVEP" annotate -i "$VCF38" -o - --gff3 "$GFF3" --fasta "$GRCH38" \
      --sa-dir "$SA_DIR" --somatic --pick --output-format json --no-progress \
    | gzip > "$JSON_OUT.part"
  mv "$JSON_OUT.part" "$JSON_OUT"
}

summarize() {
  python3 "$HERE/summarize.py" --json "$JSON_OUT" --counts "$COUNTS" --funnel "$FUNNEL" --out-dir "$WORK"
}

case $STAGE in
  prepare) prepare ;;
  annotate) annotate ;;
  summarize) summarize ;;
  all)
    prepare
    if [[ -z $SA_GIVEN ]]; then
      "$REPO/scripts/build_sa_grch38.sh" --targets "$TARGETS" --data-dir "$DATA" --out "$SA_DIR" --fastvep "$FASTVEP"
    fi
    annotate
    summarize
    ;;
esac
