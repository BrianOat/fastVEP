# MC3 oncogenicity run

Runs fastVEP's somatic oncogenicity classifier (`annotate --somatic`) over every
FILTER=PASS mutation in the TCGA MC3 public call set
(<https://gdc.cancer.gov/about-data/publications/mc3-2017>,
`mc3.v0.2.8.PUBLIC.maf.gz`).

## One-time setup

Needs samtools, bcftools, htslib, picard, vcf2maf's maf2vcf.pl, python3 and unzip on PATH. The simplest way to get the tested versions is the conda environment: conda env create -f scripts/mc3_oncogenicity/environment.yml (any conda; Miniforge recommended)

```bash
# 1. Miniforge (conda + mamba): https://github.com/conda-forge/miniforge
# 2. The tool environment
mamba env create -f scripts/mc3_oncogenicity/environment.yml
conda activate fastvep-mc3
# 3. fastVEP
cargo build --release -p fastvep-cli
```

## Run

```bash
# Smoke test: the first 2,000 PASS rows, takes minutes rather than hours
scripts/mc3_oncogenicity/run_mc3.sh all --data-dir data --work runs/mc3_smoke --limit 2000

# Everything
scripts/mc3_oncogenicity/run_mc3.sh all --data-dir data --work runs/mc3
```

`all` runs four stages, each also callable on its own:

| Stage | What it does | Output |
|---|---|---|
| `prepare` | Downloads the MAF and references, keeps PASS rows, converts with vcf2maf's `maf2vcf.pl`, lifts GRCh37 → GRCh38 with Picard, re-checks every REF against GRCh38 | `mc3.grch38.vcf`, `mc3.targets.bed`, `funnel.json` |
| SA build | `scripts/build_sa_grch38.sh --targets mc3.targets.bed` (skipped when `--sa-dir` is given) | `sa_grch38/` with `MANIFEST.json` |
| `annotate` | `fastvep annotate --somatic --pick --output-format json` | `mc3.oncogenicity.json.gz` |
| `summarize` | `summarize.py` | `report.md`, `mc3.oncogenicity.tsv.gz` |

`--data-dir` holds the downloads (FASTAs, GFF3, REVEL, MAF) and can be shared
between runs; everything a run produces goes in `--work`. Every download and
build step skips itself when its output exists, so a rerun resumes.

## Decisions worth knowing

- **Liftover.** MC3 is GRCh37; fastVEP's annotation stack is GRCh38.
  Annotating across builds does not fail, it reads the wrong positions. The
  liftover rejects anything whose REF changes, rather than repairing it, and
  `bcftools norm --check-ref e` checks REF a second time and aborts the run on
  a mismatch. Picard is used rather than CrossMap: CrossMap 0.7.4 re-reads REF
  over one lifted base, which drops every deletion, and by default overwrites
  a mismatching REF instead of rejecting it.
- **Mitochondrial variants are dropped.** UCSC's hg19 `chrM` is not the rCRS
  sequence GRCh37 uses, so the hg19 chain cannot lift them.
- **One placeholder sample.** `maf2vcf.pl` writes one VCF column per
  tumor/normal pair, about 20,000 for MC3. The classifier needs each variant
  once, so `filter_maf.awk` gives every row the same placeholder pair and keeps
  the tumor count per variant in `mc3.counts.tsv`.
- **Joining back to MC3.** The GRCh37 key (`chr_pos_ref_alt` as the MAF writes
  it) goes into the VCF ID. It survives the liftover and comes back as the JSON
  `id`.
- **The SA directory is built for this run's targets.** gnomAD is pulled only
  around MC3's variants. Outside those regions a missing gnomAD record is not
  absence, and OP4 reads a missing record as absence. `annotate` refuses an SA
  directory whose `MANIFEST.json` names different targets (`--force`
  overrides).

## What to expect

Most variants will be VUS, and that is the classifier being honest rather than
broken. Hotspots (OS3/OM3/OP3), gene role (OVS1/OM2), variants already
classified under the SOP (OS1/OM4), functional domains (OM1), functional
studies (OS2/SBS2) and tumor type (OP2) have no source in fastVEP yet, so on
the oncogenic side a variant can reach about OP1 + OP4. The report's
criterion-outcome table measures how much of the SOP each run could ask, which
is the input for deciding which source to add first. This is an active TODO

The report's build-sanity section is the part to read first: REVEL should
cover nearly every missense variant, and a low number means the VCF and the
SA databases are on different builds.
