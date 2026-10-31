#!/usr/bin/env python3
"""Summarise fastVEP --somatic output over MC3 into report.md and a TSV.

Reads fastVEP's JSON (a pretty-printed array, optionally gzipped) as a
stream, so the whole run never has to fit in memory, and joins each record
back to MC3 through its VCF ID, which run_mc3.sh sets to the GRCh37 MAF key.

Standard library only.
"""

import argparse
import collections
import gzip
import json
import os
import sys

CRITERIA = [
    "OVS1", "OS1", "OS2", "OS3", "OM1", "OM2", "OM3", "OM4", "OP1", "OP2",
    "OP3", "OP4", "SBVS1", "SBS1", "SBS2", "SBP1", "SBP2",
]
OUTCOMES = ["met", "not_met", "not_evaluable", "excluded"]
CLASSES = ["O", "LO", "VUS", "LB", "B"]

NULL_TERMS = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost",
}
LENGTH_TERMS = {"inframe_insertion", "inframe_deletion", "stop_lost", "protein_altering_variant"}

# MC3's own consequence (VEP on GENCODE v19, its own transcript choice) and
# the fastVEP term that should accompany it. Disagreement is expected at some
# rate - different release, different transcript - so this is a sanity
# signal for the liftover and the annotation, not a test.
MC3_EXPECTED = {
    "Missense_Mutation": {"missense_variant"},
    "Silent": {"synonymous_variant"},
    "Nonsense_Mutation": {"stop_gained"},
    "Frame_Shift_Del": {"frameshift_variant"},
    "Frame_Shift_Ins": {"frameshift_variant"},
    "In_Frame_Del": {"inframe_deletion"},
    "In_Frame_Ins": {"inframe_insertion"},
    "Splice_Site": {"splice_donor_variant", "splice_acceptor_variant"},
    "Nonstop_Mutation": {"stop_lost"},
    "Translation_Start_Site": {"start_lost"},
}


def iter_json_array(path):
    """Yield the objects of a JSON array without loading the array."""
    opener = gzip.open if path.endswith(".gz") else open
    decoder = json.JSONDecoder()
    buf, pos, eof = "", 0, False
    with opener(path, "rt") as f:
        while True:
            while pos < len(buf) and buf[pos] in " \t\r\n,[]":
                pos += 1
            if pos >= len(buf):
                if eof:
                    return
                buf, pos = f.read(1 << 22), 0
                eof = not buf
                continue
            try:
                obj, end = decoder.raw_decode(buf, pos)
            except json.JSONDecodeError:
                more = f.read(1 << 22)
                if not more:
                    raise
                buf, pos = buf[pos:] + more, 0
                continue
            yield obj
            pos = end


def group(terms):
    t = set(terms)
    if t & NULL_TERMS:
        return "predicted_null"
    if t & LENGTH_TERMS:
        return "inframe_or_stop_lost"
    if "missense_variant" in t:
        return "missense"
    if "synonymous_variant" in t:
        return "synonymous"
    return "other"


def pct(a, b):
    return f"{100 * a / b:.1f}%" if b else "n/a"


def table(header, rows):
    out = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--json", required=True)
    ap.add_argument("--counts", required=True, help="mc3.counts.tsv from run_mc3.sh prepare")
    ap.add_argument("--funnel", help="funnel.json from run_mc3.sh prepare")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    mc3 = {}
    with open(args.counts) as f:
        next(f)
        for line in f:
            key, n, gene, cls, hgvsp = (line.rstrip("\n").split("\t") + [""] * 5)[:5]
            mc3[key] = (int(n), gene, cls, hgvsp)

    by_class = collections.Counter()
    by_criterion = {c: collections.Counter() for c in CRITERIA}
    class_by_group = collections.defaultdict(collections.Counter)
    coverage = collections.defaultdict(collections.Counter)
    agreement = collections.defaultdict(collections.Counter)
    recurrence = collections.defaultdict(collections.Counter)
    recurrent_rows = []
    variants = no_tc = unjoined = 0

    tsv_path = os.path.join(args.out_dir, "mc3.oncogenicity.tsv.gz")
    with gzip.open(tsv_path, "wt") as tsv:
        tsv.write("\t".join([
            "grch37_key", "chrom", "start", "allele_string", "gene", "consequence",
            "mc3_gene", "mc3_class", "mc3_hgvsp", "n_tumors", "class", "points",
            "met", "excluded", "not_evaluable", "gnomad", "revel", "spliceai", "phylop",
        ]) + "\n")
        for v in iter_json_array(args.json):
            variants += 1
            key = v.get("id", "")
            n, mgene, mcls, mhgvsp = mc3.get(key, (0, "", "", ""))
            if key not in mc3:
                unjoined += 1
            tcs = [tc for tc in v.get("transcript_consequences", []) if "oncogenicity" in tc]
            if not tcs:
                no_tc += 1
                continue
            for tc in tcs:
                onc = tc["oncogenicity"]
                cls = onc["shorthand"]
                terms = tc.get("consequence_terms", [])
                g = group(terms)
                by_class[cls] += 1
                class_by_group[g][cls] += 1
                outcome = {}
                for c in onc["criteria"]:
                    outcome[c["code"]] = c["outcome"]
                    by_criterion[c["code"]][c["outcome"]] += 1

                # Positional sources (PhyloP) can sit at the variant level.
                has = {
                    k: any(s in tc or s in v for s in names)
                    for k, names in (
                        ("gnomad", ("gnomad",)),
                        ("revel", ("revel",)),
                        ("spliceai", ("spliceai", "spliceAI", "spliceAi", "splice_ai")),
                        ("phylop", ("phylop", "phyloP")),
                    )
                }
                for k, present in has.items():
                    coverage[g][k] += present
                coverage[g]["n"] += 1

                if mcls in MC3_EXPECTED:
                    agreement[mcls]["n"] += 1
                    agreement[mcls]["agree"] += bool(MC3_EXPECTED[mcls] & set(terms))

                band = "1" if n <= 1 else "2-4" if n < 5 else "5-9" if n < 10 else ">=10"
                recurrence[band][cls] += 1
                codes = lambda o: ",".join(c for c in CRITERIA if outcome.get(c) == o)
                row = [
                    key, v.get("seq_region_name", ""), v.get("start", ""),
                    v.get("allele_string", ""), tc.get("gene_symbol", ""), ",".join(terms),
                    mgene, mcls, mhgvsp, n, cls, onc["points"], codes("met"),
                    codes("excluded"), codes("not_evaluable"),
                    *(int(has[k]) for k in ("gnomad", "revel", "spliceai", "phylop")),
                ]
                tsv.write("\t".join(map(str, row)) + "\n")
                if n >= 5:
                    recurrent_rows.append(row)

    total = sum(by_class.values())
    lines = ["# fastVEP --somatic on TCGA MC3", ""]

    if args.funnel and os.path.exists(args.funnel):
        f = json.load(open(args.funnel))
        lines += ["## Funnel", "", table(["stage", "count"], [
            ["MAF rows", f["maf_rows"]],
            ["FILTER=PASS rows", f["pass_rows"]],
            ["distinct PASS variants (GRCh37 key)", f["pass_unique_variants"]],
            ["maf2vcf.pl REF-mismatch rows (skipped)", f["maf2vcf_skipped_ref_mismatch_rows"]],
            ["GRCh37 VCF records (MT removed)", f["grch37_vcf_records"]],
            ["liftover rejected", sum(f["liftover_rejected"].values())],
            ["GRCh38 VCF records (REF verified)", f["grch38_vcf_records"]],
            ["variants in fastVEP JSON", variants],
            ["  with no transcript consequence (intergenic under --pick)", no_tc],
            ["  whose ID did not join back to MC3", unjoined],
            ["classified allele-transcript rows", total],
        ]), ""]
        if f["liftover_rejected"]:
            lines += ["Liftover rejections by Picard reason: "
                      + ", ".join(f"{k} {v}" for k, v in sorted(f["liftover_rejected"].items())), ""]

    lines += ["## Classification", "", table(
        ["class", "count", "share"], [[c, by_class[c], pct(by_class[c], total)] for c in CLASSES]), ""]

    lines += ["## Criterion outcomes", "",
              "`not_evaluable` means the evidence the criterion needs is not available to this run; "
              "it is never counted as evidence either way.", "",
              table(["criterion", *OUTCOMES], [[c, *(by_criterion[c][o] for o in OUTCOMES)] for c in CRITERIA]), ""]

    groups = ["missense", "predicted_null", "inframe_or_stop_lost", "synonymous", "other"]
    lines += ["## Class by consequence", "", table(
        ["consequence", "n", *CLASSES],
        [[g, sum(class_by_group[g].values()), *(class_by_group[g][c] for c in CLASSES)] for g in groups]), ""]

    lines += ["## Annotation coverage (build sanity)", "",
              "REVEL scores every possible missense SNV, so its coverage of missense variants should be "
              "high. A low value is the signature of a build mismatch between the VCF and the SA "
              "databases. SpliceAI and PhyloP come from gnomAD records, and most somatic variants "
              "are absent from gnomAD, so their coverage is expected to be low.", "",
              table(["consequence", "n", "gnomAD", "REVEL", "SpliceAI", "PhyloP"], [
                  [g, coverage[g]["n"], *(pct(coverage[g][k], coverage[g]["n"]) for k in ("gnomad", "revel", "spliceai", "phylop"))]
                  for g in groups]), ""]
    mis = coverage["missense"]
    if mis["n"] >= 100 and mis["revel"] / mis["n"] < 0.5:
        lines += [f"**WARNING: REVEL covers only {pct(mis['revel'], mis['n'])} of missense variants. "
                  "Check that the VCF and the SA directory are on the same build.**", ""]

    lines += ["## Consequence agreement with MC3", "",
              "MC3 annotated with VEP on GENCODE v19 and its own transcript choice; this run uses "
              "Ensembl 115 and `--pick`. Some disagreement is expected; a collapse is not.", "",
              table(["MC3 class", "n", "fastVEP agrees"],
                    [[k, a["n"], pct(a["agree"], a["n"])] for k, a in sorted(agreement.items())]), ""]

    lines += ["## Recurrence vs class", "",
              "Variants seen in many tumors are hotspot-like. They should rarely be benign; any that "
              "are are listed below for review.", "",
              table(["tumors", "n", *CLASSES],
                    [[b, sum(recurrence[b].values()), *(recurrence[b][c] for c in CLASSES)]
                     for b in ("1", "2-4", "5-9", ">=10")]), ""]
    benign = [r for r in recurrent_rows if r[10] in ("LB", "B")]
    benign.sort(key=lambda r: -r[9])
    if benign:
        lines += ["### Recurrent (>=5 tumors) variants called LB/B", "", table(
            ["gene", "MC3 HGVSp", "tumors", "class", "points", "met"],
            [[r[4] or r[6], r[8], r[9], r[10], r[11], r[12]] for r in benign[:50]]), ""]
    top = sorted(recurrent_rows, key=lambda r: -r[9])[:30]
    if top:
        lines += ["### Most recurrent variants", "", table(
            ["gene", "MC3 HGVSp", "tumors", "class", "points", "met", "not evaluable"],
            [[r[4] or r[6], r[8], r[9], r[10], r[11], r[12], r[14]] for r in top]), ""]

    report = os.path.join(args.out_dir, "report.md")
    with open(report, "w") as f:
        f.write("\n".join(lines))
    print(f"wrote {report} and {tsv_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
