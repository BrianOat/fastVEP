//! End-to-end check that `--somatic` reaches the output, and that it reads
//! gnomAD the way OP4 needs.
//!
//! OP4 ("absent from controls") turns a *missing* gnomAD record into evidence,
//! so the one thing this pipeline must not do is treat a lookup that could not
//! have succeeded as a lookup that found nothing. Without a FASTA the SA loop
//! queries an indel in its raw VCF spelling, which can miss the left-aligned
//! record gnomAD stores; such a miss has to stay "not evaluable". A
//! substitution has one spelling and is safe either way. Nor may a miss count
//! when no loaded gnomAD database holds the variant's chromosome at all.
//!
//! It also pins down that `--somatic` is additive: the germline ACMG output
//! must not change when it is turned on.

use fastvep_cli::pipeline::sa_build::run_sa_build;
use fastvep_cli::pipeline::{run_annotate, AnnotateConfig};
use serde_json::Value;
use std::fs;
use std::path::Path;

// Single-exon coding transcript on contig `1`: ATG TGG CGG TAA (Met Trp Arg Ter).
const GFF3: &str = "##gff-version 3\n\
1\ttest\tgene\t1\t12\t.\t+\t.\tID=gene:ENSG_ONC1;Name=ONCTEST;gene_name=ONCTEST;biotype=protein_coding\n\
1\ttest\tmRNA\t1\t12\t.\t+\t.\tID=transcript:ENST_ONC1;Parent=gene:ENSG_ONC1;biotype=protein_coding\n\
1\ttest\texon\t1\t12\t.\t+\t.\tID=exon:ENSE_ONC1;Parent=transcript:ENST_ONC1;rank=1\n\
1\ttest\tCDS\t1\t12\t.\t+\t0\tID=CDS:ENSP_ONC1;Parent=transcript:ENST_ONC1\n";

const FASTA: &str = ">1\nATGTGGCGGTAA\n";

// `common`: Trp2Ser, in gnomAD at 2% in NFE, so SBS1.
// `absent_snv`: Arg3Gln, not in gnomAD, so OP4.
// `absent_indel`: in-frame insertion, not in gnomAD; OP4 only with a FASTA.
const VCF: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t5\tcommon\tG\tC\t.\tPASS\t.\n\
1\t6\tabsent_indel\tG\tGCGG\t.\tPASS\t.\n\
1\t8\tabsent_snv\tG\tA\t.\tPASS\t.\n";

const GNOMAD: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t5\t.\tG\tC\t.\tPASS\tAF=0.02;AC=2000;AN=100000;nhomalt=10;AF_nfe=0.02;AF_afr=0.001\n";

fn config(dir: &Path, output: &Path, format: &str, fasta: bool, somatic: bool) -> AnnotateConfig {
    AnnotateConfig {
        input: dir.join("in.vcf").to_string_lossy().into(),
        output: output.to_string_lossy().into(),
        gff3: vec![dir.join("test.gff3").to_string_lossy().into()],
        fasta: fasta.then(|| dir.join("test.fa").to_string_lossy().into()),
        output_format: format.into(),
        pick: false,
        hgvs: false,
        distance: 0,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: Some(dir.join("sa").to_string_lossy().into()),
        sa_only: false,
        acmg: false,
        acmg_config: None,
        pick_order: None,
        functional_evidence: None,
        proband: None,
        mother: None,
        father: None,
        gene_list: None,
        explicit_alleles: false,
        qc_rules: None,
        show_progress: false,
        somatic,
    }
}

fn setup() -> tempfile::TempDir {
    setup_with_gnomad(GNOMAD)
}

/// The fixture with a gnomAD database built from `gnomad_vcf`.
fn setup_with_gnomad(gnomad_vcf: &str) -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    fs::write(dir.join("test.gff3"), GFF3).unwrap();
    fs::write(dir.join("test.fa"), FASTA).unwrap();
    fs::write(dir.join("in.vcf"), VCF).unwrap();
    fs::write(dir.join("gnomad.vcf"), gnomad_vcf).unwrap();
    fs::create_dir(dir.join("sa")).unwrap();
    run_sa_build(
        "gnomad",
        dir.join("gnomad.vcf").to_str().unwrap(),
        dir.join("sa").join("gnomad").to_str().unwrap(),
        "GRCh38",
        None,
        &[],
        false,
    )
    .unwrap();
    tmp
}

/// The `oncogenicity` object of the one transcript consequence of variant `id`.
fn oncogenicity<'a>(variants: &'a [Value], id: &str) -> &'a Value {
    let v = variants
        .iter()
        .find(|v| v["id"] == id)
        .unwrap_or_else(|| panic!("no variant {id}"));
    &v["transcript_consequences"][0]["oncogenicity"]
}

fn outcome<'a>(onc: &'a Value, code: &str) -> &'a str {
    onc["criteria"]
        .as_array()
        .unwrap()
        .iter()
        .find(|c| c["code"] == code)
        .unwrap_or_else(|| panic!("no {code} in {onc}"))["outcome"]
        .as_str()
        .unwrap()
}

fn annotate_json(fasta: bool, somatic: bool) -> Vec<Value> {
    let tmp = setup();
    let out = tmp.path().join("out.json");
    run_annotate(config(tmp.path(), &out, "json", fasta, somatic)).unwrap();
    serde_json::from_str(&fs::read_to_string(&out).unwrap()).unwrap()
}

/// Annotate the fixture with `--acmg`, with or without `--somatic`, and return
/// the raw output.
fn annotate_acmg(format: &str, somatic: bool) -> String {
    let tmp = setup();
    let out = tmp.path().join(format!("out.{format}"));
    let mut cfg = config(tmp.path(), &out, format, true, somatic);
    cfg.acmg = true;
    run_annotate(cfg).unwrap();
    fs::read_to_string(&out).unwrap()
}

#[test]
fn somatic_result_reaches_json_with_every_criterion() {
    let variants = annotate_json(true, true);
    let onc = oncogenicity(&variants, "absent_snv");
    assert_eq!(onc["criteria"].as_array().unwrap().len(), 17);
    // Premise: gnomAD was loaded and consulted, so the absent SNV is absent.
    assert_eq!(outcome(onc, "OP4"), "met");
    // No cancer hotspot source exists: unknown, not negative.
    assert_eq!(outcome(onc, "OS3"), "not_evaluable");
    assert_eq!(onc["points"], 1);
    assert_eq!(onc["shorthand"], "VUS");
}

#[test]
fn common_variant_gets_benign_frequency_evidence() {
    let variants = annotate_json(true, true);
    let onc = oncogenicity(&variants, "common");
    assert_eq!(outcome(onc, "SBS1"), "met");
    assert_eq!(outcome(onc, "SBVS1"), "not_met");
    assert_eq!(outcome(onc, "OP4"), "not_met");
    assert_eq!(onc["shorthand"], "LB");
}

#[test]
fn unnormalised_indel_miss_is_not_absence() {
    let without_fasta = annotate_json(false, true);
    assert_eq!(
        outcome(oncogenicity(&without_fasta, "absent_indel"), "OP4"),
        "not_evaluable"
    );
    // A substitution needs no normalising, so it is still absent without one.
    assert_eq!(
        outcome(oncogenicity(&without_fasta, "absent_snv"), "OP4"),
        "met"
    );

    let with_fasta = annotate_json(true, true);
    assert_eq!(
        outcome(oncogenicity(&with_fasta, "absent_indel"), "OP4"),
        "met"
    );
}

#[test]
fn without_flag_nothing_is_classified() {
    let variants = annotate_json(true, false);
    assert!(variants.iter().all(|v| v["transcript_consequences"][0]
        .get("oncogenicity")
        .is_none()));
}

#[test]
fn somatic_result_reaches_vcf_csq() {
    let tmp = setup();
    let out = tmp.path().join("out.vcf");
    run_annotate(config(tmp.path(), &out, "vcf", true, true)).unwrap();
    let vcf = fs::read_to_string(&out).unwrap();
    assert!(vcf.contains("|ONCOGENICITY|ONCOGENICITY_CRITERIA"));
    let line = vcf.lines().find(|l| l.contains("\tabsent_snv\t")).unwrap();
    assert!(line.contains("|VUS|OP4"), "{line}");
}

#[test]
fn miss_on_a_chromosome_no_gnomad_database_holds_is_not_absence() {
    // gnomAD databases are built one per chromosome, and a database asked
    // about a chromosome it does not hold answers "no record". Here the only
    // gnomAD database holds contig 2 and every variant is on contig 1, so no
    // lookup could have found anything - OP4 must not read that as absence.
    let gnomad_chr2_only = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
2\t100\t.\tA\tG\t.\tPASS\tAF=0.001;AC=100;AN=100000;nhomalt=0;AF_nfe=0.001\n";
    let tmp = setup_with_gnomad(gnomad_chr2_only);
    let out = tmp.path().join("out.json");
    run_annotate(config(tmp.path(), &out, "json", true, true)).unwrap();
    let variants: Vec<Value> = serde_json::from_str(&fs::read_to_string(&out).unwrap()).unwrap();
    for id in ["absent_snv", "absent_indel", "common"] {
        let onc = oncogenicity(&variants, id);
        assert_eq!(outcome(onc, "OP4"), "not_evaluable", "{id}");
        assert_eq!(outcome(onc, "SBS1"), "not_evaluable", "{id}");
    }
}

/// Remove every `oncogenicity` key, wherever it sits.
fn strip_oncogenicity(v: &mut Value) {
    match v {
        Value::Object(map) => {
            map.remove("oncogenicity");
            map.values_mut().for_each(strip_oncogenicity);
        }
        Value::Array(items) => items.iter_mut().for_each(strip_oncogenicity),
        _ => {}
    }
}

#[test]
fn somatic_leaves_acmg_json_output_unchanged() {
    let mut without: Value = serde_json::from_str(&annotate_acmg("json", false)).unwrap();
    let mut with: Value = serde_json::from_str(&annotate_acmg("json", true)).unwrap();

    // Premise: ACMG ran and produced evidence, so equality is not vacuous.
    let acmg_met = without
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|v| v["transcript_consequences"].as_array().unwrap())
        .flat_map(|tc| tc["acmg"]["criteria"].as_array().unwrap())
        .filter(|c| c["met"] == true)
        .count();
    assert!(acmg_met > 0, "fixture produced no met ACMG criteria");
    // And --somatic did add its own field.
    assert!(with[0]["transcript_consequences"][0]
        .get("oncogenicity")
        .is_some());

    strip_oncogenicity(&mut without);
    strip_oncogenicity(&mut with);
    assert_eq!(
        without, with,
        "--somatic changed output other than its own field"
    );
}

#[test]
fn somatic_leaves_acmg_vcf_output_unchanged() {
    let without = annotate_acmg("vcf", false);
    let with = annotate_acmg("vcf", true);

    let header = without
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .unwrap();
    let names: Vec<&str> = header
        .split("Format: ")
        .nth(1)
        .unwrap()
        .trim_end_matches("\">")
        .split('|')
        .collect();
    let own: Vec<usize> = ["ONCOGENICITY", "ONCOGENICITY_CRITERIA"]
        .iter()
        .map(|n| names.iter().position(|x| x == n).unwrap())
        .collect();
    let acmg = names.iter().position(|x| *x == "ACMG").unwrap();

    // Every CSQ column but the two --somatic owns, per variant line.
    let rows = |vcf: &str| -> Vec<Vec<Vec<String>>> {
        vcf.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let info = l.split('\t').nth(7).unwrap();
                let csq = info
                    .split(';')
                    .find_map(|kv| kv.strip_prefix("CSQ="))
                    .unwrap();
                csq.split(',')
                    .map(|entry| {
                        entry
                            .split('|')
                            .enumerate()
                            .filter(|(i, _)| !own.contains(i))
                            .map(|(_, f)| f.to_string())
                            .collect()
                    })
                    .collect()
            })
            .collect()
    };
    let (a, b) = (rows(&without), rows(&with));
    assert!(
        a.iter().flatten().any(|entry| !entry[acmg].is_empty()),
        "fixture produced no ACMG classification"
    );
    assert_eq!(a, b, "--somatic changed a CSQ column other than its own");
    // The non-CSQ part of the header and records must match too.
    let non_csq = |vcf: &str| -> Vec<String> {
        vcf.lines()
            .filter(|l| !l.starts_with("##INFO=<ID=CSQ"))
            .map(|l| l.split("CSQ=").next().unwrap().to_string())
            .collect()
    };
    assert_eq!(non_csq(&without), non_csq(&with));
}
