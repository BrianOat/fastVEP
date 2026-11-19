//! A structural variant is not given a somatic oncogenicity classification.
//!
//! The SOP (Horak 2022) classifies sequence changes, and a structural variant has the terms of one
//! in its rows: `frameshift_variant` for a deletion inside an exon, `start_lost` for an inversion
//! over a start codon. The criteria would read those as what they are for a base, and a symbolic
//! allele such as `<DEL>` can never match a gnomAD record, so a gnomAD miss for it would read as
//! "absent from controls" (OP4). The ACMG pass skips structural variants for the first reason.

use std::path::{Path, PathBuf};

use fastvep_cli::pipeline::{run_annotate, AnnotateConfig, ColumnFlags, PickFlags};
use tempfile::TempDir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/vep_tab_sv")
        .join(name)
}

/// The base of `chr1` at a 1-based position, from the fixture's FASTA (lines of 60).
fn base(position: usize) -> char {
    let text = std::fs::read_to_string(fixture("genome.fa")).unwrap();
    let sequence: String = text.lines().skip(1).collect();
    sequence.as_bytes()[position - 1] as char
}

#[test]
fn a_structural_variant_has_no_oncogenicity_and_a_base_change_beside_it_has() {
    let dir = TempDir::new().unwrap();
    let (gff3, fasta) = (dir.path().join("genes.gff3"), dir.path().join("genome.fa"));
    std::fs::copy(fixture("genes.gff3"), &gff3).unwrap();
    std::fs::copy(fixture("genome.fa"), &fasta).unwrap();

    // a deletion inside the first coding exon of SVA, an inversion over its start codon, and a
    // change of one base in the same exon
    let reference = base(1060);
    let alt = if reference == 'A' { 'C' } else { 'A' };
    let vcf = dir.path().join("in.vcf");
    std::fs::write(
        &vcf,
        format!(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t1099\tdeletion\t{}\t<DEL>\t.\t.\tEND=1101\n\
             chr1\t1049\tinversion\t{}\t<INV>\t.\t.\tEND=1053\n\
             chr1\t1060\tsnv\t{reference}\t{alt}\t.\t.\t.\n",
            base(1099),
            base(1049),
        ),
    )
    .unwrap();
    let out = dir.path().join("out.json");
    run_annotate(AnnotateConfig {
        input: vcf.to_string_lossy().into(),
        output: out.to_string_lossy().into(),
        gff3: vec![gff3.to_string_lossy().into()],
        fasta: Some(fasta.to_string_lossy().into()),
        output_format: "json".into(),
        pick: PickFlags::default(),
        hgvs: true,
        columns: ColumnFlags::default(),
        distance: 1000,
        cache_dir: None,
        transcript_cache: None,
        sa_dir: None,
        sa_only: false,
        acmg: false,
        somatic: true,
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
        compat: Default::default(),
        exclude: Default::default(),
    })
    .unwrap();

    let variants: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&out).unwrap()).unwrap();
    let rows = |id: &str| -> Vec<serde_json::Value> {
        variants
            .as_array()
            .unwrap()
            .iter()
            .filter(|v| v["id"] == id)
            .flat_map(|v| v["transcript_consequences"].as_array().unwrap().clone())
            .collect()
    };
    let classified = |id: &str| rows(id).iter().any(|t| t.get("oncogenicity").is_some());
    // the terms the criteria would have read as a base change's
    let terms = |id: &str| -> String {
        rows(id)
            .iter()
            .map(|t| t["consequence_terms"].to_string())
            .collect()
    };
    assert!(terms("deletion").contains("frameshift_variant"));
    assert!(terms("inversion").contains("start_lost"));
    assert!(!classified("deletion"), "a deletion was classified");
    assert!(!classified("inversion"), "an inversion was classified");
    assert!(classified("snv"), "the change of a base is classified");
}
