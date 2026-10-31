//! Building an [`OncogenicityInput`] from what the annotation pipeline has
//! already attached to an allele.
//!
//! The record types (`GnomadData`, `RevelData`, `SpliceAiData`) are the germline
//! classifier's, reused rather than copied: they are the deserialised shape of
//! the supplementary-annotation JSON, and a second copy would drift the first
//! time a builder adds a field. What is *not* shared is interpretation. The
//! germline extractor derives PVS1, PM2 and BP7 signals from these records;
//! this one derives nothing germline-shaped, and the somatic criteria read the
//! raw records themselves.

use fastvep_classification::sa_extract::{GnomadData, RevelData, SpliceAiData};
use fastvep_core::{Consequence, SupplementaryAnnotation};

use crate::input::{OncogenicityInput, PopulationEvidence};

/// The transcript-level facts the pipeline computes for one allele.
#[derive(Debug, Clone, Copy, Default)]
pub struct AlleleTranscriptView<'a> {
    pub consequences: &'a [Consequence],
    pub gene_symbol: Option<&'a str>,
    pub is_canonical: bool,
    pub amino_acids: Option<&'a (String, String)>,
    pub protein_position: Option<u64>,
    pub hgvs_c: Option<&'a str>,
    pub escapes_nmd: Option<bool>,
}

/// Build the input the oncogenicity criteria read.
///
/// `gnomad_miss_is_absence` says whether a missing gnomAD record for this
/// allele can be trusted to mean the allele is absent from gnomAD: a gnomAD
/// database that holds this chromosome answered the lookup without error, and
/// the lookup used a representation gnomAD stores. Deciding that is the
/// caller's job, because only the caller saw the lookup. See
/// [`PopulationEvidence`].
///
/// The somatic-evidence fields of the result are left `None`; no source for
/// them exists in the pipeline yet.
pub fn extract_oncogenicity_input(
    view: &AlleleTranscriptView<'_>,
    allele_supplementary: &[(String, String)],
    variant_supplementary: &[SupplementaryAnnotation],
    gnomad_miss_is_absence: bool,
) -> OncogenicityInput {
    let mut gnomad: Option<GnomadData> = None;
    let mut revel: Option<RevelData> = None;
    let mut splice_ai: Option<SpliceAiData> = None;
    let mut phylop: Option<f64> = None;

    // Key spellings follow the germline extractor, which is where each alias
    // was found to be needed.
    for (key, json) in allele_supplementary {
        match key.as_str() {
            "gnomad" => gnomad = serde_json::from_str(json).ok(),
            "revel" => revel = serde_json::from_str(json).ok(),
            "spliceai" | "spliceAi" | "spliceAI" | "splice_ai" => {
                splice_ai = serde_json::from_str(json).ok();
            }
            "phylop" | "phyloP" if phylop.is_none() => phylop = parse_score(json),
            _ => {}
        }
    }
    // PhyloP is positional; depending on the path it arrives per allele or
    // per variant.
    if phylop.is_none() {
        phylop = variant_supplementary
            .iter()
            .find(|sa| matches!(sa.json_key.as_str(), "phylop" | "phyloP"))
            .and_then(|sa| parse_score(&sa.json_string));
    }

    let population = match (gnomad, gnomad_miss_is_absence) {
        (Some(record), _) => PopulationEvidence::Record(record),
        (None, true) => PopulationEvidence::NoRecord,
        (None, false) => PopulationEvidence::Unavailable,
    };

    OncogenicityInput {
        consequences: view.consequences.to_vec(),
        gene_symbol: view.gene_symbol.map(str::to_string),
        is_canonical: view.is_canonical,
        amino_acids: view.amino_acids.cloned(),
        protein_position: view.protein_position,
        hgvs_c: view.hgvs_c.map(str::to_string),
        escapes_nmd: view.escapes_nmd,
        population,
        revel,
        splice_ai,
        phylop,
        gene_role: None,
        hotspot: None,
        established_oncogenic: None,
        functional_domain: None,
        functional: None,
    }
}

fn parse_score(json: &str) -> Option<f64> {
    json.trim_matches('"').parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn view(consequences: &[Consequence]) -> AlleleTranscriptView<'_> {
        AlleleTranscriptView {
            consequences,
            ..Default::default()
        }
    }

    #[test]
    fn reads_the_records_the_germline_classifier_reads() {
        let sa = vec![
            (
                "gnomad".to_string(),
                r#"{"allAf":0.02,"afrAf":0.03}"#.to_string(),
            ),
            ("revel".to_string(), r#"{"score":0.91}"#.to_string()),
            ("spliceAI".to_string(), r#"{"dsAg":0.5}"#.to_string()),
            ("phylop".to_string(), "4.2".to_string()),
        ];
        let input =
            extract_oncogenicity_input(&view(&[Consequence::MissenseVariant]), &sa, &[], true);
        match &input.population {
            PopulationEvidence::Record(g) => assert_eq!(g.afr_af, Some(0.03)),
            other => panic!("expected a gnomAD record, got {other:?}"),
        }
        assert_eq!(input.revel.unwrap().score, Some(0.91));
        assert_eq!(input.splice_ai.unwrap().max_delta_score(), Some(0.5));
        assert_eq!(input.phylop, Some(4.2));
    }

    #[test]
    fn missing_gnomad_record_depends_on_whether_a_miss_can_be_trusted() {
        let v = view(&[Consequence::MissenseVariant]);
        assert!(matches!(
            extract_oncogenicity_input(&v, &[], &[], true).population,
            PopulationEvidence::NoRecord
        ));
        assert!(matches!(
            extract_oncogenicity_input(&v, &[], &[], false).population,
            PopulationEvidence::Unavailable
        ));
    }

    #[test]
    fn phylop_falls_back_to_the_variant_level_annotation() {
        let variant_sa = vec![SupplementaryAnnotation {
            json_key: "phyloP".to_string(),
            json_string: "\"1.5\"".to_string(),
            is_array: false,
        }];
        let input = extract_oncogenicity_input(
            &view(&[Consequence::SynonymousVariant]),
            &[],
            &variant_sa,
            false,
        );
        assert_eq!(input.phylop, Some(1.5));
    }

    #[test]
    fn somatic_evidence_is_never_filled_from_germline_sources() {
        // A germline-pathogenic ClinVar record and a LOF-intolerant gene are
        // present; neither may stand in for an established oncogenic variant
        // or a tumor-suppressor role.
        let sa = vec![(
            "clinvar".to_string(),
            r#"{"significance":["Pathogenic"]}"#.to_string(),
        )];
        let input = extract_oncogenicity_input(&view(&[Consequence::StopGained]), &sa, &[], true);
        assert!(input.gene_role.is_none());
        assert!(input.established_oncogenic.is_none());
        assert!(input.hotspot.is_none());
        assert!(input.functional_domain.is_none());
        assert!(input.functional.is_none());
    }
}
