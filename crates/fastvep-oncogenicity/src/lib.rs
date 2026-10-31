//! Somatic oncogenicity classification for fastVEP.
//!
//! Implements the ClinGen / CGC / VICC standard (Horak et al. 2022, PMC9081216):
//! seventeen criteria (Tables 2 and 3), summed as points,
//! into five classes from Oncogenic to Benign.
//!
//! A sibling of fastvep-classification. The SOP is
//! adapted from ACMG/AMP but asks a different question (somatic oncogenicity)
//! The two share annotation *data* (the deserialised gnomAD, REVEL and SpliceAI records)
//! and nothing else.
//!
//! TODO: Several criteria need evidence no fastVEP source provides yet: a curated
//! gene role (OVS1, OM2), cancer hotspots (OS3, OM3, OP3), variants already
//! classified under this SOP (OS1, OM4), critical functional regions (OM1),
//! curated functional studies (OS2, SBS2) and the tumor type (OP2). These are
//! typed in [`evidence`] and report `NotEvaluable` until supplied. Until then a
//! classification sums only the questions this run could ask, and every one it
//! could not is listed in `criteria` as `NotEvaluable`, so a reader can see how
//! much of the SOP the call rests on.

pub mod combiner;
pub mod config;
pub mod criteria;
pub mod evidence;
pub mod input;
pub mod sa_extract;
pub mod types;

pub use config::OncogenicityConfig;
pub use input::{OncogenicityInput, PopulationEvidence, VariantClass};
pub use sa_extract::{extract_oncogenicity_input, AlleleTranscriptView};
pub use types::{
    CriterionCode, CriterionResult, Direction, OncogenicityClassification, OncogenicityResult,
    Outcome, Strength,
};

/// Classify one allele-transcript pair.
pub fn classify(input: &OncogenicityInput, config: &OncogenicityConfig) -> OncogenicityResult {
    let criteria = criteria::evaluate_all(input, config);
    let points = criteria.iter().map(CriterionResult::points).sum();
    let classification = combiner::classify_points(points);
    OncogenicityResult {
        classification,
        shorthand: classification.shorthand().to_string(),
        points,
        criteria,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{FunctionalEffect, HotspotEvidence, OncogenicFunctionalEvidence};
    use fastvep_classification::sa_extract::RevelData;
    use fastvep_core::Consequence;

    #[test]
    fn nothing_known_is_uncertain_not_benign() {
        // A missense with no annotation at all: every criterion is either not
        // applicable or not evaluable, and none may count against it.
        let input = OncogenicityInput {
            consequences: vec![Consequence::MissenseVariant],
            amino_acids: Some(("R".into(), "C".into())),
            ..Default::default()
        };
        let result = classify(&input, &OncogenicityConfig::default());
        assert_eq!(result.points, 0);
        assert_eq!(
            result.classification,
            OncogenicityClassification::UncertainSignificance
        );
        assert!(result.criteria.iter().all(|c| c.outcome != Outcome::Met));
    }

    #[test]
    fn braf_g469e_worked_example() {
        // Horak 2022 example 2 satisfies OS2 + OP1 + OP3 + OP4 = 7 points,
        // Likely Oncogenic. Horak used FATHMM + CADD for OP1; FastVEP does not
        // use those predictors here (TODO?), so this test supplies a REVEL score that
        // satisfies FastVEP's configured OP1 computational-evidence rule.
        let input = OncogenicityInput {
            consequences: vec![Consequence::MissenseVariant],
            amino_acids: Some(("G".into(), "E".into())),
            population: PopulationEvidence::NoRecord,
            revel: Some(RevelData { score: Some(0.95) }),
            hotspot: Some(HotspotEvidence::Hotspot {
                residue_samples: 100,
                aa_change_samples: 4,
                nucleotide_resolved: false,
            }),
            functional: Some(OncogenicFunctionalEvidence {
                effect: FunctionalEffect::Oncogenic,
                strength: Strength::Strong,
                nucleotide_specific: true,
                source: None,
            }),
            ..Default::default()
        };
        let result = classify(&input, &OncogenicityConfig::default());
        let met: Vec<&str> = result
            .criteria
            .iter()
            .filter(|c| c.outcome == Outcome::Met)
            .map(|c| c.code.as_str())
            .collect();
        assert_eq!(met, ["OS2", "OP1", "OP3", "OP4"]);
        assert_eq!(result.points, 7);
        assert_eq!(result.shorthand, "LO");
    }

    #[test]
    fn result_serialises_with_outcomes_and_codes() {
        let result = classify(&OncogenicityInput::default(), &Default::default());
        let json = serde_json::to_value(&result).unwrap();
        // TODO: This test depends on criterion ordering in the serialized array.
        // Reordering criteria would break it even if serialization remains correct.
        // Prefer looking up criteria by `code` unless SOP order is part of the API contract.
        // (Design is still WIP so not editing this yet)
        assert_eq!(json["criteria"][0]["code"], "OVS1");
        assert_eq!(json["criteria"][9]["outcome"], "not_evaluable");
        assert_eq!(json["classification"], "UncertainSignificance");
    }
}
