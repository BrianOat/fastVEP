//! OS1 and OM4: comparison with variants already classified oncogenic under
//! this SOP.
//!
//! Both require the comparison variant to be oncogenic using this standard

use serde_json::json;

use crate::config::OncogenicityConfig;
use crate::input::{OncogenicityInput, VariantClass};
use crate::types::{CriterionCode, CriterionResult};

/// OS1: same amino-acid change as a previously established oncogenic variant,
/// regardless of nucleotide change.
pub fn evaluate_os1(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Os1;
    let Some((ref_aa, alt_aa)) = input.single_residue_change() else {
        return match input.variant_class() {
            VariantClass::Missense | VariantClass::PredictedNull if input.amino_acids.is_none() => {
                CriterionResult::not_evaluable(CODE, "No amino-acid change available")
            }
            _ => CriterionResult::not_met(CODE, "Not a single amino-acid substitution"),
        };
    };
    let Some(established) = &input.established_oncogenic else {
        return CriterionResult::not_evaluable(
            CODE,
            "No source of SOP-classified oncogenic variants is available",
        );
    };
    let same_change = established
        .at_residue
        .iter()
        .any(|v| !v.is_self && v.ref_aa == ref_aa && v.alt_aa == alt_aa);
    if !same_change {
        return CriterionResult::not_met(
            CODE,
            format!("No other variant producing {ref_aa}>{alt_aa} is established as oncogenic"),
        );
    }
    // Table 2: "Beware of changes that impact splicing rather than at the
    // amino acid/protein level." Surfaced rather than acted on; whether the
    // splice effect or the amino-acid change drives the phenotype is the
    // curator's call.
    let splice_score = input.splice_ai.as_ref().and_then(|s| s.max_delta_score());
    let splicing_caution = splice_score.is_some_and(|s| s >= config.spliceai_impact);
    CriterionResult::met(
        CODE,
        format!("Same amino-acid change ({ref_aa}>{alt_aa}) as an established oncogenic variant"),
    )
    .with_details(json!({
        "spliceai_max_delta": splice_score,
        "splicing_caution": splicing_caution,
    }))
}

/// OM4: missense variant at a residue where a different missense variant is
/// established as oncogenic, with an amino-acid difference at least about as
/// large.
///
/// TODO: The second half of that definition is not implemented. The SOP names
/// Grantham's distance, Epstein's coefficient and Miyata's distance as the
/// metrics, and no verified table of any of them is in this crate yet; a
/// hand-typed 20x20 matrix is exactly the kind of plausible, unchecked number
/// this tool must not ship. So a residue that does carry a different oncogenic
/// missense is reported `NotEvaluable`, with the candidates in `details`, and
/// OM4 is never awarded.
pub fn evaluate_om4(input: &OncogenicityInput) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Om4;
    if input.variant_class() != VariantClass::Missense {
        return CriterionResult::not_met(CODE, "Not a missense variant");
    }
    let Some((ref_aa, alt_aa)) = input.single_residue_change() else {
        return CriterionResult::not_evaluable(CODE, "No amino-acid change available");
    };
    let Some(established) = &input.established_oncogenic else {
        return CriterionResult::not_evaluable(
            CODE,
            "No source of SOP-classified oncogenic variants is available",
        );
    };
    let others: Vec<&str> = established
        .at_residue
        .iter()
        .filter(|v| v.ref_aa == ref_aa && v.alt_aa != alt_aa && v.alt_aa != "*")
        .map(|v| v.alt_aa.as_str())
        .collect();
    if others.is_empty() {
        return CriterionResult::not_met(
            CODE,
            "No different missense change at this residue is established as oncogenic",
        );
    }
    CriterionResult::not_evaluable(
        CODE,
        "A different oncogenic missense exists at this residue, but the amino-acid difference \
         comparison the SOP requires is not implemented",
    )
    .with_details(json!({ "ref_aa": ref_aa, "alt_aa": alt_aa, "established_alt_aas": others }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{EstablishedOncogenicVariants, EstablishedVariant};
    use crate::types::Outcome;
    use fastvep_core::Consequence;

    fn entry(r: &str, a: &str, is_self: bool) -> EstablishedVariant {
        EstablishedVariant {
            ref_aa: r.into(),
            alt_aa: a.into(),
            is_self,
        }
    }

    fn missense(r: &str, a: &str, db: Option<Vec<EstablishedVariant>>) -> OncogenicityInput {
        OncogenicityInput {
            consequences: vec![Consequence::MissenseVariant],
            amino_acids: Some((r.into(), a.into())),
            established_oncogenic: db.map(|at_residue| EstablishedOncogenicVariants { at_residue }),
            ..Default::default()
        }
    }

    #[test]
    fn os1_same_amino_acid_change() {
        let cfg = OncogenicityConfig::default();
        let i = missense("V", "L", Some(vec![entry("V", "L", false)]));
        assert_eq!(evaluate_os1(&i, &cfg).outcome, Outcome::Met);
        assert_eq!(evaluate_os1(&i, &cfg).points(), 4);
        let other = missense("V", "L", Some(vec![entry("V", "M", false)]));
        assert_eq!(evaluate_os1(&other, &cfg).outcome, Outcome::NotMet);
        let no_db = missense("V", "L", None);
        assert_eq!(evaluate_os1(&no_db, &cfg).outcome, Outcome::NotEvaluable);
    }

    #[test]
    fn os1_is_not_evidence_for_itself() {
        let i = missense("V", "L", Some(vec![entry("V", "L", true)]));
        assert_eq!(
            evaluate_os1(&i, &OncogenicityConfig::default()).outcome,
            Outcome::NotMet
        );
    }

    #[test]
    fn os1_not_asked_of_non_substitutions() {
        let i = OncogenicityInput {
            consequences: vec![Consequence::SynonymousVariant],
            amino_acids: Some(("R".into(), "R".into())),
            established_oncogenic: Some(Default::default()),
            ..Default::default()
        };
        assert_eq!(
            evaluate_os1(&i, &OncogenicityConfig::default()).outcome,
            Outcome::NotMet
        );
    }

    #[test]
    fn om4_is_never_awarded_without_a_distance_metric() {
        // The SOP's own example: p.Arg156His oncogenic, p.Arg156Cys observed.
        let i = missense("R", "C", Some(vec![entry("R", "H", false)]));
        let r = evaluate_om4(&i);
        assert_eq!(r.outcome, Outcome::NotEvaluable);
        assert_eq!(r.details["established_alt_aas"][0], "H");
        let none_other = missense("R", "C", Some(vec![entry("R", "C", false)]));
        assert_eq!(evaluate_om4(&none_other).outcome, Outcome::NotMet);
    }
}
