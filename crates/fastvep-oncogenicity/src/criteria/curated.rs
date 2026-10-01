//! OM1, OS2, SBS2, OP2: criteria that read curated evidence or clinical
//! context, and nothing the pipeline computes.

use serde_json::json;

use crate::evidence::{FunctionalDomainEvidence, FunctionalEffect};
use crate::input::OncogenicityInput;
use crate::types::{CriterionCode, CriterionResult};

/// OM1: located in a critical and well-established part of a functional
/// domain. Its exclusion by OS1 and OS3 is applied in
/// [`super::apply_exclusions`].
pub fn evaluate_om1(input: &OncogenicityInput) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Om1;
    match &input.functional_domain {
        None => CriterionResult::not_evaluable(
            CODE,
            "No curated map of critical functional regions is available",
        ),
        Some(FunctionalDomainEvidence::NotInCriticalRegion) => {
            CriterionResult::not_met(CODE, "Outside every curated critical functional region")
        }
        Some(FunctionalDomainEvidence::InCriticalRegion { region }) => {
            CriterionResult::met(CODE, format!("In a critical functional region: {region}"))
        }
    }
}

/// OS2 (studies support an oncogenic effect) and SBS2 (studies show none),
/// evaluated together because one curated result answers both.
///
/// The curator's strength is applied as given when it is a downgrade, which the
/// SOP permits for both criteria. It is capped at Strong: the SOP permits no
/// upgrade.
pub fn evaluate_functional(input: &OncogenicityInput) -> [CriterionResult; 2] {
    use CriterionCode::{Os2, Sbs2};
    let Some(f) = &input.functional else {
        let why = "No curated functional-study result for this variant";
        return [Os2, Sbs2].map(|c| CriterionResult::not_evaluable(c, why));
    };
    let (met, other) = match f.effect {
        FunctionalEffect::Oncogenic => (Os2, Sbs2),
        FunctionalEffect::NoOncogenicEffect => (Sbs2, Os2),
    };
    // Downgrade as per curator
    let strength = f.strength.min(met.default_strength());
    let source = f.source.as_deref().unwrap_or("no source given");
    let details = json!({
        "strength_assigned": f.strength.as_str(),
        "nucleotide_specific": f.nucleotide_specific,
        "source": f.source,
    });
    let met_result = CriterionResult::met(
        met,
        format!(
            "Curated functional studies ({source}) at {}",
            strength.as_str()
        ),
    )
    .with_strength(strength)
    .with_details(details.clone());
    let other_result = CriterionResult::not_met(
        other,
        format!("Curated functional studies ({source}) point the other way"),
    )
    .with_details(details);
    match met {
        Os2 => [met_result, other_result],
        _ => [other_result, met_result],
    }
}

/// OP2: somatic variant in a gene in a malignancy with a single genetic
/// etiology.
///
/// A question about the tumor, not the variant, and the pipeline carries no
/// tumor type.
/// TODO!!!
pub fn evaluate_op2(_input: &OncogenicityInput) -> CriterionResult {
    CriterionResult::not_evaluable(
        CriterionCode::Op2,
        "Requires the tumor type, which this run does not carry",
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::OncogenicFunctionalEvidence;
    use crate::types::{Outcome, Strength};

    fn functional(effect: FunctionalEffect, strength: Strength) -> OncogenicityInput {
        OncogenicityInput {
            functional: Some(OncogenicFunctionalEvidence {
                effect,
                strength,
                nucleotide_specific: true,
                source: Some("PMID:1".into()),
            }),
            ..Default::default()
        }
    }

    #[test]
    fn om1_three_states() {
        assert_eq!(
            evaluate_om1(&OncogenicityInput::default()).outcome,
            Outcome::NotEvaluable
        );

        // Source consulted; residue is not in a critical region.
        let outside = OncogenicityInput {
            functional_domain: Some(FunctionalDomainEvidence::NotInCriticalRegion),
            ..Default::default()
        };
        assert_eq!(evaluate_om1(&outside).outcome, Outcome::NotMet);

        // Source consulted; residue is in a critical region.
        let inside = OncogenicityInput {
            functional_domain: Some(FunctionalDomainEvidence::InCriticalRegion {
                region: "active site".into(),
            }),
            ..Default::default()
        };
        assert_eq!(evaluate_om1(&inside).points(), 2);
    }

    #[test]
    fn functional_evidence_answers_one_criterion_and_refuses_the_other() {
        let [os2, sbs2] =
            evaluate_functional(&functional(FunctionalEffect::Oncogenic, Strength::Strong));
        assert_eq!((os2.outcome, sbs2.outcome), (Outcome::Met, Outcome::NotMet));
        assert_eq!(os2.points(), 4);
        let [os2, sbs2] = evaluate_functional(&functional(
            FunctionalEffect::NoOncogenicEffect,
            Strength::Strong,
        ));
        assert_eq!((os2.outcome, sbs2.outcome), (Outcome::NotMet, Outcome::Met));
        assert_eq!(sbs2.points(), -4);
    }

    #[test]
    fn functional_strength_can_be_downgraded_not_upgraded() {
        // Horak 2022 example 1: SBS2 downgraded to Supporting for KRAS R164Q.
        let [_, sbs2] = evaluate_functional(&functional(
            FunctionalEffect::NoOncogenicEffect,
            Strength::Supporting,
        ));
        assert_eq!(sbs2.points(), -1);
        let [os2, _] = evaluate_functional(&functional(
            FunctionalEffect::Oncogenic,
            Strength::VeryStrong,
        ));
        assert_eq!(os2.strength, Strength::Strong);
    }

    #[test]
    fn no_functional_data_is_not_negative_evidence() {
        let r = evaluate_functional(&OncogenicityInput::default());
        assert!(r.iter().all(|c| c.outcome == Outcome::NotEvaluable));
        assert!(r.iter().all(|c| c.points() == 0));
    }

    #[test]
    fn op2_is_not_evaluable() {
        assert_eq!(
            evaluate_op2(&OncogenicityInput::default()).outcome,
            Outcome::NotEvaluable
        );
    }
}
