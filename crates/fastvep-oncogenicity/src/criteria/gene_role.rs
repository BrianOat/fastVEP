//! OVS1 and OM2: criteria that depend on the gene's role in cancer.
//!
//! Both need to know whether the gene is a tumor suppressor or an oncogene.
//! SOP asks for a "bona fide" tumor suppressor, i.e. a curated one.
//! TODO NOTE: Until a curated role reaches [`OncogenicityInput::gene_role`],
//! both criteria report `NotEvaluable` for the variants they would apply to.

use fastvep_core::Consequence;
use serde_json::json;

use crate::input::{OncogenicityInput, VariantClass};
use crate::types::{CriterionCode, CriterionResult};

/// OVS1: null variant in a bona fide tumor suppressor gene.
///
/// Table 2 attaches four cautions: a pLOF at the extreme 3' end past the NMD
/// boundary, a splice variant whose exon skip is in frame, a splice variant
/// that yields a known functional isoform, and multiple transcripts.
///
/// TODO: "Use caution" is a direction to a curator. Currently cannot exercise it, so where a
/// caution is known to apply OVS1 is left for review (`NotEvaluable`) rather
/// than awarded or refused. Where the tool cannot tell whether a caution
/// applies, OVS1 is awarded and the unassessed caution is listed in `details`.
/// TODO: Implement downgrading (as well as for OM2)
pub fn evaluate_ovs1(input: &OncogenicityInput) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Ovs1;
    if input.variant_class() != VariantClass::PredictedNull {
        return CriterionResult::not_met(CODE, "Not a null variant");
    }
    let Some(role) = input.gene_role else {
        return CriterionResult::not_evaluable(
            CODE,
            "Null variant, but no curated gene role is available to say whether the gene is \
             a tumor suppressor",
        );
    };
    if !role.tumor_suppressor {
        return CriterionResult::not_met(CODE, "Gene is not curated as a tumor suppressor");
    }

    let has = |c: Consequence| input.consequences.contains(&c);
    let truncating = has(Consequence::StopGained) || has(Consequence::FrameshiftVariant);
    let splice = has(Consequence::SpliceDonorVariant) || has(Consequence::SpliceAcceptorVariant);

    if truncating && input.escapes_nmd == Some(true) {
        return CriterionResult::not_evaluable(
            CODE,
            "Premature stop predicted to escape nonsense-mediated decay; Table 2 asks for \
             caution with such variants, which needs curator review",
        )
        .with_details(json!({ "escapes_nmd": true }));
    }

    let mut unassessed: Vec<&str> = Vec::new();
    if truncating && input.escapes_nmd.is_none() {
        unassessed.push("NMD escape not measured for this transcript");
    }
    if splice {
        unassessed.push("in-frame exon skipping and functional alternative isoforms not assessed");
    }
    if !input.is_canonical {
        // TODO: is_canonical is not sufficient to represent the
        // "multiple transcripts" OVS1 caveat. Need transcript-level
        // consequence context / alternative transcript information.
        unassessed.push("evaluated on a non-canonical transcript");
    }
    CriterionResult::met(CODE, "Null variant in a curated tumor suppressor gene")
        .with_details(json!({ "unassessed_cautions": unassessed }))
}

/// OM2: protein length change from an in-frame insertion or deletion in a
/// known oncogene or tumor suppressor, or a stop-loss in a known tumor
/// suppressor. Its exclusion by OVS1 is applied in [`super::apply_exclusions`].
pub fn evaluate_om2(input: &OncogenicityInput) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Om2;
    let has = |c: Consequence| input.consequences.contains(&c);
    let in_frame = has(Consequence::InframeInsertion) || has(Consequence::InframeDeletion);
    let stop_lost = has(Consequence::StopLost);
    if !in_frame && !stop_lost {
        if has(Consequence::ProteinAlteringVariant) {
            return CriterionResult::not_evaluable(
                CODE,
                "protein_altering_variant: whether protein length changes is not determined",
            );
        }
        return CriterionResult::not_met(CODE, "Not an in-frame indel or stop-loss");
    }
    let Some(role) = input.gene_role else {
        return CriterionResult::not_evaluable(
            CODE,
            "Protein length change, but no curated gene role is available",
        );
    };
    if in_frame && (role.oncogene || role.tumor_suppressor) {
        CriterionResult::met(
            CODE,
            "In-frame indel in a known oncogene or tumor suppressor",
        )
    } else if stop_lost && role.tumor_suppressor {
        CriterionResult::met(CODE, "Stop-loss in a known tumor suppressor")
    } else if stop_lost && !in_frame {
        CriterionResult::not_met(CODE, "Stop-loss counts only in a tumor suppressor gene")
    } else {
        CriterionResult::not_met(CODE, "Gene is not a known oncogene or tumor suppressor")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::GeneRole;
    use crate::types::Outcome;

    const TSG: GeneRole = GeneRole {
        tumor_suppressor: true,
        oncogene: false,
    };
    const ONCOGENE: GeneRole = GeneRole {
        tumor_suppressor: false,
        oncogene: true,
    };

    fn input(c: Consequence, role: Option<GeneRole>) -> OncogenicityInput {
        OncogenicityInput {
            consequences: vec![c],
            gene_role: role,
            is_canonical: true,
            escapes_nmd: Some(false),
            ..Default::default()
        }
    }

    #[test]
    fn ovs1_needs_a_curated_tumor_suppressor() {
        assert_eq!(
            evaluate_ovs1(&input(Consequence::StopGained, None)).outcome,
            Outcome::NotEvaluable
        );
        assert_eq!(
            evaluate_ovs1(&input(Consequence::StopGained, Some(ONCOGENE))).outcome,
            Outcome::NotMet
        );
        let r = evaluate_ovs1(&input(Consequence::StopGained, Some(TSG)));
        assert_eq!(r.outcome, Outcome::Met);
        assert_eq!(r.points(), 8);
    }

    #[test]
    fn ovs1_null_classes() {
        for c in [
            Consequence::StopGained,
            Consequence::FrameshiftVariant,
            Consequence::SpliceDonorVariant,
            Consequence::SpliceAcceptorVariant,
            Consequence::StartLost,
        ] {
            assert_eq!(
                evaluate_ovs1(&input(c, Some(TSG))).outcome,
                Outcome::Met,
                "{c:?}"
            );
        }
        for c in [
            Consequence::MissenseVariant,
            Consequence::SpliceRegionVariant,
            Consequence::InframeDeletion,
        ] {
            assert_eq!(
                evaluate_ovs1(&input(c, Some(TSG))).outcome,
                Outcome::NotMet,
                "{c:?}"
            );
        }
    }

    #[test]
    fn ovs1_nmd_escape_goes_to_review_not_to_a_verdict() {
        let mut i = input(Consequence::StopGained, Some(TSG));
        i.escapes_nmd = Some(true);
        assert_eq!(evaluate_ovs1(&i).outcome, Outcome::NotEvaluable);
    }

    #[test]
    fn ovs1_lists_cautions_it_could_not_assess() {
        let mut i = input(Consequence::SpliceDonorVariant, Some(TSG));
        i.is_canonical = false;
        let r = evaluate_ovs1(&i);
        assert_eq!(r.outcome, Outcome::Met);
        let cautions = r.details["unassessed_cautions"].as_array().unwrap();
        assert_eq!(cautions.len(), 2);
    }

    #[test]
    fn om2_in_frame_in_either_role_stop_loss_in_tsg_only() {
        for role in [TSG, ONCOGENE] {
            assert_eq!(
                evaluate_om2(&input(Consequence::InframeDeletion, Some(role))).outcome,
                Outcome::Met
            );
        }
        assert_eq!(
            evaluate_om2(&input(Consequence::StopLost, Some(TSG))).outcome,
            Outcome::Met
        );
        assert_eq!(
            evaluate_om2(&input(Consequence::StopLost, Some(ONCOGENE))).outcome,
            Outcome::NotMet
        );
        assert_eq!(
            evaluate_om2(&input(Consequence::InframeInsertion, None)).outcome,
            Outcome::NotEvaluable
        );
        assert_eq!(
            evaluate_om2(&input(Consequence::MissenseVariant, None)).outcome,
            Outcome::NotMet
        );
    }
}
