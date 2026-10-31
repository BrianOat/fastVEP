//! OP1, SBP1, SBP2: computational evidence.
//!
//! OP1 and SBP1 turn on "all utilized lines of computational evidence
//! (conservation/evolutionary, splicing impact, etc.)", and Table 2 adds that
//! algorithms sharing inputs must not count as independent lines. The SOP does
//! not say which lines apply to which variant, so this module decides, and the
//! decision is recorded here because it is not SOP text:
//!
//! - **Missense.** The protein line (REVEL) and the splicing line (SpliceAI).
//!   OP1 needs the protein line to predict impact; a splicing prediction of no
//!   effect does not count against it, because a missense variant need not act
//!   through splicing to be oncogenic. Both worked examples in Horak 2022 award
//!   OP1 to missense variants on missense predictors alone. SBP1 is the
//!   opposite claim, *no* impact, so it needs every mechanism ruled out: REVEL
//!   and SpliceAI must both predict none. PhyloP is not a separate line here:
//!   REVEL is an ensemble whose inputs include conservation scores so counting
//!   PhyloP beside it would be the double counting the SOP warns about.
//! - **Synonymous and other changes that alter no residue.** Splicing
//!   (SpliceAI) and conservation (PhyloP), and both must agree. With no protein
//!   effect to weigh, splicing is the mechanism and conservation is the only
//!   independent check on it.
//! - (TODO) **Predicted null and in-frame length changes.** Not evaluated. No line in
//!   the pipeline is calibrated for them: REVEL scores missense only, and
//!   SpliceAI and PhyloP say nothing about a truncated or shortened protein.
//!   Their effect is what OVS1 and OM2 assess.
//!
//! The thresholds are in [`OncogenicityConfig`]. A score between the impact
//! and no-impact thresholds supports neither criterion.

use serde_json::json;

use crate::config::OncogenicityConfig;
use crate::input::{OncogenicityInput, VariantClass};
use crate::types::{CriterionCode, CriterionResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Verdict {
    Impact,
    NoImpact,
    Indeterminate,
}

fn banded(score: f64, impact_at: f64, no_impact_at: f64) -> Verdict {
    if score >= impact_at {
        Verdict::Impact
    } else if score <= no_impact_at {
        Verdict::NoImpact
    } else {
        Verdict::Indeterminate
    }
}

/// The lines available for one variant. `None` for a line means its data is
/// missing from this run.
struct Lines {
    revel: Option<(f64, Verdict)>,
    splicing: Option<(f64, Verdict)>,
    conservation: Option<(f64, Verdict)>,
}

impl Lines {
    fn read(input: &OncogenicityInput, config: &OncogenicityConfig) -> Self {
        let revel = input
            .revel
            .as_ref()
            .and_then(|r| r.score)
            .map(|s| (s, banded(s, config.revel_impact, config.revel_no_impact)));
        let splicing = input
            .splice_ai
            .as_ref()
            .and_then(|s| s.max_delta_score())
            .map(|s| {
                (
                    s,
                    banded(s, config.spliceai_impact, config.spliceai_no_impact),
                )
            });
        // Conservation is a single cut: conserved or not. The SOP's SBP2 asks
        // only whether the nucleotide is "highly conserved".
        let conservation = input.phylop.map(|p| {
            let v = if p >= config.phylop_conserved {
                Verdict::Impact
            } else {
                Verdict::NoImpact
            };
            (p, v)
        });
        Self {
            revel,
            splicing,
            conservation,
        }
    }

    fn details(&self) -> serde_json::Value {
        json!({
            "revel": self.revel.map(|l| l.0),
            "spliceai_max_delta": self.splicing.map(|l| l.0),
            "phylop": self.conservation.map(|l| l.0),
        })
    }
}

/// Return the reason why OP1/SBP1 should not be evaluated for this variant class.
/// Return None if the class is supported.
fn unsupported_class(class: VariantClass) -> Option<&'static str> {
    match class {
        // OP1/SBP1 are not meant to judge frameshifts, stop-gains, etc
        VariantClass::PredictedNull => Some(
            "No computational line in this run is calibrated for predicted-null variants \
             (assessed under OVS1)",
        ),
        VariantClass::ProteinLengthChange => Some(
            "No computational line in this run is calibrated for in-frame length changes \
             (assessed under OM2)",
        ),
        _ => None,
    }
}

/// OP1: all utilized lines of computational evidence support an oncogenic
/// effect.
pub fn evaluate_op1(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Op1;
    let class = input.variant_class();
    if let Some(why) = unsupported_class(class) {
        return CriterionResult::not_evaluable(CODE, why);
    }
    let lines = Lines::read(input, config);
    let details = lines.details();
    let result = if class == VariantClass::Missense {
        match lines.revel {
            None => {
                CriterionResult::not_evaluable(CODE, "No REVEL score for this missense variant")
            }
            Some((s, Verdict::Impact)) => {
                CriterionResult::met(CODE, format!("REVEL {s:.3} predicts a damaging effect"))
            }
            Some((s, _)) => CriterionResult::not_met(
                CODE,
                format!("REVEL {s:.3} does not predict a damaging effect"),
            ),
        }
    } else {
        match (lines.splicing, lines.conservation) {
            (Some((s, Verdict::Impact)), Some((p, Verdict::Impact))) => CriterionResult::met(
                CODE,
                format!("SpliceAI {s:.2} predicts a splicing effect at a conserved base (PhyloP {p:.2})"),
            ),
            (Some(_), Some(_)) => CriterionResult::not_met(
                CODE,
                "Splicing and conservation lines do not both support an effect",
            ),
            _ => CriterionResult::not_evaluable(
                CODE,
                "Needs both SpliceAI and PhyloP for a variant that changes no residue",
            ),
        }
    };
    result.with_details(details)
}

/// SBP1: all utilized lines of computational evidence suggest no impact.
pub fn evaluate_sbp1(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Sbp1;
    let class = input.variant_class();
    if let Some(why) = unsupported_class(class) {
        return CriterionResult::not_evaluable(CODE, why);
    }
    let lines = Lines::read(input, config);
    let details = lines.details();
    let (first, second, needs) = if class == VariantClass::Missense {
        (lines.revel, lines.splicing, "REVEL and SpliceAI")
    } else {
        (lines.splicing, lines.conservation, "SpliceAI and PhyloP")
    };
    let result = match (first, second) {
        (Some((_, Verdict::NoImpact)), Some((_, Verdict::NoImpact))) => {
            CriterionResult::met(CODE, format!("{needs} both predict no impact"))
        }
        (Some(_), Some(_)) => {
            CriterionResult::not_met(CODE, format!("{needs} do not both predict no impact"))
        }
        _ => CriterionResult::not_evaluable(CODE, format!("Needs both {needs}")),
    };
    result.with_details(details)
}

/// SBP2: a synonymous variant with no predicted splicing impact at a base that
/// is not highly conserved.
pub fn evaluate_sbp2(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Sbp2;
    if input.variant_class() != VariantClass::Synonymous {
        return CriterionResult::not_met(CODE, "Not a synonymous variant");
    }
    let lines = Lines::read(input, config);
    let details = lines.details();
    // SpliceAI's maximum delta spans acceptor and donor, gain and loss, which
    // is SBP2's "no impact to the splice consensus sequence nor the creation
    // of a new splice site".
    let result = match (lines.splicing, lines.conservation) {
        (Some((s, Verdict::NoImpact)), Some((p, Verdict::NoImpact))) => CriterionResult::met(
            CODE,
            format!("Synonymous; SpliceAI {s:.2} predicts no splicing impact; PhyloP {p:.2} is not highly conserved"),
        ),
        (Some((s, _)), Some((p, _))) => CriterionResult::not_met(
            CODE,
            format!("Synonymous, but SpliceAI {s:.2} / PhyloP {p:.2} do not both indicate no impact"),
        ),
        _ => CriterionResult::not_evaluable(CODE, "Needs both SpliceAI and PhyloP"),
    };
    result.with_details(details)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Outcome;
    use fastvep_classification::sa_extract::{RevelData, SpliceAiData};
    use fastvep_core::Consequence;

    fn input(
        consequence: Consequence,
        revel: Option<f64>,
        spliceai: Option<f64>,
        phylop: Option<f64>,
    ) -> OncogenicityInput {
        OncogenicityInput {
            consequences: vec![consequence],
            revel: revel.map(|s| RevelData { score: Some(s) }),
            splice_ai: spliceai.map(|s| SpliceAiData {
                ds_ag: Some(s),
                ..Default::default()
            }),
            phylop,
            ..Default::default()
        }
    }

    fn outcomes(i: &OncogenicityInput) -> (Outcome, Outcome, Outcome) {
        let c = OncogenicityConfig::default();
        (
            evaluate_op1(i, &c).outcome,
            evaluate_sbp1(i, &c).outcome,
            evaluate_sbp2(i, &c).outcome,
        )
    }

    #[test]
    fn missense_op1_rests_on_the_protein_line() {
        let damaging = input(Consequence::MissenseVariant, Some(0.9), Some(0.0), None);
        assert_eq!(outcomes(&damaging).0, Outcome::Met);
        let middling = input(Consequence::MissenseVariant, Some(0.5), None, None);
        assert_eq!(outcomes(&middling).0, Outcome::NotMet);
        let no_revel = input(Consequence::MissenseVariant, None, Some(0.9), Some(5.0));
        assert_eq!(outcomes(&no_revel).0, Outcome::NotEvaluable);
    }

    #[test]
    fn missense_sbp1_needs_splicing_ruled_out_too() {
        let benign = input(Consequence::MissenseVariant, Some(0.1), Some(0.01), None);
        assert_eq!(outcomes(&benign).1, Outcome::Met);
        let splice_hit = input(Consequence::MissenseVariant, Some(0.1), Some(0.5), None);
        assert_eq!(outcomes(&splice_hit).1, Outcome::NotMet);
        let no_splice_data = input(Consequence::MissenseVariant, Some(0.1), None, None);
        assert_eq!(outcomes(&no_splice_data).1, Outcome::NotEvaluable);
    }

    #[test]
    fn op1_and_sbp1_never_both_met() {
        for revel in [0.0, 0.2, 0.5, 0.7, 1.0] {
            for sai in [0.0, 0.15, 0.5] {
                for c in [Consequence::MissenseVariant, Consequence::IntronVariant] {
                    let i = input(c, Some(revel), Some(sai), Some(sai * 10.0));
                    let (op1, sbp1, _) = outcomes(&i);
                    assert!(
                        !(op1 == Outcome::Met && sbp1 == Outcome::Met),
                        "{c:?} revel={revel} spliceai={sai}"
                    );
                }
            }
        }
    }

    #[test]
    fn non_coding_lines_must_agree() {
        let both = input(Consequence::IntronVariant, None, Some(0.8), Some(4.0));
        assert_eq!(outcomes(&both).0, Outcome::Met);
        let unconserved = input(Consequence::IntronVariant, None, Some(0.8), Some(0.5));
        assert_eq!(outcomes(&unconserved).0, Outcome::NotMet);
        let quiet = input(Consequence::IntronVariant, None, Some(0.02), Some(0.5));
        assert_eq!(outcomes(&quiet).1, Outcome::Met);
    }

    #[test]
    fn null_and_length_changes_are_not_asked() {
        for c in [
            Consequence::StopGained,
            Consequence::FrameshiftVariant,
            Consequence::InframeDeletion,
        ] {
            let i = input(c, Some(0.99), Some(0.9), Some(8.0));
            let (op1, sbp1, _) = outcomes(&i);
            assert_eq!(op1, Outcome::NotEvaluable, "{c:?}");
            assert_eq!(sbp1, Outcome::NotEvaluable, "{c:?}");
        }
    }

    #[test]
    fn sbp2_synonymous_only() {
        let syn = input(Consequence::SynonymousVariant, None, Some(0.01), Some(0.3));
        assert_eq!(outcomes(&syn).2, Outcome::Met);
        let conserved = input(Consequence::SynonymousVariant, None, Some(0.01), Some(6.0));
        assert_eq!(outcomes(&conserved).2, Outcome::NotMet);
        let no_phylop = input(Consequence::SynonymousVariant, None, Some(0.01), None);
        assert_eq!(outcomes(&no_phylop).2, Outcome::NotEvaluable);
        let intronic = input(Consequence::IntronVariant, None, Some(0.01), Some(0.3));
        assert_eq!(outcomes(&intronic).2, Outcome::NotMet);
    }

    #[test]
    fn indeterminate_spliceai_supports_neither_direction() {
        let i = input(Consequence::SynonymousVariant, None, Some(0.15), Some(0.3));
        let (op1, sbp1, sbp2) = outcomes(&i);
        assert_eq!(op1, Outcome::NotMet);
        assert_eq!(sbp1, Outcome::NotMet);
        assert_eq!(sbp2, Outcome::NotMet);
    }
}
