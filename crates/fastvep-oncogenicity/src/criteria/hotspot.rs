//! OS3, OM3, OP3: cancerhotspots.org recurrence.
//!
//! The three are one piece of evidence at three strengths, split by two counts
//! Table 2 states exactly: samples with any somatic variant at the residue
//! (50) and samples with this amino-acid change (10). The bands are disjoint by
//! construction, so at most one of them is ever met.
//!
//! Horak says that for a tumor type not well covered by CancerHotspots, 
//! resources such as COSMIC or a tumor-type-specific study could be used
//! 
//! TODO: Address COSMIC/Tumor Type Implementation

use serde_json::json;

use crate::evidence::HotspotEvidence;
use crate::input::{OncogenicityInput, VariantClass};
use crate::types::{CriterionCode, CriterionResult};

const RESIDUE_SAMPLES_STRONG: u32 = 50;
const AA_CHANGE_SAMPLES: u32 = 10;

/// Evaluate OS3, OM3 and OP3 together, in that order.
/// Tier should only be one of the 3, mutually exclusive conditions
pub fn evaluate_hotspot(input: &OncogenicityInput) -> [CriterionResult; 3] {
    use CriterionCode::{Om3, Op3, Os3};
    let (residue, change) = match &input.hotspot {
        None => {
            let why = "No cancer hotspot source is available";
            return [Os3, Om3, Op3].map(|c| CriterionResult::not_evaluable(c, why));
        }
        Some(HotspotEvidence::NotHotspot) => {
            let why = "Residue is not a cancerhotspots.org hotspot";
            return [Os3, Om3, Op3].map(|c| CriterionResult::not_met(c, why));
        }
        Some(HotspotEvidence::Hotspot {
            residue_samples,
            aa_change_samples,
            ..
        }) => (*residue_samples, *aa_change_samples),
    };

    let details = json!({
        "residue_samples": residue,
        "aa_change_samples": change,
        // Table 2 cautions about hotspots driven by truncating variants; a
        // predicted-null variant at a hotspot is where that applies.
        "truncating_caution": input.variant_class() == VariantClass::PredictedNull,
    });
    let tier = if change >= AA_CHANGE_SAMPLES && residue >= RESIDUE_SAMPLES_STRONG {
        Os3
    } else if change >= AA_CHANGE_SAMPLES {
        Om3
    } else {
        // Includes zero: Table 2 asks only that the change count be "below 10".
        Op3
    };
    let summary = format!(
        "Hotspot residue with {residue} samples; this amino-acid change in {change} samples"
    );
    [Os3, Om3, Op3].map(|c| {
        let r = if c == tier {
            CriterionResult::met(c, summary.as_str())
        } else {
            CriterionResult::not_met(c, format!("{summary}; meets {}", tier.as_str()))
        };
        r.with_details(details.clone())
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Outcome;

    fn at(residue_samples: u32, aa_change_samples: u32) -> [Outcome; 3] {
        let i = OncogenicityInput {
            hotspot: Some(HotspotEvidence::Hotspot {
                residue_samples,
                aa_change_samples,
                nucleotide_resolved: false,
            }),
            ..Default::default()
        };
        evaluate_hotspot(&i).map(|r| r.outcome)
    }

    const MET_NOT_NOT: [Outcome; 3] = [Outcome::Met, Outcome::NotMet, Outcome::NotMet];
    const NOT_MET_NOT: [Outcome; 3] = [Outcome::NotMet, Outcome::Met, Outcome::NotMet];
    const NOT_NOT_MET: [Outcome; 3] = [Outcome::NotMet, Outcome::NotMet, Outcome::Met];

    #[test]
    fn table_2_boundaries() {
        assert_eq!(at(50, 10), MET_NOT_NOT);
        assert_eq!(at(49, 10), NOT_MET_NOT);
        assert_eq!(at(500, 9), NOT_NOT_MET);
        assert_eq!(at(10, 10), NOT_MET_NOT);
        // BRAF p.Gly469Glu, Horak 2022 example 2: change count 4 gives OP3.
        assert_eq!(at(200, 4), NOT_NOT_MET);
    }

    #[test]
    fn without_hotspot_source_criteria_are_not_evaluable() {
        let r = evaluate_hotspot(&OncogenicityInput::default());

        assert!(r.iter().all(|c| c.outcome == Outcome::NotEvaluable));
    }

    #[test]
    fn non_hotspot_residue_does_not_meet_hotspot_criteria() {
        let input = OncogenicityInput {
            hotspot: Some(HotspotEvidence::NotHotspot),
            ..Default::default()
        };

        let r = evaluate_hotspot(&input);

        assert!(r.iter().all(|c| c.outcome == Outcome::NotMet));
    }
}
