//! The seventeen criteria, and the exclusions between them.

mod computational;
mod curated;
mod established;
mod gene_role;
mod hotspot;
mod population;

use crate::config::OncogenicityConfig;
use crate::evidence::HotspotEvidence;
use crate::input::OncogenicityInput;
use crate::types::{CriterionCode, CriterionResult, Outcome};

use CriterionCode::*;

/// Evaluate every criterion, then apply the exclusions. Returned in SOP table
/// order.
pub fn evaluate_all(
    input: &OncogenicityInput,
    config: &OncogenicityConfig,
) -> Vec<CriterionResult> {
    let [os3, om3, op3] = hotspot::evaluate_hotspot(input);
    let [os2, sbs2] = curated::evaluate_functional(input);
    let mut criteria = vec![
        gene_role::evaluate_ovs1(input),
        established::evaluate_os1(input, config),
        os2,
        os3,
        curated::evaluate_om1(input),
        gene_role::evaluate_om2(input),
        om3,
        established::evaluate_om4(input),
        computational::evaluate_op1(input, config),
        curated::evaluate_op2(input),
        op3,
        population::evaluate_op4(input, config),
        population::evaluate_sbvs1(input, config),
        population::evaluate_sbs1(input, config),
        sbs2,
        computational::evaluate_sbp1(input, config),
        computational::evaluate_sbp2(input, config),
    ];
    apply_exclusions(&mut criteria, input);
    criteria
}

/// `(criterion, criteria whose applicability normally excludes it)`.
///
/// OS2 and OS3 have SOP-defined exceptions to their OS1 exclusion; those are
/// handled by `os1_exclusion_waived`.
///
/// FastVEP currently interprets "applicable" as "the criterion's underlying
/// conditions are met", even if that criterion is subsequently excluded.
/// Exclusions are therefore resolved against the raw pre-exclusion outcomes.
///
/// TODO: Validate this interpretation of "applicable" against the intended SOP
/// semantics before treating it as final.
/// The final two entries are FastVEP anti-double-counting policies, not
/// explicit Table 3 exclusions.
const EXCLUSIONS: &[(CriterionCode, &[CriterionCode])] = &[
    // "If OS1 is applicable, this rule can be used only if functional studies
    // are based on the particular nucleotide change of the variant."
    (Os2, &[Os1]),
    // "This rule cannot be used if OS1 is applicable, unless it is possible to
    // observe hotspots based on the particular nucleotide change."
    (Os3, &[Os1]),
    // "This rule cannot be used if OS1 or OS3 is applicable."
    (Om1, &[Os1, Os3]),
    // "This rule cannot be used if OVS1 is applicable."
    (Om2, &[Ovs1]),
    // "This rule cannot be used if OS1 or OS3 or OM1 is applicable."
    (Om4, &[Os1, Os3, Om1]),
    // "This rule cannot be used if OM1 or OM4 is applicable."
    (Om3, &[Om1, Om4]),
    // Not in Table 3. SBVS1 and SBS1 are one frequency measured against two
    // thresholds, so anything over 5% is also over 1%. Counting both would give
    // a common variant -12 for a single observation.
    (Sbs1, &[Sbvs1]),
    // Not in Table 3. For a synonymous variant SBP2 and SBP1 read the same two
    // predictions (SpliceAI, PhyloP), and the SOP says of SBP1 that algorithms
    // using the same inputs must not count as independent evidence.
    (Sbp1, &[Sbp2]),
];

/// Whether Table 2's escape clause lifts the OS1 exclusion for `code`.
fn os1_exclusion_waived(code: CriterionCode, input: &OncogenicityInput) -> bool {
    match code {
        Os2 => input
            .functional
            .as_ref()
            .is_some_and(|f| f.nucleotide_specific),
        Os3 => matches!(
            input.hotspot,
            Some(HotspotEvidence::Hotspot {
                nucleotide_resolved: true,
                ..
            })
        ),
        _ => false,
    }
}

pub(crate) fn apply_exclusions(criteria: &mut [CriterionResult], input: &OncogenicityInput) {
    // Decided on the raw outcomes, before any exclusion, per the reading of
    // "applicable" above.
    let met: Vec<CriterionCode> = criteria
        .iter()
        .filter(|c| c.outcome == Outcome::Met)
        .map(|c| c.code)
        .collect();
    for (target, excluders) in EXCLUSIONS {
        let Some(c) = criteria
            .iter_mut()
            .find(|c| c.code == *target && c.outcome == Outcome::Met)
        else {
            continue;
        };
        let by: Vec<CriterionCode> = excluders
            .iter()
            .copied()
            .filter(|e| met.contains(e))
            .filter(|e| !(*e == Os1 && os1_exclusion_waived(*target, input)))
            .collect();
        if !by.is_empty() {
            c.outcome = Outcome::Excluded;
            c.excluded_by = by;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{
        EstablishedOncogenicVariants, EstablishedVariant, FunctionalDomainEvidence,
        FunctionalEffect, GeneRole, OncogenicFunctionalEvidence,
    };
    use crate::input::PopulationEvidence;
    use crate::types::Strength;
    use fastvep_classification::sa_extract::GnomadData;
    use fastvep_core::Consequence;

    fn outcome(criteria: &[CriterionResult], code: CriterionCode) -> &CriterionResult {
        criteria.iter().find(|c| c.code == code).unwrap()
    }

    /// A missense with an established identical change (OS1), in a critical
    /// domain (OM1), at a hotspot. Adjust per test.
    fn missense_with_os1() -> OncogenicityInput {
        OncogenicityInput {
            consequences: vec![Consequence::MissenseVariant],
            amino_acids: Some(("G".into(), "D".into())),
            established_oncogenic: Some(EstablishedOncogenicVariants {
                at_residue: vec![EstablishedVariant {
                    ref_aa: "G".into(),
                    alt_aa: "D".into(),
                    is_self: false,
                }],
            }),
            ..Default::default()
        }
    }

    fn hotspot(residue: u32, change: u32, nucleotide_resolved: bool) -> Option<HotspotEvidence> {
        Some(HotspotEvidence::Hotspot {
            residue_samples: residue,
            aa_change_samples: change,
            nucleotide_resolved,
        })
    }

    #[test]
    fn all_seventeen_in_table_order() {
        let criteria = evaluate_all(&OncogenicityInput::default(), &Default::default());
        let codes: Vec<&str> = criteria.iter().map(|c| c.code.as_str()).collect();
        assert_eq!(
            codes,
            [
                "OVS1", "OS1", "OS2", "OS3", "OM1", "OM2", "OM3", "OM4", "OP1", "OP2", "OP3",
                "OP4", "SBVS1", "SBS1", "SBS2", "SBP1", "SBP2"
            ]
        );
    }

    #[test]
    fn os1_excludes_os3_unless_nucleotide_resolved() {
        let mut i = missense_with_os1();
        i.hotspot = hotspot(100, 20, false);
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Os3).outcome, Outcome::Excluded);
        assert_eq!(outcome(&c, Os3).excluded_by, [Os1]);
        assert_eq!(outcome(&c, Os3).points(), 0);

        i.hotspot = hotspot(100, 20, true);
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Os3).outcome, Outcome::Met);
    }

    #[test]
    fn os1_excludes_os2_unless_studies_used_this_nucleotide_change() {
        let mut i = missense_with_os1();
        i.functional = Some(OncogenicFunctionalEvidence {
            effect: FunctionalEffect::Oncogenic,
            strength: Strength::Strong,
            nucleotide_specific: false,
            source: None,
        });
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Os2).outcome, Outcome::Excluded);

        i.functional.as_mut().unwrap().nucleotide_specific = true;
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Os2).outcome, Outcome::Met);
    }

    #[test]
    fn om1_excluded_by_os3_and_om3_by_om1() {
        let mut i = missense_with_os1();
        i.established_oncogenic = Some(Default::default());
        i.functional_domain = Some(FunctionalDomainEvidence::InCriticalRegion {
            region: "P-loop".into(),
        });

        // Strong hotspot: OS3 counts, OM1 does not.
        i.hotspot = hotspot(100, 20, false);
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Os3).outcome, Outcome::Met);
        assert_eq!(outcome(&c, Om1).excluded_by, [Os3]);

        // Moderate hotspot: OM1 counts, OM3 does not.
        i.hotspot = hotspot(30, 20, false);
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Om1).outcome, Outcome::Met);
        assert_eq!(outcome(&c, Om3).excluded_by, [Om1]);
    }

    #[test]
    fn exclusion_reads_raw_applicability() {
        // OS1 rules out OM1; OM1's conditions still hold, so OM3 stays out.
        let mut i = missense_with_os1();
        i.functional_domain = Some(FunctionalDomainEvidence::InCriticalRegion {
            region: "P-loop".into(),
        });
        i.hotspot = hotspot(30, 20, false);
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Om1).excluded_by, [Os1]);
        assert_eq!(outcome(&c, Om3).excluded_by, [Om1]);
        let total: i32 = c.iter().map(|c| c.points()).sum();
        assert_eq!(total, 4, "OS1 alone");
    }

    #[test]
    fn ovs1_excludes_om2() {
        // A frameshift and an in-frame term on the same allele-transcript.
        let i = OncogenicityInput {
            consequences: vec![Consequence::FrameshiftVariant, Consequence::StopLost],
            gene_role: Some(GeneRole {
                tumor_suppressor: true,
                oncogene: false,
            }),
            escapes_nmd: Some(false),
            ..Default::default()
        };
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Ovs1).outcome, Outcome::Met);
        assert_eq!(outcome(&c, Om2).excluded_by, [Ovs1]);
    }

    #[test]
    fn sbvs1_subsumes_sbs1() {
        let i = OncogenicityInput {
            population: PopulationEvidence::Record(GnomadData {
                afr_af: Some(0.3),
                ..Default::default()
            }),
            ..Default::default()
        };
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Sbvs1).outcome, Outcome::Met);
        assert_eq!(outcome(&c, Sbs1).excluded_by, [Sbvs1]);
        assert_eq!(c.iter().map(|c| c.points()).sum::<i32>(), -8);
    }

    #[test]
    fn sbp2_subsumes_sbp1_for_synonymous() {
        use fastvep_classification::sa_extract::SpliceAiData;
        let i = OncogenicityInput {
            consequences: vec![Consequence::SynonymousVariant],
            splice_ai: Some(SpliceAiData {
                ds_ag: Some(0.0),
                ..Default::default()
            }),
            phylop: Some(0.1),
            ..Default::default()
        };
        let c = evaluate_all(&i, &Default::default());
        assert_eq!(outcome(&c, Sbp2).outcome, Outcome::Met);
        assert_eq!(outcome(&c, Sbp1).excluded_by, [Sbp2]);
    }
}
