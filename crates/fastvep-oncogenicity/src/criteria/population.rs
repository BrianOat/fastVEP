//! OP4, SBVS1, SBS1: population frequency in gnomAD.

use fastvep_classification::sa_extract::GnomadData;
use serde_json::json;

use crate::config::OncogenicityConfig;
use crate::input::{OncogenicityInput, PopulationEvidence};
use crate::types::{CriterionCode, CriterionResult};

type PopulationAf = fn(&GnomadData) -> Option<f64>;

/// Table 3 restricts SBVS1 and SBS1 to "any of 5 general continental
/// populations: African, East Asian, European (Non-Finnish), Latino, and South
/// Asian". The SOP text adds that Finnish and Ashkenazi Jewish frequencies may
/// be used "only if founder effects have been ruled out for the gene", which is
/// a per-gene judgement this tool cannot make, so they are left out, along
/// with every group the SOP does not name.
const CONTINENTAL_POPULATIONS: [(&str, PopulationAf); 5] = [
    ("afr", |g| g.afr_af),
    ("eas", |g| g.eas_af),
    ("nfe", |g| g.nfe_af),
    ("amr", |g| g.amr_af),
    ("sas", |g| g.sas_af),
];

/// Highest allele frequency across the five continental populations, and which
/// population it came from.
fn continental_max(g: &GnomadData) -> Option<(&'static str, f64)> {
    CONTINENTAL_POPULATIONS
        .iter()
        .filter_map(|(name, af)| af(g).map(|v| (*name, v)))
        .max_by(|a, b| a.1.total_cmp(&b.1))
}

/// OP4: absent from controls, or at an extremely low frequency, in gnomAD.
pub fn evaluate_op4(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    const CODE: CriterionCode = CriterionCode::Op4;
    let g = match &input.population {
        PopulationEvidence::Unavailable => {
            return CriterionResult::not_evaluable(CODE, "No usable gnomAD lookup for this allele");
        }
        PopulationEvidence::NoRecord => {
            return CriterionResult::met(CODE, "Absent from gnomAD: no record of this allele")
                .with_details(json!({ "gnomad": "no_record" }));
        }
        PopulationEvidence::Record(g) => g,
    };
    // A non-PASS record says nothing reliable about the real frequency, in
    // either direction: the germline classifier's reasoning, and it applies
    // unchanged.
    if let Some(filter) = g.failed_filter() {
        return CriterionResult::not_evaluable(
            CODE,
            format!("gnomAD record failed the {filter} filter; frequency not usable"),
        );
    }
    if g.all_ac == Some(0) || g.all_af == Some(0.0) {
        return CriterionResult::met(CODE, "Absent from gnomAD: record carries zero alleles")
            .with_details(json!({ "gnomad_allAc": g.all_ac, "gnomad_allAf": g.all_af }));
    }
    let Some(max_af) = g.max_pop_af() else {
        return CriterionResult::not_evaluable(CODE, "gnomAD record carries no allele frequency");
    };
    let details = json!({ "gnomad_max_pop_af": max_af, "op4_max_af": config.op4_max_af });
    match config.op4_max_af {
        Some(t) if max_af <= t => CriterionResult::met(
            CODE,
            format!("gnomAD max population AF {max_af:.2e} is at or below {t:.2e}"),
        )
        .with_details(details),
        Some(t) => CriterionResult::not_met(
            CODE,
            format!("gnomAD max population AF {max_af:.2e} is above {t:.2e}"),
        )
        .with_details(details),
        None => CriterionResult::not_met(
            CODE,
            format!("Present in gnomAD (max population AF {max_af:.2e})"),
        )
        .with_details(details),
    }
}

/// SBVS1 (>5%) and SBS1 (>1%) share every step but the threshold.
fn evaluate_benign_frequency(
    code: CriterionCode,
    threshold: f64,
    input: &OncogenicityInput,
) -> CriterionResult {
    let g = match &input.population {
        PopulationEvidence::Unavailable => {
            return CriterionResult::not_evaluable(code, "No usable gnomAD lookup for this allele");
        }
        PopulationEvidence::NoRecord => {
            return CriterionResult::not_met(code, "Absent from gnomAD");
        }
        PopulationEvidence::Record(g) => g,
    };
    if let Some(filter) = g.failed_filter() {
        return CriterionResult::not_evaluable(
            code,
            format!("gnomAD record failed the {filter} filter; frequency not usable"),
        );
    }
    let Some((pop, af)) = continental_max(g) else {
        return CriterionResult::not_evaluable(
            code,
            "gnomAD record has no frequency for any of the five continental populations",
        );
    };
    let details = json!({ "population": pop, "af": af, "threshold": threshold });
    if af > threshold {
        // Table 3's caveat is a cut-off question for a curator, not a switch
        // this tool can throw, so it travels with the result.
        CriterionResult::met(
            code,
            format!(
                "gnomAD {pop} AF {af:.4} exceeds {threshold}; if this gene causes hereditary \
                 cancer, its ClinGen expert panel cut-off should be consulted"
            ),
        )
        .with_details(details)
    } else {
        CriterionResult::not_met(
            code,
            format!("Highest continental gnomAD AF {af:.4} ({pop}) does not exceed {threshold}"),
        )
        .with_details(details)
    }
}

pub fn evaluate_sbvs1(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    evaluate_benign_frequency(CriterionCode::Sbvs1, config.sbvs1_af, input)
}

pub fn evaluate_sbs1(input: &OncogenicityInput, config: &OncogenicityConfig) -> CriterionResult {
    evaluate_benign_frequency(CriterionCode::Sbs1, config.sbs1_af, input)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Outcome;

    fn with_population(population: PopulationEvidence) -> OncogenicityInput {
        OncogenicityInput {
            population,
            ..Default::default()
        }
    }

    fn record(g: GnomadData) -> OncogenicityInput {
        with_population(PopulationEvidence::Record(g))
    }

    #[test]
    fn op4_is_not_evaluable_without_gnomad() {
        let r = evaluate_op4(&OncogenicityInput::default(), &Default::default());
        assert_eq!(r.outcome, Outcome::NotEvaluable);
    }

    #[test]
    fn op4_met_when_gnomad_loaded_and_no_record() {
        let r = evaluate_op4(
            &with_population(PopulationEvidence::NoRecord),
            &Default::default(),
        );
        assert_eq!(r.outcome, Outcome::Met);
        assert_eq!(r.points(), 1);
    }

    #[test]
    fn op4_met_on_zero_allele_record() {
        let r = evaluate_op4(
            &record(GnomadData {
                all_ac: Some(0),
                all_af: Some(0.0),
                ..Default::default()
            }),
            &Default::default(),
        );
        assert_eq!(r.outcome, Outcome::Met);
    }

    #[test]
    fn op4_default_requires_absence() {
        let rare = record(GnomadData {
            all_ac: Some(1),
            all_af: Some(4e-6),
            ..Default::default()
        });
        assert_eq!(
            evaluate_op4(&rare, &Default::default()).outcome,
            Outcome::NotMet
        );
        let cfg = OncogenicityConfig {
            op4_max_af: Some(1e-5),
            ..Default::default()
        };
        assert_eq!(evaluate_op4(&rare, &cfg).outcome, Outcome::Met);
    }

    #[test]
    fn filtered_record_is_not_evidence_either_way() {
        let g = GnomadData {
            all_af: Some(0.2),
            afr_af: Some(0.2),
            filter_vqsr: true,
            ..Default::default()
        };
        let input = record(g);
        assert_eq!(
            evaluate_op4(&input, &Default::default()).outcome,
            Outcome::NotEvaluable
        );
        assert_eq!(
            evaluate_sbvs1(&input, &Default::default()).outcome,
            Outcome::NotEvaluable
        );
    }

    #[test]
    fn sbvs1_and_sbs1_thresholds() {
        let cfg = OncogenicityConfig::default();
        let at = |af: f64| {
            record(GnomadData {
                nfe_af: Some(af),
                ..Default::default()
            })
        };
        // Table 3 says "is >5%" / "is >1%": the boundary itself does not count.
        assert_eq!(evaluate_sbvs1(&at(0.05), &cfg).outcome, Outcome::NotMet);
        assert_eq!(evaluate_sbvs1(&at(0.051), &cfg).outcome, Outcome::Met);
        assert_eq!(evaluate_sbvs1(&at(0.051), &cfg).points(), -8);
        assert_eq!(evaluate_sbs1(&at(0.01), &cfg).outcome, Outcome::NotMet);
        assert_eq!(evaluate_sbs1(&at(0.02), &cfg).outcome, Outcome::Met);
        assert_eq!(evaluate_sbs1(&at(0.02), &cfg).points(), -4);
    }

    #[test]
    fn benign_frequency_ignores_populations_the_sop_does_not_name() {
        // Common only in Finnish and Ashkenazi: founder populations Table 3
        // leaves out.
        let input = record(GnomadData {
            fin_af: Some(0.2),
            asj_af: Some(0.2),
            nfe_af: Some(0.001),
            ..Default::default()
        });
        let cfg = OncogenicityConfig::default();
        assert_eq!(evaluate_sbvs1(&input, &cfg).outcome, Outcome::NotMet);
        assert_eq!(evaluate_sbs1(&input, &cfg).outcome, Outcome::NotMet);
    }

    #[test]
    fn benign_frequency_without_continental_data_is_not_evaluable() {
        // An overall AF alone cannot say which population carries it.
        let input = record(GnomadData {
            all_af: Some(0.2),
            ..Default::default()
        });
        assert_eq!(
            evaluate_sbs1(&input, &Default::default()).outcome,
            Outcome::NotEvaluable
        );
    }

    #[test]
    fn absent_variant_is_not_benign_by_frequency() {
        let input = with_population(PopulationEvidence::NoRecord);
        assert_eq!(
            evaluate_sbs1(&input, &Default::default()).outcome,
            Outcome::NotMet
        );
    }
}
