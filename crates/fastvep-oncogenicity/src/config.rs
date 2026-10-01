//! Thresholds for the criteria whose SOP text states no number.
//!
//! OP1, SBP1 and SBP2 ask whether computational evidence "supports an oncogenic
//! effect" or "suggests no impact", and OP4 asks whether a variant is absent
//! or "at an extremely low frequency". Horak 2022 gives none of these a
//! numeric threshold. The defaults below are therefore choices, not SOP text:
//! the predictor cut-offs are the ClinGen SVI calibrations the germline
//! classifier already applies (`AcmgConfig::default()` in
//! fastvep-classification), and OP4 defaults to absence only. They sit in one
//! place so they can be changed without touching criterion logic.
//!
//! SBVS1's 5% and SBS1's 1% are SOP text (Table 3).

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct OncogenicityConfig {
    /// REVEL at or above which the missense line supports an effect. The
    /// germline PP3 Supporting default (ClinGen SVI calibration).
    pub revel_impact: f64,
    /// REVEL at or below which the missense line suggests no impact. The
    /// germline BP4 Supporting default (ClinGen SVI calibration).
    pub revel_no_impact: f64,
    /// Maximum SpliceAI delta at or above which the splicing line supports an
    /// effect. Walker 2023, as in the germline default.
    pub spliceai_impact: f64,
    /// Maximum SpliceAI delta at or below which the splicing line suggests no
    /// impact. Walker 2023, as in the germline default.
    pub spliceai_no_impact: f64,
    /// PhyloP at or above which a nucleotide counts as highly conserved. The
    /// germline `phylop_conserved` default, which carries no citation there
    /// either.
    pub phylop_conserved: f64,
    /// Highest gnomAD allele frequency OP4 accepts as "extremely low". `None`
    /// accepts only absence, which is the only reading the SOP text settles.
    pub op4_max_af: Option<f64>,
    /// SBVS1: allele frequency above which a variant is benign, Very Strong.
    pub sbvs1_af: f64,
    /// SBS1: allele frequency above which a variant is benign, Strong.
    pub sbs1_af: f64,
}

impl Default for OncogenicityConfig {
    fn default() -> Self {
        Self {
            revel_impact: 0.644,
            revel_no_impact: 0.290,
            spliceai_impact: 0.2,
            spliceai_no_impact: 0.1,
            phylop_conserved: 2.0,
            op4_max_af: None,
            sbvs1_af: 0.05,
            sbs1_af: 0.01,
        }
    }
}
