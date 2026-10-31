use serde::{Deserialize, Serialize};

/// Evidence strength of an oncogenicity criterion, with the point weight
/// Horak 2022 Tables 2 and 3 assign to it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Strength {
    Supporting,
    Moderate,
    Strong,
    VeryStrong,
}

impl Strength {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Supporting => "Supporting",
            Self::Moderate => "Moderate",
            Self::Strong => "Strong",
            Self::VeryStrong => "Very_Strong",
        }
    }

    /// Unsigned weight: Supporting 1, Moderate 2, Strong 4, Very Strong 8.
    pub fn points(&self) -> i32 {
        match self {
            Self::Supporting => 1,
            Self::Moderate => 2,
            Self::Strong => 4,
            Self::VeryStrong => 8,
        }
    }
}

/// Whether a criterion argues for or against oncogenicity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Direction {
    Oncogenic,
    Benign,
}

/// The seventeen criteria of Horak 2022 (Table 2: oncogenic, Table 3: somatic
/// benign).
///
/// An enum rather than a string because the SOP's exclusions are relations
/// between criteria, and a misspelt code in an exclusion table would silently
/// never fire.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CriterionCode {
    #[serde(rename = "OVS1")]
    Ovs1,
    #[serde(rename = "OS1")]
    Os1,
    #[serde(rename = "OS2")]
    Os2,
    #[serde(rename = "OS3")]
    Os3,
    #[serde(rename = "OM1")]
    Om1,
    #[serde(rename = "OM2")]
    Om2,
    #[serde(rename = "OM3")]
    Om3,
    #[serde(rename = "OM4")]
    Om4,
    #[serde(rename = "OP1")]
    Op1,
    #[serde(rename = "OP2")]
    Op2,
    #[serde(rename = "OP3")]
    Op3,
    #[serde(rename = "OP4")]
    Op4,
    #[serde(rename = "SBVS1")]
    Sbvs1,
    #[serde(rename = "SBS1")]
    Sbs1,
    #[serde(rename = "SBS2")]
    Sbs2,
    #[serde(rename = "SBP1")]
    Sbp1,
    #[serde(rename = "SBP2")]
    Sbp2,
}

impl CriterionCode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Ovs1 => "OVS1",
            Self::Os1 => "OS1",
            Self::Os2 => "OS2",
            Self::Os3 => "OS3",
            Self::Om1 => "OM1",
            Self::Om2 => "OM2",
            Self::Om3 => "OM3",
            Self::Om4 => "OM4",
            Self::Op1 => "OP1",
            Self::Op2 => "OP2",
            Self::Op3 => "OP3",
            Self::Op4 => "OP4",
            Self::Sbvs1 => "SBVS1",
            Self::Sbs1 => "SBS1",
            Self::Sbs2 => "SBS2",
            Self::Sbp1 => "SBP1",
            Self::Sbp2 => "SBP2",
        }
    }

    pub fn direction(&self) -> Direction {
        match self {
            Self::Sbvs1 | Self::Sbs1 | Self::Sbs2 | Self::Sbp1 | Self::Sbp2 => Direction::Benign,
            _ => Direction::Oncogenic,
        }
    }

    /// Strength the SOP assigns the criterion before any downgrade.
    pub fn default_strength(&self) -> Strength {
        match self {
            Self::Ovs1 | Self::Sbvs1 => Strength::VeryStrong,
            Self::Os1 | Self::Os2 | Self::Os3 | Self::Sbs1 | Self::Sbs2 => Strength::Strong,
            Self::Om1 | Self::Om2 | Self::Om3 | Self::Om4 => Strength::Moderate,
            Self::Op1 | Self::Op2 | Self::Op3 | Self::Op4 | Self::Sbp1 | Self::Sbp2 => {
                Strength::Supporting
            }
        }
    }
}

/// What evaluating a criterion established.
///
/// Four states because three different facts hide behind "did not fire".
/// `NotMet` is an answer: the evidence was consulted and does not satisfy the
/// criterion. `NotEvaluable` is the absence of one: the evidence the criterion
/// needs is not available to this run, and nothing may be inferred from that
/// about the variant. `Excluded` means the criterion's own conditions hold but
/// the SOP forbids counting it alongside another criterion that also holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    Met,
    NotMet,
    NotEvaluable,
    Excluded,
}

/// The evaluation of one criterion for one allele-transcript pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriterionResult {
    pub code: CriterionCode,
    pub direction: Direction,
    /// Strength the SOP assigns the criterion.
    pub default_strength: Strength,
    /// Strength applied. Differs from `default_strength` only by a downgrade,
    /// which the SOP permits for OVS1, OS2, OM1, OM2, OM4 and SBS2.
    pub strength: Strength,
    pub outcome: Outcome,
    /// For `Outcome::Excluded`, the met criteria that excluded this one.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub excluded_by: Vec<CriterionCode>,
    /// Why the outcome is what it is, in terms a reader can check.
    pub summary: String,
    /// The data points the outcome rests on.
    #[serde(default, skip_serializing_if = "serde_json::Value::is_null")]
    pub details: serde_json::Value,
}

impl CriterionResult {
    pub(crate) fn new(code: CriterionCode, outcome: Outcome, summary: impl Into<String>) -> Self {
        Self {
            code,
            direction: code.direction(),
            default_strength: code.default_strength(),
            strength: code.default_strength(),
            outcome,
            excluded_by: Vec::new(),
            summary: summary.into(),
            details: serde_json::Value::Null,
        }
    }

    pub(crate) fn met(code: CriterionCode, summary: impl Into<String>) -> Self {
        Self::new(code, Outcome::Met, summary)
    }

    pub(crate) fn not_met(code: CriterionCode, summary: impl Into<String>) -> Self {
        Self::new(code, Outcome::NotMet, summary)
    }

    pub(crate) fn not_evaluable(code: CriterionCode, summary: impl Into<String>) -> Self {
        Self::new(code, Outcome::NotEvaluable, summary)
    }

    pub(crate) fn with_strength(mut self, strength: Strength) -> Self {
        self.strength = strength;
        self
    }

    pub(crate) fn with_details(mut self, details: serde_json::Value) -> Self {
        self.details = details;
        self
    }

    /// Signed points this criterion contributes. Only a `Met` criterion
    /// contributes; an excluded one holds but is not counted.
    pub fn points(&self) -> i32 {
        if self.outcome != Outcome::Met {
            return 0;
        }
        match self.direction {
            Direction::Oncogenic => self.strength.points(),
            Direction::Benign => -self.strength.points(),
        }
    }
}

/// Five-tier oncogenicity classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum OncogenicityClassification {
    Oncogenic,
    LikelyOncogenic,
    UncertainSignificance,
    LikelyBenign,
    Benign,
}

impl OncogenicityClassification {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Oncogenic => "Oncogenic",
            Self::LikelyOncogenic => "Likely_oncogenic",
            Self::UncertainSignificance => "Uncertain_significance",
            Self::LikelyBenign => "Likely_benign",
            Self::Benign => "Benign",
        }
    }

    pub fn shorthand(&self) -> &'static str {
        match self {
            Self::Oncogenic => "O",
            Self::LikelyOncogenic => "LO",
            Self::UncertainSignificance => "VUS",
            Self::LikelyBenign => "LB",
            Self::Benign => "B",
        }
    }
}

/// Classification of one allele-transcript pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OncogenicityResult {
    pub classification: OncogenicityClassification,
    pub shorthand: String,
    /// Sum of the points of every `Met` criterion.
    pub points: i32,
    /// Every criterion, in SOP table order, whatever its outcome.
    pub criteria: Vec<CriterionResult>,
}
