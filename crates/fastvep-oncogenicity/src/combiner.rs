//! Points to classification.
//!
//! Horak 2022 combines evidence by summing points (Supporting 1, Moderate 2,
//! Strong 4, Very Strong 8; negative for benign). The bands below are from
//! the paper's Figure 1.

use crate::types::OncogenicityClassification;

pub fn classify_points(points: i32) -> OncogenicityClassification {
    match points {
        10.. => OncogenicityClassification::Oncogenic,
        6..=9 => OncogenicityClassification::LikelyOncogenic,
        0..=5 => OncogenicityClassification::UncertainSignificance,
        -6..=-1 => OncogenicityClassification::LikelyBenign,
        _ => OncogenicityClassification::Benign,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use OncogenicityClassification::*;

    #[test]
    fn worked_examples_from_the_paper() {
        // Example 2, BRAF p.Gly469Glu: OS2 + OP1 + OP3 + OP4 = 7.
        assert_eq!(classify_points(7), LikelyOncogenic);
        // Example 1, KRAS p.Arg164Gln: SBS2 (downgraded to -1) + OP1 = 0.
        assert_eq!(classify_points(0), UncertainSignificance);
    }

    #[test]
    fn bands() {
        assert_eq!(classify_points(10), Oncogenic);
        assert_eq!(classify_points(9), LikelyOncogenic);
        assert_eq!(classify_points(6), LikelyOncogenic);
        assert_eq!(classify_points(5), UncertainSignificance);
        assert_eq!(classify_points(-1), LikelyBenign);
        assert_eq!(classify_points(-6), LikelyBenign);
        assert_eq!(classify_points(-7), Benign);
    }
}
