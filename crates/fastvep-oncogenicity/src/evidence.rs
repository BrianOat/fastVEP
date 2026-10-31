//! Somatic evidence the annotation pipeline does not produce.
//!
//! Every type here answers a question no current fastVEP source can: whether a
//! gene is a tumor suppressor, whether a residue is a cancerhotspots.org
//! hotspot, which variants have already been classified oncogenic under the
//! SOP, where a protein's critical functional regions are, and what curated
//! functional studies found. gnomAD constraint says a gene is intolerant of loss of
//! function, not that it suppresses tumors. A COSMIC allele count is not a hotspot call.
//!
//! The criteria that read these fields are implemented and tested against them.
//! [`crate::extract_oncogenicity_input`] leaves them `None`, and each such
//! criterion reports `NotEvaluable`, until a caller has a real source to fill
//! them from.
//!
//! Where a lookup can be run and find nothing, the type has an explicit
//! negative variant, so that "consulted, found nothing" stays distinct from
//! "not consulted" (`None`).

use serde::{Deserialize, Serialize};

use crate::types::Strength;

/// The gene's role in cancer, as curated for the tumor type under
/// classification.
///
/// Both flags can be set (as both can be true depending on tumor type).
///
/// The SOP's EZH2 example is a gene whose role depends
/// on tumor type, which is why this belongs to the classification context
/// rather than to the gene.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct GeneRole {
    pub tumor_suppressor: bool,
    pub oncogene: bool,
}

/// What cancerhotspots.org says about the residue this variant changes.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum HotspotEvidence {
    /// The hotspot index was consulted and the residue is not a hotspot.
    NotHotspot,
    Hotspot {
        /// Samples with any somatic variant at this amino-acid position.
        residue_samples: u32,
        /// Samples with this exact amino-acid change.
        aa_change_samples: u32,
        /// Whether the hotspot evidence is resolved to this variant's exact nucleotide
        /// change. OS3 may be counted alongside OS1 only when it is (Table 2).
        /// Therefore, this should be true only when the hotspot evidence specifically
        /// supports the nucleotide change being classified.
        nucleotide_resolved: bool,
    },
}

/// Variants at this residue already classified Oncogenic or Likely Oncogenic
/// under this SOP, from whatever store holds such classifications.
///
/// OS1: another variant at this residue produces the SAME amino-acid change,
/// regardless of whether the nucleotide change is the same.
///
/// OM4: another oncogenic variant exists at this residue but produces a
/// DIFFERENT amino-acid change; the physicochemical difference must then be
/// evaluated according to the SOP.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct EstablishedOncogenicVariants {
    /// Entries at the residue this variant changes. Empty means the store was
    /// consulted and holds none.
    pub at_residue: Vec<EstablishedVariant>,
}

/// One previously classified oncogenic variant.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EstablishedVariant {
    /// Reference amino acid, one-letter code.
    pub ref_aa: String,
    /// Alternate amino acid, one-letter code (`*` for a stop).
    pub alt_aa: String,
    /// True when the entry is the variant being classified: the same nucleotide
    /// change. It is not evidence for itself; OS1 asks about another variant
    /// producing the same amino-acid change.
    pub is_self: bool,
}

/// Whether the variant lies in a critical, well-established part of a
/// functional domain (OM1), such as an enzyme's active site.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum FunctionalDomainEvidence {
    /// The curated domain map was consulted and the residue is outside every
    /// critical region in it.
    NotInCriticalRegion,
    InCriticalRegion {
        /// The region, for the summary: e.g. "kinase activation loop".
        region: String,
    },
}

/// A curated result from well-established in vitro or in vivo functional
/// studies (OS2 / SBS2).
///
/// Separate from the germline `FunctionalEvidence` in fastvep-classification on
/// purpose. OS2/SBS2 ask whether the variant has an oncogenic effect.
/// For an activating oncogene variant the answers point in opposite directions,
/// so one file cannot serve both.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OncogenicFunctionalEvidence {
    pub effect: FunctionalEffect,
    /// Strength the curator assigned. The SOP lets curators downgrade OS2 and
    /// SBS2 when evidence is limited or partly conflicting, never upgrade.
    pub strength: Strength,
    /// Whether the studies used this variant's exact nucleotide change. When
    /// OS1 also applies, OS2 may be counted only if they did (Table 2).
    pub nucleotide_specific: bool,
    /// Citation or note surfaced in the summary.
    pub source: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FunctionalEffect {
    /// Studies support an oncogenic effect (OS2).
    Oncogenic,
    /// Studies show no oncogenic effect (SBS2).
    NoOncogenicEffect,
}
