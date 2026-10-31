use fastvep_classification::sa_extract::{GnomadData, RevelData, SpliceAiData};
use fastvep_core::Consequence;

use crate::evidence::{
    EstablishedOncogenicVariants, FunctionalDomainEvidence, GeneRole, HotspotEvidence,
    OncogenicFunctionalEvidence,
};

/// Everything the oncogenicity criteria read, for one allele against one
/// transcript.
///
/// Two kinds of field live here. The first block is what the annotation
/// pipeline already computes and [`crate::extract_oncogenicity_input`] fills:
/// consequence, protein position, NMD, gnomAD, REVEL, SpliceAI, PhyloP. The
/// second block is somatic evidence no fastVEP source produces yet (see
/// [`crate::evidence`]); it stays `None` until a caller supplies it, and the
/// criteria that need it report `NotEvaluable`.
///
/// Tumor type is absent altogether. The SOP classifies "in the context of
/// relevant tumor type(s)", and OP2 cannot be asked without one; nothing in the
/// pipeline carries it.
/// //! TODO: Address Tumor Type Implementation
#[derive(Debug, Clone, Default)]
pub struct OncogenicityInput {
    // ── From the annotation pipeline ────────────────────────────────────
    pub consequences: Vec<Consequence>,
    pub gene_symbol: Option<String>,
    pub is_canonical: bool,
    /// Reference and alternate peptide as the consequence predictor writes
    /// them: one-letter codes, `*` for a stop, `-` for none.
    pub amino_acids: Option<(String, String)>,
    /// First affected residue.
    pub protein_position: Option<u64>,
    pub hgvs_c: Option<String>,
    /// Whether a premature stop here escapes nonsense-mediated decay by the
    /// 50-nt rule. `None` when the pipeline could not measure it.
    pub escapes_nmd: Option<bool>,
    pub population: PopulationEvidence,
    pub revel: Option<RevelData>,
    pub splice_ai: Option<SpliceAiData>,
    pub phylop: Option<f64>,

    // TODO: Missing evidence source integration
    // ── Somatic evidence with no pipeline source yet ────────────────────
    pub gene_role: Option<GeneRole>,
    pub hotspot: Option<HotspotEvidence>,
    pub established_oncogenic: Option<EstablishedOncogenicVariants>,
    pub functional_domain: Option<FunctionalDomainEvidence>,
    pub functional: Option<OncogenicFunctionalEvidence>,
}

/// What gnomAD says about this allele.
///
/// Three states because the absence of a gnomAD record means two different
/// things. Where the lookup is trustworthy, no record is the finding OP4 asks
/// for. Where it is not, a miss says nothing, and there are four ways for it
/// not to be:
///
/// - no gnomAD database is loaded, so every variant misses;
/// - no loaded gnomAD database holds this chromosome - per-chromosome
///   databases are the usual layout, and a database asked about a chromosome
///   it doesn't hold answers "no record" just as it does for an absent allele;
/// - the lookup failed;
/// - an indel was looked up without the reference normalisation gnomAD's keys
///   need, and can miss a record that exists.
///
/// Reading any of these as "absent from controls" would award OP4 on no
/// evidence.
///
/// One case is not caught: holding a chromosome is not covering every
/// position on it. A database built only from target regions, or a source
/// such as gnomAD exomes that was never sequenced outside its capture, can
/// miss a position it has no data for and still count as having answered.
/// Catching that would need a coverage track at annotation time. The germline PM2 has a config switch for the same ambiguity; here
/// the caller says which case it is instead.
#[derive(Debug, Clone, Default)]
// One per input and never collected, so the unused space costs nothing; boxing
// would cost an allocation per allele-transcript on the annotation hot path.
#[allow(clippy::large_enum_variant)]
pub enum PopulationEvidence {
    /// No record, and nothing to say the lookup could have found one.
    #[default]
    Unavailable,
    /// A gnomAD database holding this chromosome answered, in a spelling
    /// gnomAD stores, and has no record of this allele.
    NoRecord,
    Record(GnomadData),
}

/// The shape of change, which decides which criteria and which computational
/// lines apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantClass {
    /// Nonsense, frameshift, canonical ±1/2 splice site, or initiation codon:
    /// the null classes (OVS1)
    PredictedNull,
    /// In-frame insertion or deletion, stop loss, or another length-changing
    /// protein alteration (OM2).
    ProteinLengthChange,
    Missense,
    Synonymous,
    /// Everything that changes no amino acid and is not synonymous: intronic,
    /// splice region, UTR, non-coding, intergenic.
    NoProteinChange,
}

impl OncogenicityInput {
    pub fn variant_class(&self) -> VariantClass {
        if self.consequences.iter().any(is_null_consequence) {
            VariantClass::PredictedNull
        } else if self.consequences.iter().any(is_protein_length_change) {
            VariantClass::ProteinLengthChange
        } else if self.consequences.contains(&Consequence::MissenseVariant) {
            VariantClass::Missense
        } else if self.consequences.contains(&Consequence::SynonymousVariant) {
            VariantClass::Synonymous
        } else {
            VariantClass::NoProteinChange
        }
    }

    /// The variant's change as a single-residue substitution, `(ref, alt)`,
    /// when it is one. Multi-residue and indel peptides return `None`.
    pub fn single_residue_change(&self) -> Option<(&str, &str)> {
        let (r, a) = self.amino_acids.as_ref()?;
        (r.len() == 1 && a.len() == 1 && r != "-" && a != "-" && r != a)
            .then_some((r.as_str(), a.as_str()))
    }
}

/// Consequence terms for the null variant classes OVS1 lists.
/// `splice_donor_variant` / `splice_acceptor_variant` are assigned on the two
/// intronic bases at each end of an intron (fastvep-consequence `splice.rs`),
/// which is OVS1's "canonical ±1 or 2 splice sites"; the wider splice-region
/// terms are separate and not null.
///
/// OVS1 also names single- and multi-exon deletions. Those arrive through the
/// structural-variant path, not as these terms, and are not assessed here.
pub(crate) fn is_null_consequence(c: &Consequence) -> bool {
    matches!(
        c,
        Consequence::StopGained
            | Consequence::FrameshiftVariant
            | Consequence::SpliceDonorVariant
            | Consequence::SpliceAcceptorVariant
            | Consequence::StartLost
    )
}

pub(crate) fn is_protein_length_change(c: &Consequence) -> bool {
    matches!(
        c,
        Consequence::InframeInsertion
            | Consequence::InframeDeletion
            | Consequence::StopLost
            | Consequence::ProteinAlteringVariant
    )
}
