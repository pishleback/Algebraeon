use algebraeon_sets::structure::{OrdSignature, SetSignature, Signature};

pub mod complex_plane;
pub mod open_box;
pub mod open_interval;
pub mod real_line;

pub trait TopologicalSpaceOpenSubsetsSignature: SetSignature {}

pub trait TopologicalSpaceApproximatePointsSignature: SetSignature {
    type Precision: OrdSignature;
    type OpenSubsetsStructure: TopologicalSpaceOpenSubsetsSignature;

    /// An open subset containing the point
    fn open_neighbourhood(
        &self,
        approx_point: &Self::Set,
    ) -> <Self::OpenSubsetsStructure as SetSignature>::Set;
    /// A value which strictly decreases with refinements.
    fn precision(&self, approx_point: &Self::Set) -> <Self::Precision as SetSignature>::Set;
    /// Refine approx_point until its precision is at most the given value.
    fn refine_to(
        &self,
        approx_point: Self::Set,
        precision: &<Self::Precision as SetSignature>::Set,
    ) -> Self::Set;
}
