use crate::structure::{OrdSignature, SetSignature};

/// A set of subsets of some set.
pub trait SubsetsSignature: SetSignature {}

/// A set of points in some space, defined in terms of a sequence of subsets whose intersection is the point being approximated.
pub trait ApproximatePointsSignature: SetSignature {
    type Precision: OrdSignature;
    type OpenSubsetsStructure: SubsetsSignature;

    /// An open subset containing the point
    fn open_neighbourhood(
        &self,
        approx_point: &Self::Set,
    ) -> <Self::OpenSubsetsStructure as SetSignature>::Set;
    /// A value which decreases with refinements of the subset containing the point.
    fn precision(&self, approx_point: &Self::Set) -> <Self::Precision as SetSignature>::Set;
    /// Refine approx_point until its precision is at most the given value.
    fn refine_to(
        &self,
        approx_point: Self::Set,
        precision: &<Self::Precision as SetSignature>::Set,
    ) -> Self::Set;
}
