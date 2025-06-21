use algebraeon_sets::structure::BorrowedStructure;

use crate::structure::{IdealsSignature, RingSignature};

pub trait LocalRingSignature: RingSignature {}

pub trait LocalRingIdealsSignature<Ring: LocalRingSignature, RingB: BorrowedStructure<Ring>>:
    IdealsSignature<Ring, RingB>
{
    /// A local ring has a unique maximal ideal so the corresponding `IdealsSignature`
    /// for such a ring has this to single out which of `Self::Set` representing all the ideals
    /// is that unique maximal ideal.
    /// `Self::Set` might encode ideals in a redundant manner so this just has to be one of them
    fn unique_maximal(&self) -> Self::Set;
}
