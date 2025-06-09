use algebraeon_sets::structure::BorrowedStructure;

use crate::structure::{IdealsSignature, RingSignature};

pub trait LocalRingSignature: RingSignature {}

pub trait LocalRingIdealsSignature<Ring: LocalRingSignature, RingB: BorrowedStructure<Ring>>:
    IdealsSignature<Ring, RingB>
{
    fn unique_maximal(&self) -> Self::Set;
}
