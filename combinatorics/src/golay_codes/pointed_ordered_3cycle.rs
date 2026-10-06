use algebraeon_macros::CanonicalStructure;
use algebraeon_structures::*;
use cantor::Finite;

/// A 4-element set consisting of an ordered 3-cycle and a fourth distinguished point
#[derive(CanonicalStructure, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Finite)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub enum PointedOrdered3Cycle {
    P,
    C1,
    C2,
    C3,
}

impl ConstSizeFiniteSetSignature<4> for PointedOrdered3CycleCanonicalStructure {}
