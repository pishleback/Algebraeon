use algebraeon_macros::CanonicalStructure;
use algebraeon_structures::*;
use cantor::Finite;

/// Given an ordered syntheme `ab cd ef` specify either the left or right element of a pair i.e. `{a, c, e}` vs `{b, d, f}`
#[derive(CanonicalStructure, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Finite)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub enum OrderedSynthemeSide {
    Left,
    Right,
}

impl ConstSizeFiniteSetSignature<2> for OrderedSynthemeSideCanonicalStructure {}

/// Given an ordered syntheme `ab cd ef` specify a pair i.e. `ab`, `cd`, or `ef`.
#[derive(CanonicalStructure, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Finite)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub enum OrderedSynthemePair {
    Left,
    Middle,
    Right,
}

impl ConstSizeFiniteSetSignature<3> for OrderedSynthemePairCanonicalStructure {}

/// Given an ordered syntheme `ab cd ef` specify one of its points i.e. `a`, `b`, `c`, `d`, `e`, or `f`.
#[derive(CanonicalStructure, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Finite)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub struct OrderedSynthemePoint {
    pub pair: OrderedSynthemePair,
    pub side: OrderedSynthemeSide,
}

impl ConstSizeFiniteSetSignature<6> for OrderedSynthemePointCanonicalStructure {}
