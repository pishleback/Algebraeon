use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

pub mod factorization;
pub mod functions;

impl SemiRingSignature for NaturalCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Natural::ZERO
    }
    fn one(&self) -> Self::Set {
        Natural::ONE
    }
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}
