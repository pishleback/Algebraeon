use crate::structure::structure::*;
use algebraeon_sets::structure::*;

pub mod factorization;
pub mod functions;
pub mod primes;

use algebraeon_nzq::integer::*;
use algebraeon_nzq::natural::*;

impl SemiRingStructure for CannonicalStructure<Natural> {
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
