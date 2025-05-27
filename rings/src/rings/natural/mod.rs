use crate::structure::*;
use algebraeon_nzq::{traits::DivMod, *};
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

impl CharacteristicSignature for NaturalCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl SemiRingUnitsSignature for NaturalCanonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match *a {
            Natural::ZERO => Err(RingDivisionError::DivideByZero),
            Natural::ONE => Ok(Natural::ONE),
            _ => Err(RingDivisionError::NotDivisible),
        }
    }
}

impl EuclideanDivisionSignature  for NaturalCanonicalStructure {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if elem == &Natural::ZERO {
            None
        } else {
            Some(elem.clone())
        }
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        if b == &Natural::ZERO {
            None
        } else {
            Some(a.div_mod(b))
        }
    }
}
