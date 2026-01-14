use crate::structure::*;
use algebraeon_nzq::{traits::DivMod, *};
use algebraeon_sets::structure::*;

pub mod factorization;
pub mod functions;

impl RinglikeSpecializationSignature for NaturalCanonicalStructure {}

impl ZeroSignature for NaturalCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Natural::ZERO
    }
}

impl AdditionSignature for NaturalCanonicalStructure {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
}

impl CancellativeAdditionSignature for NaturalCanonicalStructure {
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        a.try_sub(b)
    }
}

impl TryNegateSignature for NaturalCanonicalStructure {
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        let z = self.zero();
        if a == &z { Some(self.zero()) } else { None }
    }
}

impl AdditiveMonoidSignature for NaturalCanonicalStructure {}

impl OneSignature for NaturalCanonicalStructure {
    fn one(&self) -> Self::Set {
        Natural::ONE
    }
}

impl MultiplicationSignature for NaturalCanonicalStructure {
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl MultiplicativeMonoidSignature for NaturalCanonicalStructure {}

impl FavoriteAssociateSignature for NaturalCanonicalStructure {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        (Natural::ONE, a.clone())
    }
}

impl SemiRingSignature for NaturalCanonicalStructure {}

impl CharacteristicSignature for NaturalCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl MultiplicativeMonoidUnitsSignature for NaturalCanonicalStructure {
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set> {
        match *a {
            Natural::ZERO => None,
            Natural::ONE => Some(Natural::ONE),
            _ => None,
        }
    }
}

impl EuclideanDivisionSignature for NaturalCanonicalStructure {
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

// pub trait NaturalFns {
//     fn is_prime(&self) -> bool;
//     fn factor(self) -> Option<Vec<(Natural, Natural)>>;
// }

// impl NaturalFns for Natural {
//     fn is_prime(&self) -> bool {
//         factorization::primes::is_prime(self)
//     }
//     fn factor(self) -> Option<Vec<(Natural, Natural)>> {
//         factorization::factor(self)
//     }
// }
