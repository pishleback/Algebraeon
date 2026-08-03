use std::sync::Arc;

use crate::structure::*;
use algebraeon_structures::*;

impl RinglikeSpecializationSignature for NaturalCanonicalStructure {}

impl ZeroSignature for NaturalCanonicalStructure {
    fn zero(self: &Arc<Self>) -> Self::Elem {
        Natural::ZERO
    }
}

impl AdditionSignature for NaturalCanonicalStructure {
    fn add(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a + b
    }
}

impl CancellativeAdditionSignature for NaturalCanonicalStructure {
    fn try_sub(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        a.try_sub(b)
    }
}

impl TryNegateSignature for NaturalCanonicalStructure {
    fn try_neg(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        let z = self.zero();
        if a == &z { Some(self.zero()) } else { None }
    }
}

impl AdditiveMonoidSignature for NaturalCanonicalStructure {}

impl OneSignature for NaturalCanonicalStructure {
    fn one(self: &Arc<Self>) -> Self::Elem {
        Natural::ONE
    }
}

impl MultiplicationSignature for NaturalCanonicalStructure {
    fn mul(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a * b
    }
}

impl CommutativeMultiplicationSignature for NaturalCanonicalStructure {}

impl MultiplicativeMonoidSignature for NaturalCanonicalStructure {}

impl FavoriteAssociateSignature for NaturalCanonicalStructure {
    fn factor_fav_assoc(self: &Arc<Self>, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        (Natural::ONE, a.clone())
    }
}

impl MultiplicativeAbsorptionMonoidSignature for NaturalCanonicalStructure {}

impl LeftDistributiveMultiplicationOverAddition for NaturalCanonicalStructure {}

impl RightDistributiveMultiplicationOverAddition for NaturalCanonicalStructure {}

impl SemiRingSignature for NaturalCanonicalStructure {}

impl CharacteristicSignature for NaturalCanonicalStructure {
    fn characteristic(self: &Arc<Self>) -> Natural {
        Natural::ZERO
    }
}

impl TryReciprocalSignature for NaturalCanonicalStructure {
    fn try_reciprocal(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        match *a {
            Natural::ZERO => None,
            Natural::ONE => Some(Natural::ONE),
            _ => None,
        }
    }
}

impl EuclideanDivisionSignature for NaturalCanonicalStructure {
    fn norm(self: &Arc<Self>, elem: &Self::Elem) -> Option<Natural> {
        if elem == &Natural::ZERO {
            None
        } else {
            Some(elem.clone())
        }
    }

    fn quorem(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<(Self::Elem, Self::Elem)> {
        if b == &Natural::ZERO {
            None
        } else {
            Some(a.div_mod(b))
        }
    }
}
