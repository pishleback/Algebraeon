use crate::structure::*;
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegerIdealsStructure {}

impl IntegerIdealsStructure {
    fn new() -> Arc<Self> {
        Self {}.into()
    }
}

impl RingToIdealsSignature for IntegerCanonicalStructure {
    type Ideals = IntegerIdealsStructure;

    fn ideals(self: &Arc<Self>) -> Arc<Self::Ideals> {
        IntegerIdealsStructure::new()
    }
}

impl Signature for IntegerIdealsStructure {}

impl SetSignature for IntegerIdealsStructure {
    type Elem = Natural;
    fn validate_element(self: &Arc<Self>, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl IdealsSignature<IntegerCanonicalStructure> for IntegerIdealsStructure {
    fn ring(self: &Arc<Self>) -> Arc<IntegerCanonicalStructure> {
        Integer::structure()
    }
}

impl EqSignature for IntegerIdealsStructure {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl RinglikeSpecializationSignature for IntegerIdealsStructure {}

impl ZeroSignature for IntegerIdealsStructure {
    fn zero(self: &Arc<Self>) -> Self::Elem {
        Natural::ZERO
    }
}

impl AdditionSignature for IntegerIdealsStructure {
    fn add(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        gcd(a.clone(), b.clone())
    }
}

impl AdditiveMonoidSignature for IntegerIdealsStructure {}

impl TryNegateSignature for IntegerIdealsStructure {
    fn try_neg(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        if self.is_zero(a) {
            Some(self.zero())
        } else {
            None
        }
    }
}

impl OneSignature for IntegerIdealsStructure {
    fn one(self: &Arc<Self>) -> Self::Elem {
        Natural::ONE
    }
}

impl MultiplicationSignature for IntegerIdealsStructure {
    fn mul(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a * b
    }
}

impl CommutativeMultiplicationSignature for IntegerIdealsStructure {}

impl MultiplicativeMonoidSignature for IntegerIdealsStructure {}

impl MultiplicativeAbsorptionMonoidSignature for IntegerIdealsStructure {}

impl LeftDistributiveMultiplicationOverAddition for IntegerIdealsStructure {}

impl RightDistributiveMultiplicationOverAddition for IntegerIdealsStructure {}

impl SemiRingSignature for IntegerIdealsStructure {}

impl IdealsArithmeticSignature<IntegerCanonicalStructure> for IntegerIdealsStructure {
    fn principal_ideal(self: &Arc<Self>, a: &Integer) -> Self::Elem {
        Abs::abs(a)
    }

    fn contains_ideal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        b % a == Natural::ZERO
    }

    fn intersect(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        lcm(a.clone(), b.clone())
    }

    fn quotient(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        if b == &Natural::ZERO {
            Natural::ONE
        } else {
            a / gcd(a.clone(), b.clone())
        }
    }
}

impl PrincipalIdealsSignature<IntegerCanonicalStructure> for IntegerIdealsStructure {
    fn ideal_generator(self: &Arc<Self>, ideal: &Natural) -> Integer {
        Integer::from(ideal)
    }
}

impl DedekindDomainIdealsSignature<IntegerCanonicalStructure> for IntegerIdealsStructure {}

impl TryReciprocalSignature for IntegerIdealsStructure {
    fn try_reciprocal(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.factorization_exponents().try_reciprocal(a)
    }
}

impl FavoriteAssociateSignature for IntegerIdealsStructure {
    fn factor_fav_assoc(self: &Arc<Self>, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        self.factorization_exponents().factor_fav_assoc(a)
    }
}

impl UniqueFactorizationMonoidSignature for IntegerIdealsStructure {
    type FactoredExponent = NaturalCanonicalStructure;

    fn factorization_exponents(self: &Arc<Self>) -> Arc<Self::FactoredExponent> {
        Natural::structure()
    }

    fn factorization_pow(self: &Arc<Self>, a: &Self::Elem, k: &Natural) -> Self::Elem {
        self.factorization_exponents().nat_pow(a, k)
    }

    fn try_is_irreducible(self: &Arc<Self>, a: &Self::Elem) -> Option<bool> {
        self.factorization_exponents().try_is_irreducible(a)
    }
}

impl FactoringMonoidSignature for IntegerIdealsStructure {
    fn factor_unchecked(self: &Arc<Self>, ideal: &Natural) -> Factored<Natural, Natural> {
        Natural::structure().factor_unchecked(ideal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_ideals() {
        let ideals = Integer::structure().ideals();

        assert!(ideals.equal(&Natural::from(3u32), &Natural::from(3u32)));

        assert!(!ideals.equal(&Natural::from(2u32), &Natural::from(3u32)));

        assert!(ideals.contains_ideal(&Natural::from(2u32), &Natural::from(6u32)));

        assert!(!ideals.contains_ideal(&Natural::from(6u32), &Natural::from(2u32)));

        assert!(!ideals.contains_ideal(&Natural::from(5u32), &Natural::from(7u32)));

        assert!(ideals.equal(
            &ideals.add(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(3u32)
        ));

        assert!(
            ideals.equal(
                &Integer::structure()
                    .ideals()
                    .intersect(&Natural::from(6u32), &Natural::from(15u32)),
                &Natural::from(30u32)
            )
        );

        assert!(
            ideals.equal(
                &Integer::structure()
                    .ideals()
                    .mul(&Natural::from(6u32), &Natural::from(15u32)),
                &Natural::from(90u32)
            )
        );

        assert_eq!(ideals.generated_ideal(vec![-15, 6]), 3u32.into());
    }

    #[test]
    fn factor_integer_ideal() {
        let ideals = Integer::structure().ideals();

        let f = ideals.factor(
            &Integer::structure()
                .ideals()
                .principal_ideal(&Integer::from(0)),
        );
        println!("{:?}", f);
        assert!(f.is_zero());

        let f = ideals.factor(
            &Integer::structure()
                .ideals()
                .principal_ideal(&Integer::from(18)),
        );
        println!("{:?}", f);
        assert!(!f.is_zero());
    }
}
