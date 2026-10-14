use crate::algebraic_number_field::AlgebraicIntegerRingSignature;
use crate::structure::*;
use algebraeon_structures::*;
use std::collections::HashSet;
use std::sync::Arc;

impl RinglikeSpecializationSignature for IntegerCanonicalStructure {
    fn try_ring_restructure(
        self: Arc<Self>,
    ) -> Option<Arc<impl EqSignature<Elem = Self::Elem> + RingSignature>> {
        Some(self)
    }

    fn try_char_zero_ring_restructure(
        self: Arc<Self>,
    ) -> Option<Arc<impl EqSignature<Elem = Self::Elem> + CharZeroRingSignature>> {
        Some(self)
    }
}

impl ZeroSignature for IntegerCanonicalStructure {
    fn zero(self: &Arc<Self>) -> Self::Elem {
        Integer::ZERO
    }
}

impl AdditionSignature for IntegerCanonicalStructure {
    fn add<'a>(
        self: &Arc<Self>,
        a: impl Into<Arg<'a, Self::Elem>>,
        b: impl Into<Arg<'a, Self::Elem>>,
    ) -> Self::Elem
    where
        Self: 'a,
        Self::Elem: 'a,
    {
        match (a.into(), b.into()) {
            (Arg::Borrowed(a), Arg::Borrowed(b)) => a + b,
            (Arg::Borrowed(a), Arg::Owned(b)) => a + b,
            (Arg::Owned(a), Arg::Borrowed(b)) => a + b,
            (Arg::Owned(a), Arg::Owned(b)) => a + b,
        }
    }
}

impl CancellativeAdditionSignature for IntegerCanonicalStructure {
    fn try_sub(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl TryNegateSignature for IntegerCanonicalStructure {
    fn try_neg(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl AdditiveMonoidSignature for IntegerCanonicalStructure {}

impl AdditiveGroupSignature for IntegerCanonicalStructure {
    fn neg(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        -a
    }

    fn sub(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a - b
    }
}

impl OneSignature for IntegerCanonicalStructure {
    fn one(self: &Arc<Self>) -> Self::Elem {
        Integer::ONE
    }
}

impl MultiplicationSignature for IntegerCanonicalStructure {
    fn mul(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a * b
    }
}

impl CommutativeMultiplicationSignature for IntegerCanonicalStructure {}

impl MultiplicativeMonoidSignature for IntegerCanonicalStructure {}

impl MultiplicativeAbsorptionMonoidSignature for IntegerCanonicalStructure {}

impl LeftDistributiveMultiplicationOverAddition for IntegerCanonicalStructure {}

impl RightDistributiveMultiplicationOverAddition for IntegerCanonicalStructure {}

impl SemiRingSignature for IntegerCanonicalStructure {}

impl RingSignature for IntegerCanonicalStructure {
    fn is_reduced(self: &Arc<Self>) -> Result<bool, String> {
        Ok(true)
    }
}

impl CharacteristicSignature for IntegerCanonicalStructure {
    fn characteristic(self: &Arc<Self>) -> Natural {
        Natural::ZERO
    }
}

impl TryReciprocalSignature for IntegerCanonicalStructure {
    fn try_reciprocal(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.try_divide(&self.one(), a)
    }
}

impl CancellativeMultiplicationSignature for IntegerCanonicalStructure {
    fn try_divide(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        match self.quorem(a, b) {
            Some((q, r)) => {
                if r == self.zero() {
                    Some(q)
                } else {
                    None
                }
            }
            None => None,
        }
    }
}

impl MultiplicativeIntegralMonoidSignature for IntegerCanonicalStructure {}

impl IntegralDomainSignature for IntegerCanonicalStructure {}

impl OrderedRingSignature for IntegerCanonicalStructure {}

impl CountableSetSignature for MultiplicativeMonoidUnitsStructure<IntegerCanonicalStructure> {
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl FiniteSetSignature for MultiplicativeMonoidUnitsStructure<IntegerCanonicalStructure> {
    fn list_all_elements(self: &Arc<Self>) -> Vec<Self::Elem> {
        vec![Integer::ONE, -Integer::ONE]
    }
}

impl FavoriteAssociateSignature for IntegerCanonicalStructure {
    fn factor_fav_assoc(self: &Arc<Self>, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        #[allow(clippy::comparison_chain)]
        if a == &Integer::ZERO {
            (Integer::ONE, Integer::ZERO)
        } else if a < &Integer::ZERO {
            (-Integer::ONE, -a)
        } else {
            (Integer::ONE, a.clone())
        }
    }
}

impl UniqueFactorizationMonoidSignature for IntegerCanonicalStructure {
    type FactoredExponent = NaturalCanonicalStructure;

    fn factorization_exponents(self: &Arc<Self>) -> Arc<Self::FactoredExponent> {
        Natural::structure()
    }

    fn try_is_irreducible(self: &Arc<Self>, a: &Self::Elem) -> Option<bool> {
        Some(Abs::abs(a).is_irreducible())
    }

    fn factorization_pow(self: &Arc<Self>, a: &Self::Elem, k: &Natural) -> Self::Elem {
        self.nat_pow(a, k)
    }
}

impl FactoringMonoidSignature for IntegerCanonicalStructure {
    fn factor_unchecked(self: &Arc<Self>, a: &Self::Elem) -> Factored<Integer, Natural> {
        if a == &Integer::ZERO {
            Factored::Zero
        } else {
            let unit = if a < &Integer::ZERO {
                Integer::from(-1)
            } else {
                Integer::from(1)
            };
            let f = Abs::abs(a).factor();
            Integer::structure()
                .factorizations()
                .new_unit_and_powers_unchecked(
                    unit,
                    f.into_powers()
                        .unwrap()
                        .into_iter()
                        .map(|(p, k)| (Integer::from(p), k))
                        .collect(),
                )
        }
    }
}

impl EuclideanDivisionSignature for IntegerCanonicalStructure {
    fn norm(self: &Arc<Self>, elem: &Self::Elem) -> Option<Natural> {
        if elem == &Integer::ZERO {
            None
        } else {
            Some(Abs::abs(elem))
        }
    }

    fn quorem(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<(Self::Elem, Self::Elem)> {
        if b == &Integer::ZERO {
            None
        } else {
            Some(a.div_mod(b.clone()))
        }
    }
}

impl GreatestCommonDivisorSignature for IntegerCanonicalStructure {
    fn gcd(self: &Arc<Self>, x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        Integer::structure().euclidean_gcd(x.clone(), y.clone())
    }
}

impl BezoutDomainSignature for IntegerCanonicalStructure {
    fn xgcd(
        self: &Arc<Self>,
        x: &Self::Elem,
        y: &Self::Elem,
    ) -> (Self::Elem, Self::Elem, Self::Elem) {
        Integer::euclidean_xgcd(x.clone(), y.clone())
    }
}

impl DedekindDomainSignature for IntegerCanonicalStructure {}

impl CharZeroRingSignature for IntegerCanonicalStructure {
    fn try_to_int(self: &Arc<Self>, x: &Integer) -> Option<Integer> {
        Some(x.clone())
    }
}

impl AlgebraicIntegerRingSignature<RationalCanonicalStructure> for IntegerCanonicalStructure {
    fn anf(self: &Arc<Self>) -> Arc<RationalCanonicalStructure> {
        Rational::structure()
    }

    fn to_anf(self: &Arc<Self>, x: &Integer) -> Rational {
        Rational::from(x)
    }

    fn try_from_anf(self: &Arc<Self>, y: &Rational) -> Option<Integer> {
        Integer::try_from_rat(y)
    }

    fn integral_basis(self: &Arc<Self>) -> Vec<Integer> {
        vec![Integer::ONE]
    }
}

impl ComplexSubsetSignature for IntegerCanonicalStructure {
    fn as_f64_real_and_imaginary_parts(self: &Arc<Self>, z: &Self::Elem) -> (f64, f64) {
        (self.as_f64(z), 0.0)
    }

    fn as_f32_real_and_imaginary_parts(self: &Arc<Self>, z: &Self::Elem) -> (f32, f32) {
        (self.as_f32(z), 0.0)
    }
}

impl RealSubsetSignature for IntegerCanonicalStructure {
    fn as_f64(self: &Arc<Self>, x: &Self::Elem) -> f64 {
        x.into()
    }

    fn as_f32(self: &Arc<Self>, x: &Self::Elem) -> f32 {
        x.into()
    }
}

impl MultiplicativeMonoidSquareOpsSignature for IntegerCanonicalStructure {
    fn sqrt_if_square(self: &Arc<Self>, a: &Integer) -> Option<Integer> {
        a.sqrt_if_square().map(|n: Natural| Integer::from(n))
    }

    fn is_square(self: &Arc<Self>, a: &Integer) -> bool {
        a.is_square()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum IntegerInitialRingGeneratorNeverType {}

impl FreeRingSignature for IntegerCanonicalStructure {
    type Generator = IntegerInitialRingGeneratorNeverType;

    fn free_generators(self: &Arc<Self>) -> HashSet<Self::Generator> {
        HashSet::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_gcd() {
        assert_eq!(
            Integer::euclidean_gcd(Integer::from(0), Integer::from(0)),
            Integer::from(0)
        );

        assert_eq!(
            Integer::euclidean_gcd(Integer::from(12), Integer::from(0)),
            Integer::from(12)
        );

        assert_eq!(
            Integer::euclidean_gcd(Integer::from(0), Integer::from(12)),
            Integer::from(12)
        );

        assert_eq!(
            Integer::euclidean_gcd(Integer::from(12), Integer::from(18)),
            Integer::from(6)
        );

        assert_eq!(
            Integer::gcd_by_factor(&Integer::from(0), &Integer::from(0)),
            Integer::from(0)
        );

        assert_eq!(
            Integer::gcd_by_factor(&Integer::from(12), &Integer::from(0)),
            Integer::from(12)
        );

        assert_eq!(
            Integer::gcd_by_factor(&Integer::from(0), &Integer::from(12)),
            Integer::from(12)
        );

        assert_eq!(
            Integer::gcd_by_factor(&Integer::from(12), &Integer::from(18)),
            Integer::from(6)
        );

        assert_eq!(
            Integer::lcm_by_factor(&Integer::from(12), &Integer::from(18)),
            Some(Integer::from(36))
        );

        assert_eq!(
            Integer::lcm_by_factor(&Integer::from(0), &Integer::from(18)),
            None
        );

        assert_eq!(
            Integer::lcm_by_factor(&Integer::from(12), &Integer::from(0)),
            None
        );

        assert_eq!(
            Integer::lcm_by_factor(&Integer::from(0), &Integer::from(0)),
            None
        );
    }
}
