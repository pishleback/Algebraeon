use crate::algebraic_number_field::AlgebraicIntegerRingSignature;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::traits::DivMod;
use algebraeon_nzq::*;
use algebraeon_sets::structure::BorrowedStructure;
use algebraeon_sets::structure::CountableSetSignature;
use algebraeon_sets::structure::EqSignature;
use algebraeon_sets::structure::FiniteSetSignature;
use algebraeon_sets::structure::MetaType;
use std::collections::HashSet;

pub mod berlekamp_zassenhaus;
pub mod ideal;
pub mod modulo;
pub mod polynomial;
pub mod zimmermann_polys;

impl RinglikeSpecializationSignature for IntegerCanonicalStructure {
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Set = Self::Set> + RingSignature> {
        Some(self.clone())
    }

    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Set = Self::Set> + CharZeroRingSignature> {
        Some(self.clone())
    }
}

impl ZeroSignature for IntegerCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }
}

impl AdditionSignature for IntegerCanonicalStructure {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
}

impl CancellativeAdditionSignature for IntegerCanonicalStructure {
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl TryNegateSignature for IntegerCanonicalStructure {
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }
}

impl AdditiveMonoidSignature for IntegerCanonicalStructure {}

impl AdditiveGroupSignature for IntegerCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a - b
    }
}

impl OneSignature for IntegerCanonicalStructure {
    fn one(&self) -> Self::Set {
        Integer::ONE
    }
}

impl MultiplicationSignature for IntegerCanonicalStructure {
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
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
    fn is_reduced(&self) -> Result<bool, String> {
        Ok(true)
    }
}

impl CharacteristicSignature for IntegerCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl TryReciprocalSignature for IntegerCanonicalStructure {
    fn try_reciprocal(&self, a: &Self::Set) -> Option<Self::Set> {
        self.try_divide(&self.one(), a)
    }
}

impl CancellativeMultiplicationSignature for IntegerCanonicalStructure {
    fn try_divide(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
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

impl OrderedRingSignature for IntegerCanonicalStructure {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<IntegerCanonicalStructure, B>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        self.list_all_elements().into_iter()
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<IntegerCanonicalStructure, B>
{
    fn list_all_elements(&self) -> Vec<Self::Set> {
        vec![Integer::ONE, -Integer::ONE]
    }
}

impl FavoriteAssociateSignature for IntegerCanonicalStructure {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
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

    fn factorization_exponents(&self) -> &Self::FactoredExponent {
        Natural::structure_ref()
    }

    fn into_factorization_exponents(self) -> Self::FactoredExponent {
        Natural::structure()
    }

    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool> {
        Some(Abs::abs(a).is_irreducible())
    }

    fn factorization_pow(&self, a: &Self::Set, k: &Natural) -> Self::Set {
        self.nat_pow(a, k)
    }
}

impl FactoringMonoidSignature for IntegerCanonicalStructure {
    fn factor_unchecked(&self, a: &Self::Set) -> Factored<Integer, Natural> {
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
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if elem == &Integer::ZERO {
            None
        } else {
            Some(Abs::abs(elem))
        }
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        if b == &Integer::ZERO {
            None
        } else {
            Some(a.div_mod(b.clone()))
        }
    }
}

impl GreatestCommonDivisorSignature for IntegerCanonicalStructure {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        Integer::structure().euclidean_gcd(x.clone(), y.clone())
    }
}

impl BezoutDomainSignature for IntegerCanonicalStructure {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        Integer::euclidean_xgcd(x.clone(), y.clone())
    }
}

impl DedekindDomainSignature for IntegerCanonicalStructure {}

impl CharZeroRingSignature for IntegerCanonicalStructure {
    fn try_to_int(&self, x: &Integer) -> Option<Integer> {
        Some(x.clone())
    }
}

impl AlgebraicIntegerRingSignature<RationalCanonicalStructure> for IntegerCanonicalStructure {
    fn anf(&self) -> &RationalCanonicalStructure {
        Rational::structure_ref()
    }

    fn to_anf(&self, x: &Integer) -> Rational {
        Rational::from(x)
    }

    fn try_from_anf(&self, y: &Rational) -> Option<Integer> {
        Integer::try_from_rat(y)
    }

    fn integral_basis(&self) -> Vec<Integer> {
        vec![Integer::ONE]
    }
}

impl ComplexSubsetSignature for IntegerCanonicalStructure {
    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
        (self.as_f64(z), 0.0)
    }

    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32) {
        (self.as_f32(z), 0.0)
    }
}

impl RealSubsetSignature for IntegerCanonicalStructure {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        x.into()
    }

    fn as_f32(&self, x: &Self::Set) -> f32 {
        x.into()
    }
}

impl MultiplicativeMonoidSquareOpsSignature for IntegerCanonicalStructure {
    fn sqrt_if_square(&self, a: &Integer) -> Option<Integer> {
        a.sqrt_if_square().map(|n: Natural| Integer::from(n))
    }

    fn is_square(&self, a: &Integer) -> bool {
        a.is_square()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum IntegerInitialRingGeneratorNeverType {}

impl FreeRingSignature for IntegerCanonicalStructure {
    type Generator = IntegerInitialRingGeneratorNeverType;

    fn free_generators(&self) -> HashSet<Self::Generator> {
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
