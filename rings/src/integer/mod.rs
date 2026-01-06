use super::natural::factorization::NaturalCanonicalFactorizationStructure;
use crate::algebraic_number_field::AlgebraicIntegerRingSignature;
use crate::natural::NaturalFns;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::traits::DivMod;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::collections::HashSet;

pub mod berlekamp_zassenhaus;
pub mod ideal;
pub mod modulo;
pub mod polynomial;
pub mod zimmermann_polys;

impl AdditiveMonoidSignature for IntegerCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }

    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl AdditiveGroupSignature for IntegerCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a - b
    }
}

impl SemiRingSignature for IntegerCanonicalStructure {
    fn one(&self) -> Self::Set {
        Integer::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

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

impl SemiRingUnitsSignature for IntegerCanonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }
}

impl IntegralDomainSignature for IntegerCanonicalStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.quorem(a, b) {
            Some((q, r)) => {
                if r == self.zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }
}

impl OrderedRingSignature for IntegerCanonicalStructure {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FiniteUnitsSignature for IntegerCanonicalStructure {
    fn all_units(&self) -> Vec<Self::Set> {
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

impl UniqueFactorizationDomainSignature for IntegerCanonicalStructure {
    // type FactorOrdering = Self;
    type Factorizations<SelfB: BorrowedStructure<Self>> = FactoredRingElementStructure<Self, SelfB>;

    fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self> {
        FactoredRingElementStructure::new(self)
    }

    fn into_factorizations(self) -> Self::Factorizations<Self> {
        FactoredRingElementStructure::new(self)
    }

    // fn factor_ordering(&self) -> Cow<Self::FactorOrdering> {
    //     Cow::Borrowed(self)
    // }

    fn debug_try_is_irreducible(&self, a: &Self::Set) -> Option<bool> {
        Some(a.abs().is_prime())
    }
}

impl FactorableSignature for IntegerCanonicalStructure {
    fn factor(&self, a: &Self::Set) -> Option<FactoredRingElement<Integer>> {
        if a == &Integer::ZERO {
            None
        } else {
            let unit = if a < &Integer::ZERO {
                Integer::from(-1)
            } else {
                Integer::from(1)
            };
            let f = a.abs().factor().unwrap();
            Some(
                Integer::structure()
                    .factorizations()
                    .from_unit_and_factor_powers_unchecked(
                        unit,
                        Natural::structure()
                            .factorizations()
                            .into_powers(f)
                            .into_iter()
                            .map(|(p, k)| (Integer::from(p), k))
                            .collect(),
                    ),
            )
        }
    }
}

impl EuclideanDivisionSignature for IntegerCanonicalStructure {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if elem == &Integer::ZERO {
            None
        } else {
            Some(elem.abs())
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
    }
}
