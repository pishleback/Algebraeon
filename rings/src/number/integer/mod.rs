use std::collections::HashSet;

use super::natural::factorization::factor;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::traits::DivMod;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

pub mod berlekamp_zassenhaus;
pub mod ideal;
pub mod modulo;
pub mod polynomial;
pub mod zimmermann_polys;

impl SemiRingStructure for IntegerCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }

    fn one(&self) -> Self::Set {
        Integer::ONE
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl RingStructure for IntegerCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }
}

impl UnitsStructure for IntegerCanonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }
}

impl IntegralDomainStructure for IntegerCanonicalStructure {
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

impl OrderedRingStructure for IntegerCanonicalStructure {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FiniteUnitsStructure for IntegerCanonicalStructure {
    fn all_units(&self) -> Vec<Self::Set> {
        vec![Integer::ONE, -Integer::ONE]
    }
}

impl FavoriteAssociateStructure for IntegerCanonicalStructure {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        if a == &Integer::ZERO {
            (Integer::ONE, Integer::ZERO)
        } else if a < &Integer::ZERO {
            (-Integer::ONE, -a)
        } else {
            (Integer::ONE, a.clone())
        }
    }
}

impl UniqueFactorizationStructure for IntegerCanonicalStructure {
    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool> {
        Some(self.is_irreducible(a))
    }
}

impl FactorableStructure for IntegerCanonicalStructure {
    fn factor(&self, a: &Self::Set) -> Option<Factored<Self>> {
        if a == &Integer::ZERO {
            None
        } else {
            let unit;
            if a < &Integer::ZERO {
                unit = Integer::from(-1);
            } else {
                unit = Integer::from(1);
            }
            let f = factor(a.abs()).unwrap();
            Some(Factored::from_unit_and_factor_powers_unchecked(
                self.clone().into(),
                unit,
                f.into_powers()
                    .into_iter()
                    .map(|(p, k)| (Integer::from(p), Natural::from(k)))
                    .collect(),
            ))
        }
    }
}

impl EuclideanDivisionStructure for IntegerCanonicalStructure {
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

impl GreatestCommonDivisorStructure for IntegerCanonicalStructure {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        Integer::structure().euclidean_gcd(x.clone(), y.clone())
    }
}

impl BezoutDomainStructure for IntegerCanonicalStructure {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        Integer::euclidean_xgcd(x.clone(), y.clone())
    }
}

impl CharZeroRingStructure for IntegerCanonicalStructure {
    fn try_to_int(&self, x: &Integer) -> Option<Integer> {
        Some(x.clone())
    }
}

impl ComplexSubsetStructure for IntegerCanonicalStructure {}

impl RealSubsetStructure for IntegerCanonicalStructure {}

impl RealToFloatStructure for IntegerCanonicalStructure {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        x.into()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum IntegerInitialRingGeneratorNeverType {}

impl FreeRingStructure for IntegerCanonicalStructure {
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
