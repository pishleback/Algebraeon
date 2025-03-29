use super::natural::factorization::factor;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::traits::DivMod;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

pub mod berlekamp_zassenhaus;
pub mod modulo;
pub mod polynomial;
pub mod zimmermann_polys;

impl SemiRingStructure for CannonicalStructure<Integer> {
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

impl RingStructure for CannonicalStructure<Integer> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }
}

impl UnitsStructure for CannonicalStructure<Integer> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }
}

impl IntegralDomainStructure for CannonicalStructure<Integer> {
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

impl OrderedRingStructure for CannonicalStructure<Integer> {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FiniteUnitsStructure for CannonicalStructure<Integer> {
    fn all_units(&self) -> Vec<Self::Set> {
        vec![Integer::ONE, -Integer::ONE]
    }
}

impl FavoriteAssociateStructure for CannonicalStructure<Integer> {
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

impl UniqueFactorizationStructure for CannonicalStructure<Integer> {}

impl FactorableStructure for CannonicalStructure<Integer> {
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
            Some(Factored::new_unchecked(
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

impl EuclideanDivisionStructure for CannonicalStructure<Integer> {
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

impl GreatestCommonDivisorStructure for CannonicalStructure<Integer> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        Integer::structure().euclidean_gcd(x.clone(), y.clone())
    }
}

impl BezoutDomainStructure for CannonicalStructure<Integer> {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        Integer::euclidean_xgcd(x.clone(), y.clone())
    }
}

impl CharZeroStructure for CannonicalStructure<Integer> {}

impl ComplexSubsetStructure for CannonicalStructure<Integer> {}

impl RealSubsetStructure for CannonicalStructure<Integer> {}

impl RealToFloatStructure for CannonicalStructure<Integer> {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        x.into()
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
