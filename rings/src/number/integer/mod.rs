use factor::factor_by_try_primes;
use itertools::Itertools;
use malachite_base::num::arithmetic::traits::DivMod;
use malachite_base::num::arithmetic::traits::UnsignedAbs;
use malachite_base::num::basic::traits::One;
use malachite_base::num::basic::traits::Two;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use crate::number::natural::*;
use crate::ring_structure::cannonical::*;
use crate::ring_structure::factorization::*;
use crate::ring_structure::structure::*;

use algebraeon_structure::*;

pub mod modulo;
pub mod polynomial;

impl RingStructure for CannonicalStructure<Integer> {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }

    fn one(&self) -> Self::Set {
        Integer::ONE
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
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
        if a == &0 {
            (Integer::ONE, Integer::ZERO)
        } else if a < &0 {
            (-Integer::ONE, -a)
        } else {
            (Integer::ONE, a.clone())
        }
    }
}

impl UniqueFactorizationStructure for CannonicalStructure<Integer> {
    fn factor(&self, a: &Self::Set) -> Option<Factored<Self>> {
        if a == &0 {
            None
        } else {
            let unit;
            if a < &0 {
                unit = Integer::from(-1);
            } else {
                unit = Integer::from(1);
            }

            Some(Factored::new_unchecked(
                self.clone().into(),
                unit,
                factor_by_try_primes(a.unsigned_abs())
                    .unwrap()
                    .into_powers()
                    .into_iter()
                    .map(|(p, k)| (Integer::from(p), k))
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
            Some(elem.unsigned_abs())
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
        if x < &0 {
            -self.as_f64(&-x)
        } else {
            let limbs = x.clone().into_twos_complement_limbs_asc();
            let mut flt = 0.0;
            for (i, k) in limbs.into_iter().enumerate() {
                flt += (k as f64) * (2.0 as f64).powf(i as f64 * 64.0);
            }
            flt
        }
    }
}
