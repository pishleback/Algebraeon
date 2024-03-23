use std::collections::HashMap;
use std::rc::Rc;
use std::thread::sleep;
use std::time::Duration;

use malachite_base::num::arithmetic::traits::DivMod;
use malachite_base::num::arithmetic::traits::UnsignedAbs;
use malachite_base::num::basic::traits::One;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use crate::ring_structure::quotient::QuotientStructure;

use super::super::super::structure::*;
use super::super::polynomial::polynomial::*;
use super::super::ring_structure::cannonical::*;
use super::super::ring_structure::factorization::*;
use super::super::ring_structure::structure::*;

use super::natural::*;

impl StructuredType for Integer {
    type Structure = CannonicalStructure<Self>;

    fn structure() -> Rc<Self::Structure> {
        Self::Structure::new().into()
    }
}

impl EqualityStructure for CannonicalStructure<Integer> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

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
    fn factor(
        &self,
        a: &Self::Set,
    ) -> Option<crate::ring_structure::factorization::Factored<Self>> {
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
                factor_nat(a.unsigned_abs())
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

impl RealStructure for CannonicalStructure<Integer> {
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

impl GreatestCommonDivisorStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
        //TODO: use a better algorithm here
        self.factorize_by_kroneckers_method(p)
    }
}

impl Polynomial<Integer> {
    fn find_factor_primitive_sqfree_by_zassenhaus_algorithm(
        self,
    ) -> Option<(Polynomial<Integer>, Polynomial<Integer>)> {
        let f = self;
        let f_deg = f.degree().unwrap();
        debug_assert_ne!(f_deg, 0);
        println!("zassenhaus: {}", f);
        if f_deg == 1 {
            None
        } else {
            let prime_gen = NaturalPrimeGenerator::new();
            for p in prime_gen.take(10) {
                let mod_p = QuotientStructure::new_field(Integer::structure(), Integer::from(&p));
                let poly_mod_p = PolynomialStructure::new(mod_p.into());

                println!("f mod {} = {}", p, poly_mod_p.elem_to_string(&f));
                println!("{}", poly_mod_p.factor(&f).unwrap());

                // let f_mod_p = f.apply_map_ref(|c| {
                //     UniversalEuclideanQuotient::<true, _>::new(c.clone(), Integer::from(&p))
                // });

                // println!("{}", f_mod_p);

                sleep(Duration::from_millis(100));
            }

            todo!()
        }
    }

    pub fn factorize_by_zassenhaus_algorithm(
        self,
    ) -> Option<Factored<PolynomialStructure<CannonicalStructure<Integer>>>> {
        if self == Self::zero() {
            None
        } else {
            Some(
                Self::structure().factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(
                    self,
                    &|f| {
                        factorize_by_find_factor(&Self::structure(), f, &|f| {
                            Self::find_factor_primitive_sqfree_by_zassenhaus_algorithm(f)
                        })
                    },
                ),
            )
        }
    }
}

impl FiniteUnitsStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn all_units(&self) -> Vec<Self::Set> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < self.modulus().unsigned_abs() {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl FiniteFieldStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.modulus().unsigned_abs(), Natural::ONE)
    }
}

impl UniqueFactorizationStructure
    for PolynomialStructure<QuotientStructure<CannonicalStructure<Integer>, true>>
{
    fn factor(
        &self,
        a: &Self::Set,
    ) -> Option<crate::ring_structure::factorization::Factored<Self>> {
        self.factorize_by_berlekamps_algorithm(a.clone())
    }
}

// #[cfg(test)]
// mod tests {
//     use crate::{rings::ring_structure::elements::RingElement, ComRing};

//     use super::*;

//     #[test]
//     fn run() {
//         let a = Integer::from(-4);
//         let b = Integer::from(5);

//         let rs = Integer::ring_structure();
//         println!("{:?}", rs.add(&a, &b));

//         println!("{:?}", a.factor_fav_assoc());
//     }
// }
