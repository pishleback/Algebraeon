#[allow(dead_code)]
use std::collections::HashMap;

use super::poly::*;
use super::ring::*;
use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

impl ComRing for Integer {
    fn to_string(&self) -> String {
        ToString::to_string(&self)
    }

    // fn equal( a: &Self, b: &Self) -> bool {
    //     a == b
    // }

    fn zero() -> Self {
        Self::from(0)
    }
    fn one() -> Self {
        Self::from(1)
    }

    fn neg_mut(&mut self) {
        *self *= Self::from(-1)
    }
    fn neg(self) -> Self {
        -self
    }

    fn add_mut(elem: &mut Self, x: &Self) {
        *elem += x;
    }
    fn add(a: Self, b: Self) -> Self {
        a + b
    }
    fn add_ref(a: Self, b: &Self) -> Self {
        a + b
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        a + b
    }

    fn mul_mut(elem: &mut Self, x: &Self) {
        *elem *= x;
    }
    fn mul(a: Self, b: Self) -> Self {
        a * b
    }
    fn mul_ref(a: Self, b: &Self) -> Self {
        a * b
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        a * b
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
        match Self::quorem(a, b) {
            Some((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }

    fn div_lref(a: &Self, b: Self) -> Result<Self, RingDivisionError> {
        match Self::quorem_lref(a, b) {
            Some((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }

    fn div_rref(a: Self, b: &Self) -> Result<Self, RingDivisionError> {
        match Self::quorem_rref(a, b) {
            Some((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }

    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        match Self::quorem_refs(a, b) {
            Some((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }
}

impl CharacteristicZero for Integer {}

impl FiniteUnits for Integer {
    fn all_units() -> Vec<Self> {
        vec![Self::from(1), Self::from(-1)]
    }
}

impl IntegralDomain for Integer {}

impl FavoriteAssociate for Integer {
    fn factor_fav_assoc(self) -> (Self, Self) {
        if self == 0 {
            (Self::one(), Self::zero())
        } else if self < 0 {
            (Integer::from(-1), Self::neg(self))
        } else {
            (Integer::from(1), self)
        }
    }
}

impl UniqueFactorizationDomain for Integer {
    fn factor(&self) -> Option<Factored<Self>> {
        if self == &0 {
            None
        } else {
            let unit;
            if self < &0 {
                unit = Integer::from(-1);
            } else {
                unit = Integer::from(1);
            }

            fn factor_nat(mut n: Natural) -> HashMap<Natural, Natural> {
                //TODO: more efficient implementations
                assert_ne!(n, 0);
                let mut fs = HashMap::new();
                let mut p = Natural::from(2u8);
                while n > 1 && p <= n {
                    while &n % &p == 0 {
                        *fs.entry(p.clone()).or_insert(Natural::from(0u8)) += Natural::from(1u8);
                        n /= &p;
                    }
                    p += Natural::from(1u8);
                }
                fs
            }

            Some(Factored::new_unchecked(
                unit,
                factor_nat(self.unsigned_abs())
                    .into_iter()
                    .map(|(p, k)| (Integer::from(p), k))
                    .collect(),
            ))
        }
    }
}

impl EuclideanDomain for Integer {
    fn norm(elem: &Self) -> Option<Natural> {
        if elem == &Integer::from(0) {
            None
        } else {
            Some(elem.unsigned_abs())
        }
    }

    fn quorem(a: Self, b: Self) -> Option<(Self, Self)> {
        if b == Integer::from(0) {
            None
        } else {
            Some(a.div_mod(b.clone()))
        }
    }
}

impl Real for Integer {
    fn as_f64(x: Self) -> f64 {
        if x < 0 {
            -Self::as_f64(-x)
        } else {
            let limbs = x.into_twos_complement_limbs_asc();
            let mut flt = 0.0;
            for (i, k) in limbs.into_iter().enumerate() {
                flt += (k as f64) * (2.0 as f64).powf(i as f64 * 64.0);
            }
            flt
        }
    }
}

impl DenseReal for Rational {
    fn from_f64_approx(x: f64) -> Self {
        let mut x = Rational::from_sci_string_simplest(x.to_string().as_str()).unwrap();
        malachite_q::arithmetic::traits::ApproximateAssign::approximate_assign(
            &mut x,
            &Natural::from(100usize),
        );
        x
    }
}

impl UniqueFactorizationDomain for Polynomial<Integer> {
    fn factor(&self) -> Option<Factored<Self>> {
        //TODO: use a more efficient algorithm: reduce mod large prime and lift, knapsack alg using LLL to pick factor subsets efficiently
        self.factorize_by_kroneckers_method()
    }
}

// impl InterpolatablePolynomials for ZZ {
//     fn interpolate(points: &Vec<(Self, Self)>) -> Option<Polynomial<Self>> {
//         Polynomial::interpolate_by_lagrange_basis(points)
//     }
// }

impl ComRing for Rational {
    fn to_string(&self) -> String {
        ToString::to_string(&self)
    }

    // fn equal( a: &Self, b: &Self) -> bool {
    //     a == b
    // }

    fn zero() -> Self {
        Self::from(0)
    }
    fn one() -> Self {
        Self::from(1)
    }

    fn neg_mut(&mut self) {
        *self *= Self::from(-1);
    }
    fn neg_ref(&self) -> Self {
        -self
    }
    fn neg(self) -> Self {
        -self
    }

    fn add_mut(elem: &mut Self, x: &Self) {
        *elem += x;
    }
    fn add(a: Self, b: Self) -> Self {
        a + b
    }
    fn add_ref(a: Self, b: &Self) -> Self {
        a + b
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        a + b
    }

    fn mul_mut(elem: &mut Self, x: &Self) {
        *elem *= x;
    }
    fn mul(a: Self, b: Self) -> Self {
        a * b
    }
    fn mul_ref(a: Self, b: &Self) -> Self {
        a * b
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        a * b
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
        if b == Rational::from(0) {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}
impl IntegralDomain for Rational {}

impl Field for Rational {}

impl FieldOfFractions for Rational {
    type R = Integer;

    fn numerator(elem: &Self) -> Self::R {
        if elem >= &0 {
            Integer::from(elem.numerator_ref())
        } else {
            -Integer::from(elem.numerator_ref())
        }
    }

    fn denominator(elem: &Self) -> Self::R {
        Integer::from(elem.denominator_ref())
    }

    fn from_base_ring(elem: Self::R) -> Self {
        Rational::from(elem)
    }
}

// impl FiniteUnits for EuclideanQuotient<true, IntegerRing> {
//     fn all_units(&self) -> Vec<Self> {
//         let mut units = vec![];
//         let mut u = self.one();
//         while u < self.get_n() {
//             units.push(u.clone());
//             u += self.one();
//         }
//         units
//     }
// }

pub struct NaturalPrimeGenerator {
    n: Natural,
    primes: Vec<Natural>,
}

impl NaturalPrimeGenerator {
    pub fn new() -> Self {
        Self {
            n: Natural::from(2u8),
            primes: vec![],
        }
    }
}

impl Iterator for NaturalPrimeGenerator {
    type Item = Natural;

    fn next(&mut self) -> Option<Self::Item> {
        'next_loop: loop {
            //todo: only check primes up to sqrt n
            for p in &self.primes {
                if &self.n % p == 0 {
                    self.n += Natural::from(1u8);
                    continue 'next_loop;
                }
            }
            let next_p = self.n.clone();
            self.n += Natural::from(1u8);
            self.primes.push(next_p.clone());
            return Some(next_p);
        }
    }
}

pub fn factorial(n: Natural) -> Natural {
    let mut k = Natural::from(1u8);
    let mut i = Natural::from(1u8);
    while i <= n {
        k *= &i;
        i += Natural::from(1u8);
    }
    k
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int() {
        //happy div
        {
            let a = Integer::from(18);
            let b = Integer::from(6);
            let c = Integer::div(a, b);
            match c {
                Ok(_) => {}
                Err(_e) => panic!(),
            }
        }

        //sad div
        {
            let a = Integer::from(18);
            let b = Integer::from(7);
            let c = Integer::div(a, b);
            match c {
                Ok(_) => panic!(),
                Err(e) => match e {
                    RingDivisionError::DivideByZero => panic!(),
                    RingDivisionError::NotDivisible => {}
                },
            }
        }

        //euclidean div
        {
            let a = Integer::from(18);
            let b = Integer::from(7);
            let (q, r) = Integer::quorem_refs(&a, &b).unwrap();
            assert!(Integer::norm(&r) < Integer::norm(&b));
            assert_eq!(a, b * q + r);
        }

        //xgcd
        {
            let a = Integer::from(31);
            let b = Integer::from(57);
            let (g, x, y) = Integer::xgcd(a.clone(), b.clone());
            assert_eq!(x * a + y * b, g);
        }
    }

    #[test]
    fn test_rational_numerator_and_denominator() {
        let x = Rational::from_signeds(-22, 7);
        let (n, d) = (Rational::numerator(&x), Rational::denominator(&x));
        assert_eq!(n, Integer::from(-22));
        assert_eq!(d, Integer::from(7));

        let x = Rational::from_signeds(22, -7);
        let (n, d) = (Rational::numerator(&x), Rational::denominator(&x));
        assert_eq!(n, Integer::from(-22));
        assert_eq!(d, Integer::from(7));

        let x = Rational::from_signeds(22, 7);
        let (n, d) = (Rational::numerator(&x), Rational::denominator(&x));
        assert_eq!(n, Integer::from(22));
        assert_eq!(d, Integer::from(7));

        let x = Rational::from_signeds(-22, -7);
        let (n, d) = (Rational::numerator(&x), Rational::denominator(&x));
        assert_eq!(n, Integer::from(22));
        assert_eq!(d, Integer::from(7));

        let x = Rational::from_signeds(0, 1);
        let (n, d) = (Rational::numerator(&x), Rational::denominator(&x));
        assert_eq!(n, Integer::from(0));
        assert_eq!(d, Integer::from(1));
    }

    // #[test]
    // fn test_factor_int() {}

    #[test]
    fn test_factorial() {
        debug_assert_eq!(factorial(Natural::from(0u8)), Natural::from(1u8));
        debug_assert_eq!(factorial(Natural::from(1u8)), Natural::from(1u8));
        debug_assert_eq!(factorial(Natural::from(2u8)), Natural::from(2u8));
        debug_assert_eq!(factorial(Natural::from(3u8)), Natural::from(6u8));
        debug_assert_eq!(factorial(Natural::from(4u8)), Natural::from(24u8));
        debug_assert_eq!(factorial(Natural::from(5u8)), Natural::from(120u8));
        debug_assert_eq!(factorial(Natural::from(6u8)), Natural::from(720u16));
    }
}
