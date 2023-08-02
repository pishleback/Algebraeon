use std::collections::HashMap;

use super::{poly::*, ring::*};
use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ZZ {}

impl ComRing for ZZ {
    type ElemT = Integer;

    fn zero() -> Self::ElemT {
        Self::ElemT::from(0)
    }
    fn one() -> Self::ElemT {
        Self::ElemT::from(1)
    }

    fn neg_mut(elem: &mut Self::ElemT) {
        *elem *= Self::ElemT::from(-1)
    }
    fn neg(elem: Self::ElemT) -> Self::ElemT {
        -elem
    }

    fn add_mut(elem: &mut Self::ElemT, x: &Self::ElemT) {
        *elem += x;
    }
    fn add(a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        a + b
    }
    fn add_ref(a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a + b
    }
    fn add_refs(a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a + b
    }

    fn mul_mut(elem: &mut Self::ElemT, x: &Self::ElemT) {
        *elem *= x;
    }
    fn mul(a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        a * b
    }
    fn mul_ref(a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a * b
    }
    fn mul_refs(a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a * b
    }

    fn div(a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        match <Self as EuclideanDomain>::quorem(a, b) {
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

    fn div_lref(a: &Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        match <Self as EuclideanDomain>::quorem_lref(a, b) {
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

    fn div_rref(a: Self::ElemT, b: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        match <Self as EuclideanDomain>::quorem_rref(a, b) {
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

    fn div_refs(a: &Self::ElemT, b: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        match <Self as EuclideanDomain>::quorem_refs(a, b) {
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

impl CharacteristicZero for ZZ {}

impl FiniteUnits for ZZ {
    fn all_units() -> Vec<Self::ElemT> {
        vec![Self::ElemT::from(1), Self::ElemT::from(-1)]
    }
}

impl IntegralDomain for ZZ {}

impl FavoriteAssociate for ZZ {
    fn factor_fav_assoc(elem: Self::ElemT) -> (Self::ElemT, Self::ElemT) {
        if elem == 0 {
            (Self::one(), Self::zero())
        } else if elem < 0 {
            (Integer::from(-1), Self::neg(elem))
        } else {
            (Integer::from(1), elem)
        }
    }
}

impl UniqueFactorizationDomain for ZZ {}

impl UniquelyFactorable for ZZ {
    fn factor(elem: &Self::ElemT) -> Option<Factored<Self::ElemT>> {
        if elem == &0 {
            None
        } else {
            let unit;
            if elem < &0 {
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
                elem.clone(),
                unit,
                factor_nat(elem.unsigned_abs())
                    .into_iter()
                    .map(|(p, k)| (Integer::from(p), k))
                    .collect(),
            ))
        }
    }
}

impl EuclideanDomain for ZZ {
    fn norm(elem: &Self::ElemT) -> Option<Natural> {
        if elem == &Integer::from(0) {
            None
        } else {
            Some(elem.unsigned_abs())
        }
    }

    fn quorem(a: Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if b == Integer::from(0) {
            None
        } else {
            Some(a.div_mod(b.clone()))
        }
    }
}

// impl InterpolatablePolynomials for ZZ {
//     fn interpolate(points: &Vec<(Self::ElemT, Self::ElemT)>) -> Option<Polynomial<Self>> {
//         Polynomial::interpolate_by_lagrange_basis(points)
//     }
// }

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct QQ {}

impl ComRing for QQ {
    type ElemT = Rational;

    fn zero() -> Self::ElemT {
        Self::ElemT::from(0)
    }
    fn one() -> Self::ElemT {
        Self::ElemT::from(1)
    }

    fn neg_mut(elem: &mut Self::ElemT) {
        *elem *= Self::ElemT::from(-1);
    }
    fn neg_ref(elem: &Self::ElemT) -> Self::ElemT {
        -elem
    }
    fn neg(elem: Self::ElemT) -> Self::ElemT {
        -elem
    }

    fn add_mut(elem: &mut Self::ElemT, x: &Self::ElemT) {
        *elem += x;
    }
    fn add(a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        a + b
    }
    fn add_ref(a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a + b
    }
    fn add_refs(a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a + b
    }

    fn mul_mut(elem: &mut Self::ElemT, x: &Self::ElemT) {
        *elem *= x;
    }
    fn mul(a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        a * b
    }
    fn mul_ref(a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a * b
    }
    fn mul_refs(a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        a * b
    }

    fn div(a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        if b == Rational::from(0) {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}
impl IntegralDomain for QQ {}

impl Field for QQ {
    // fn inv(a: Self) -> Result<Self, OppErr> {
    //     if a.numerator_ref() == &Natural::from(0u8) {
    //         Err(OppErr::DivideByZero)
    //     } else {
    //         Ok(a.reciprocal())
    //     }
    // }
}

impl FieldOfFractions for QQ {
    type R = ZZ;
}

// #[derive(Debug, Clone)]
// struct ModularRing {
//     a: Integer,
//     n: Natural,
// }

// impl ToString for ModularRing {
//     fn to_string(&self) -> String {
//         self.a.to_string()
//     }
// }

// impl PartialEq for ModularRing {
//     fn eq(&self, other: &Self) -> bool {
//         if self.n != other.n {
//             panic!()
//         }
//         self.a == other.a
//     }
// }

// impl Eq for ModularRing {

// }

// impl ComRing for ModularRing {
//     fn zero() -> Self {
//         Self {a : Integer::zero(), n}
//     }
//     fn one() -> Self;
//     fn neg_mut(&mut self);
//     fn neg_ref(&self) -> Self {
//         self.clone().neg()
//     }
//     fn neg(mut self) -> Self {
//         self.neg_mut();
//         self
//     }

//     fn add_mut(&mut self, x: &Self);
//     fn add(mut a: Self, b: Self) -> Self {
//         a.add_mut(&b);
//         a
//     }
//     fn add_ref(mut a: Self, b: &Self) -> Self {
//         a.add_mut(b);
//         a
//     }
//     fn add_refs(a: &Self, b: &Self) -> Self {
//         let mut new_a = a.clone();
//         new_a.add_mut(b);
//         new_a
//     }

//     fn mul_mut(&mut self, x: &Self);
//     fn mul(mut a: Self, b: Self) -> Self {
//         a.mul_mut(&b);
//         a
//     }
//     fn mul_ref(mut a: Self, b: &Self) -> Self {
//         a.mul_mut(b);
//         a
//     }
//     fn mul_refs(a: &Self, b: &Self) -> Self {
//         let mut new_a = a.clone();
//         new_a.mul_mut(b);
//         new_a
//     }

//     fn div(a: Self, b: Self) -> Result<Self, RingDivisionError>;
//     fn div_lref(a: &Self, b: Self) -> Result<Self, RingDivisionError> {
//         Self::div(a.clone(), b)
//     }
//     fn div_rref(a: Self, b: &Self) -> Result<Self, RingDivisionError> {
//         Self::div(a, b.clone())
//     }
//     fn div_refs(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
//         Self::div(a.clone(), b.clone())
//     }

//     fn divisible(a: Self, b: Self) -> bool {
//         match Self::div(a, b) {
//             Ok(_q) => true,
//             Err(RingDivisionError::NotDivisible) => false,
//             Err(RingDivisionError::DivideByZero) => false,
//         }
//     }

//     fn sum(elems: Vec<Self>) -> Self {
//         let mut ans = Self::zero();
//         for elem in elems {
//             ans = Self::add(ans, elem);
//         }
//         ans
//     }

//     fn product(elems: Vec<Self>) -> Self {
//         let mut ans = Self::one();
//         for elem in elems {
//             ans = Self::mul(ans, elem);
//         }
//         ans
//     }

//     fn nat_pow(&self, n: &Natural) -> Self {
//         if *n == 0 {
//             Self::one()
//         } else if *n == 1 {
//             self.clone()
//         } else {
//             debug_assert!(*n >= 2);
//             let (q, r) = n.div_rem(Natural::from(2u8));
//             Self::mul(self.nat_pow(&q), self.nat_pow(&(&q + r)))
//         }
//     }

//     fn int_pow(&self, n: &Integer) -> Option<Self> {
//         if *n == 0 {
//             Some(Self::one())
//         } else if self == &Self::zero() {
//             Some(Self::zero())
//         } else if *n > 0 {
//             Some(self.nat_pow(&n.unsigned_abs()))
//         } else {
//             match self.clone().inv() {
//                 Ok(self_inv) => Some(self_inv.nat_pow(&(-n).unsigned_abs())),
//                 Err(RingDivisionError::NotDivisible) => None,
//                 Err(RingDivisionError::DivideByZero) => panic!(),
//             }
//         }
//     }

//     fn from_int(x: &Integer) -> Self {
//         if *x < 0 {
//             Self::from_int(&-x).neg()
//         } else if *x == 0 {
//             Self::zero()
//         } else if *x == 1 {
//             Self::one()
//         } else {
//             let two = Self::add(Self::one(), Self::one());
//             debug_assert!(*x >= 2);
//             let bits: Vec<bool> = x.unsigned_abs().bits().collect();
//             let mut ans = Self::zero();
//             let mut v = Self::one();
//             for i in 0..bits.len() {
//                 if bits[i] {
//                     ans.add_mut(&v);
//                 }
//                 v.mul_mut(&two);
//             }
//             ans
//         }
//     }

//     fn is_unit(self) -> bool {
//         match Self::div(Self::one(), self) {
//             Ok(_inv) => true,
//             Err(RingDivisionError::DivideByZero) => false,
//             Err(RingDivisionError::NotDivisible) => false,
//             // Err(_) => panic!(),
//         }
//     }

//     fn inv(self) -> Result<Self, RingDivisionError> {
//         Self::div(Self::one(), self)
//     }

//     fn inv_ref(a: &Self) -> Result<Self, RingDivisionError> {
//         Self::div_rref(Self::one(), a)
//     }
// }

// struct ModularField {
//     a: Integer,
//     p: Natural,
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int() {
        //happy div
        {
            let a = Integer::from(18);
            let b = Integer::from(6);
            let c = ZZ::div(a, b);
            match c {
                Ok(_) => {}
                Err(_e) => panic!(),
            }
        }

        //sad div
        {
            let a = Integer::from(18);
            let b = Integer::from(7);
            let c = ZZ::div(a, b);
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
            let (q, r) = ZZ::quorem_refs(&a, &b).unwrap();
            assert!(ZZ::norm(&r) < ZZ::norm(&b));
            assert_eq!(a, b * q + r);
        }

        //xgcd
        {
            let a = Integer::from(31);
            let b = Integer::from(57);
            let (g, x, y) = ZZ::xgcd(a.clone(), b.clone());
            assert_eq!(x * a + y * b, g);
        }
    }
}
