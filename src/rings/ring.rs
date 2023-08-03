#![allow(dead_code)]

use std::{collections::HashMap, fmt::Debug, hash::Hash};

use malachite_base::num::{
    arithmetic::traits::{DivRem, UnsignedAbs},
    logic::traits::BitIterable,
};
use malachite_nz::{integer::Integer, natural::Natural};

#[derive(Debug)]
pub enum RingDivisionError {
    DivideByZero,
    NotDivisible,
}

pub trait ComRing: Sized + Clone + PartialEq + Eq + Hash + Debug {
    //todo: remove sized here
    type ElemT: Sized + Clone + PartialEq + Eq + Hash + Debug;

    fn to_string(&self, elem : &Self::ElemT) -> String;

    fn zero(&self) -> Self::ElemT;
    fn one(&self) -> Self::ElemT;

    fn neg_mut(&self, elem: &mut Self::ElemT);
    fn neg_ref(&self, elem: &Self::ElemT) -> Self::ElemT {
        self.neg(elem.clone())
    }
    fn neg(&self, mut elem: Self::ElemT) -> Self::ElemT {
        self.neg_mut(&mut elem);
        elem
    }

    fn add_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT);
    fn add(&self, mut a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        self.add_mut(&mut a, &b);
        a
    }
    fn add_ref(&self, mut a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        self.add_mut(&mut a, b);
        a
    }
    fn add_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        let mut new_a = a.clone();
        self.add_mut(&mut new_a, b);
        new_a
    }

    fn mul_mut(&self, elem: &mut Self::ElemT, mul: &Self::ElemT);
    fn mul(&self, mut a: Self::ElemT, b: Self::ElemT) -> Self::ElemT {
        self.mul_mut(&mut a, &b);
        a
    }
    fn mul_ref(&self, mut a: Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        self.mul_mut(&mut a, b);
        a
    }
    fn mul_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        let mut new_a = a.clone();
        self.mul_mut(&mut new_a, b);
        new_a
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError>;
    fn div_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div(a.clone(), b)
    }
    fn div_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div(a, b.clone())
    }
    fn div_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div(a.clone(), b.clone())
    }

    fn divisible(&self, a: Self::ElemT, b: Self::ElemT) -> bool {
        match self.div(a, b) {
            Ok(_q) => true,
            Err(RingDivisionError::NotDivisible) => false,
            Err(RingDivisionError::DivideByZero) => false,
        }
    }

    fn sum(&self, elems: Vec<Self::ElemT>) -> Self::ElemT {
        let mut ans = self.zero();
        for elem in elems {
            ans = self.add(ans, elem);
        }
        ans
    }

    fn product(&self, elems: Vec<Self::ElemT>) -> Self::ElemT {
        let mut ans = self.one();
        for elem in elems {
            ans = self.mul(ans, elem);
        }
        ans
    }

    fn nat_pow(&self, elem: &Self::ElemT, n: &Natural) -> Self::ElemT {
        if *n == 0 {
            self.one()
        } else if *n == 1 {
            elem.clone()
        } else {
            debug_assert!(*n >= 2);
            let (q, r) = n.div_rem(Natural::from(2u8));
            self.mul(self.nat_pow(elem, &q), self.nat_pow(elem, &(&q + r)))
        }
    }

    fn int_pow(&self, elem: &Self::ElemT, n: &Integer) -> Option<Self::ElemT> {
        if *n == 0 {
            Some(self.one())
        } else if elem == &self.zero() {
            Some(self.zero())
        } else if *n > 0 {
            Some(self.nat_pow(elem, &n.unsigned_abs()))
        } else {
            match self.inv(elem.clone()) {
                Ok(elem_inv) => Some(self.nat_pow(&elem_inv, &(-n).unsigned_abs())),
                Err(RingDivisionError::NotDivisible) => None,
                Err(RingDivisionError::DivideByZero) => panic!(),
            }
        }
    }

    fn from_int(&self, x: &Integer) -> Self::ElemT {
        if *x < 0 {
            self.neg(self.from_int(&-x))
        } else if *x == 0 {
            self.zero()
        } else if *x == 1 {
            self.one()
        } else {
            let two = self.add(self.one(), self.one());
            debug_assert!(*x >= 2);
            let bits: Vec<bool> = x.unsigned_abs().bits().collect();
            let mut ans = self.zero();
            let mut v = self.one();
            for i in 0..bits.len() {
                if bits[i] {
                    self.add_mut(&mut ans, &v);
                }
                self.mul_mut(&mut v, &two);
            }
            ans
        }
    }

    fn is_unit(&self, elem: Self::ElemT) -> bool {
        match self.div(self.one(), elem) {
            Ok(_inv) => true,
            Err(RingDivisionError::DivideByZero) => false,
            Err(RingDivisionError::NotDivisible) => false,
            // Err(_) => panic!(),
        }
    }

    fn inv(&self, elem: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div(self.one(), elem)
    }

    fn inv_ref(&self, a: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div_rref(self.one(), a)
    }
}

pub trait InfiniteRing: ComRing {
    //generate an infinite sequence of distinct elements
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::ElemT>>;
}

pub trait CharacteristicZero: ComRing {
    //promise that the integers are distinct in the ring
}

impl<R: CharacteristicZero> InfiniteRing for R {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::ElemT>> {
        struct IntegerIterator {
            next: Integer,
        }

        impl Iterator for IntegerIterator {
            type Item = Integer;

            fn next(&mut self) -> Option<Self::Item> {
                let mut next = self.next.clone();
                if 0 < next {
                    next = -next;
                } else {
                    next = Integer::from(1) - next;
                }
                Some(next)
            }
        }
        Box::new(
            IntegerIterator {
                next: Integer::from(0),
            }
            .map(|n| self.from_int(&n)),
        )
    }
}

pub trait FiniteUnits: ComRing {
    //a commutative ring with finitely many units
    fn all_units(&self) -> Vec<Self::ElemT>;
}

pub trait IntegralDomain: ComRing {
    //promise that mul(a, b) == 0 implies a == 0 or b == 0
}

pub trait FavoriteAssociate: IntegralDomain {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)

    //it happens that usually the product of favorite associates is another favorite associate. Should this be a requirement?

    fn factor_fav_assoc(&self, elem: Self::ElemT) -> (Self::ElemT, Self::ElemT);
    fn factor_fav_assoc_ref(&self, elem: &Self::ElemT) -> (Self::ElemT, Self::ElemT) {
        self.factor_fav_assoc(elem.clone())
    }
    fn is_fav_assoc(&self, elem: &Self::ElemT) -> bool {
        let (_u, a) = self.factor_fav_assoc(elem.clone());
        elem == &a
    }
}

pub trait GreatestCommonDivisorDomain: FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(&self, x: Self::ElemT, y: Self::ElemT) -> Self::ElemT;
    fn gcd_list(&self, elems: Vec<&Self::ElemT>) -> Self::ElemT {
        let mut ans = self.zero();
        for x in elems {
            ans = self.gcd(ans, x.clone());
        }
        ans
    }
}

pub trait UniqueFactorizationDomain: FavoriteAssociate {
    //unique factorizations exist
}

pub trait UniquelyFactorable: UniqueFactorizationDomain {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self, elem: &Self::ElemT) -> Option<Factored<Self::ElemT>>;

    fn is_irreducible(&self, elem: &Self::ElemT) -> Option<bool> {
        match self.factor(elem) {
            None => None,
            Some(factored) => Some(factored.is_irreducible()),
        }
    }

    fn one_factored(&self) -> Factored<Self::ElemT> {
        Factored {
            elem: self.one(),
            unit: self.one(),
            factors: HashMap::new(),
        }
    }

    fn irreducible_factored_unchecked(&self, elem: Self::ElemT) -> Factored<Self::ElemT> {
        let (unit, assoc) = self.factor_fav_assoc(elem);
        Factored {
            elem: assoc.clone(),
            unit,
            factors: HashMap::from([(assoc, Natural::from(1u8))]),
        }
    }

    fn check_invariants(&self, f: Factored<Self::ElemT>) -> Result<(), &'static str> {
        if !self.is_unit(f.unit.clone()) {
            return Err("unit must be a unit");
        }
        let mut prod = f.unit.clone();
        for (p, k) in &f.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if p == &self.zero() {
                return Err("prime factor must not be zero");
            }
            if !self.is_fav_assoc(p) {
                return Err("prime factor must be their fav assoc");
            }
            if !self.is_irreducible(p).unwrap() {
                return Err("prime factor must not be reducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                self.mul_mut(&mut prod, p);
                i += Natural::from(1u8);
            }
        }
        if f.elem != prod {
            return Err("product is incorrect");
        }
        Ok(())
    }
}

#[derive(Debug)]
pub struct Factored<ElemT: PartialEq + Eq + Hash + Clone> {
    elem: ElemT,
    unit: ElemT,
    factors: HashMap<ElemT, Natural>,
}

impl<ElemT: PartialEq + Eq + Hash + Clone> PartialEq for Factored<ElemT> {
    fn eq(&self, other: &Self) -> bool {
        self.unit == other.unit && self.factors == other.factors
    }
}

impl<ElemT: PartialEq + Eq + Hash + Clone> Eq for Factored<ElemT> {}

impl<ElemT: PartialEq + Eq + Hash + Clone + ToString> ToString for Factored<ElemT> {
    fn to_string(&self) -> String {
        let mut s = self.unit.to_string();
        for (f, k) in &self.factors {
            s += " * (";
            s += &f.to_string();
            s += ")";
            if k != &Natural::from(1u8) {
                s += "^";
                s += &k.to_string();
            }
        }
        s
    }
}

impl<ElemT: PartialEq + Eq + Hash + Clone> Factored<ElemT> {
    pub fn new_unchecked(elem: ElemT, unit: ElemT, factors: HashMap<ElemT, Natural>) -> Self {
        Self {
            elem,
            unit,
            factors,
        }
    }

    pub fn new_unit_unchecked(unit: ElemT) -> Self {
        Self {
            elem: unit.clone(),
            unit,
            factors: HashMap::new(),
        }
    }

    pub fn unit(&self) -> &ElemT {
        &self.unit
    }

    pub fn factors(&self) -> &HashMap<ElemT, Natural> {
        &self.factors
    }

    pub fn is_irreducible(&self) -> bool {
        if self.factors.len() == 1 {
            let (_p, k) = self.factors.iter().next().unwrap();
            if k == &1 {
                return true;
            }
        }
        false
    }

    // pub fn mul(a: Self, b: Self) -> Self {
    //     let mut mul_factors = a.factors;
    //     for (p, k) in b.factors {
    //         *mul_factors.entry(p.clone()).or_insert(Natural::from(0u8)) += k;
    //     }
    //     Self::new_unchecked(R::mul(a.elem, b.elem), R::mul(a.unit, b.unit), mul_factors)
    // }

    // pub fn divisors<'a>(&self) -> Box<dyn Iterator<Item = R::ElemT> + 'a>
    // where
    //     R: 'a,
    // {
    //     if self.factors.len() == 0 {
    //         Box::new(vec![R::one()].into_iter())
    //     } else {
    //         let mut factor_powers = vec![];
    //         for (p, k) in &self.factors {
    //             let j = factor_powers.len();
    //             factor_powers.push(vec![]);
    //             let mut p_pow = R::one();
    //             let mut i = Natural::from(0u8);
    //             while &i <= k {
    //                 factor_powers[j].push(p_pow.clone());
    //                 p_pow = R::mul_ref(p_pow, &p);
    //                 i += Natural::from(1u8);
    //             }
    //         }

    //         Box::new(
    //             itertools::Itertools::multi_cartesian_product(
    //                 factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
    //             )
    //             .map(|prime_power_factors| R::product(prime_power_factors).clone()),
    //         )
    //     }
    // }

    pub fn count_divisors(&self) -> Option<Natural> {
        let mut count = Natural::from(1u8);
        for (_p, k) in &self.factors {
            count *= k + Natural::from(1u8);
        }
        Some(count)
    }
}

pub trait PrincipalIdealDomain: GreatestCommonDivisorDomain {
    //any gcds should be the standard associate representative
    fn xgcd(&self, x: Self::ElemT, y: Self::ElemT) -> (Self::ElemT, Self::ElemT, Self::ElemT);
    fn xgcd_list(&self, elems: Vec<&Self::ElemT>) -> (Self::ElemT, Vec<Self::ElemT>) {
        match elems.len() {
            0 => (self.zero(), vec![]),
            1 => {
                let (unit, assoc) = self.factor_fav_assoc_ref(elems[0]);
                (assoc, vec![self.inv(unit).unwrap()])
            }
            2 => {
                let (g, x, y) = self.xgcd(elems[0].clone(), elems[1].clone());
                (g, vec![x, y])
            }
            n => {
                let k = n / 2;
                let (g1, coeffs1) = self.xgcd_list((0..k).map(|i| elems[i]).collect());
                let (g2, coeffs2) = self.xgcd_list((k..n).map(|i| elems[i]).collect());
                let (g, x, y) = self.xgcd(g1, g2);
                let mut coeffs = vec![];
                for c in coeffs1 {
                    coeffs.push(self.mul_refs(&x, &c));
                }
                for c in coeffs2 {
                    coeffs.push(self.mul_refs(&y, &c));
                }
                (g, coeffs)
            }
        }

        // if all(elem == 0 for elem in elems):
        //     return cls.int(0), [cls.int(0) for elem in elems]
        // elems = list(elems)
        // assert len(elems) >= 1
        // if len(elems) == 1:
        //     return elems[0], [cls.int(1)]
        // elif len(elems) == 2:
        //     g, x, y = cls.xgcd(elems[0], elems[1])
        //     return g, [x, y]
        // else:
        //     n = len(elems) // 2
        //     g1, coeffs1 = cls.xgcd_list(elems[:n])
        //     g2, coeffs2 = cls.xgcd_list(elems[n:])
        //     g, x, y = cls.xgcd(g1, g2)
        //     return g, [x * c for c in coeffs1] + [y * c for c in coeffs2]
    }
}

pub trait EuclideanDomain: IntegralDomain {
    //should return None for 0, and Some(norm) for everything else
    fn norm(&self, elem: &Self::ElemT) -> Option<Natural>;
    fn quorem(&self, a: Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)>;
    fn quorem_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        self.quorem(a.clone(), b)
    }
    fn quorem_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        self.quorem(a, b.clone())
    }
    fn quorem_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        self.quorem(a.clone(), b.clone())
    }

    fn quo(&self, a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        match self.quorem(a, b) {
            Some((q, _r)) => Some(q),
            None => None,
        }
    }
    fn quo_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        self.quo(a.clone(), b)
    }
    fn quo_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        self.quo(a, b.clone())
    }
    fn quo_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        self.quo(a.clone(), b.clone())
    }

    fn rem(&self, a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        match self.quorem(a, b) {
            Some((_q, r)) => Some(r),
            None => None,
        }
    }
    fn rem_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        self.rem(a.clone(), b)
    }
    fn rem_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        self.rem(a, b.clone())
    }
    fn rem_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        self.rem(a.clone(), b.clone())
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> GreatestCommonDivisorDomain for R {
    fn gcd(&self, mut x: Self::ElemT, mut y: Self::ElemT) -> Self::ElemT {
        //Euclidean algorithm
        while y != self.zero() {
            let r = self.rem_rref(x, &y).unwrap();
            (x, y) = (y, r)
        }
        let (_unit, assoc) = self.factor_fav_assoc(x);
        assoc
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> PrincipalIdealDomain for R {
    fn xgcd(&self, mut x: Self::ElemT, mut y: Self::ElemT) -> (Self::ElemT, Self::ElemT, Self::ElemT) {
        let mut pa = self.one();
        let mut a = self.zero();
        let mut pb = self.zero();
        let mut b = self.one();

        while y != self.zero() {
            let (q, r) = self.quorem_rref(x, &y).unwrap();
            let new_a = self.add(pa, self.neg(self.mul_refs(&q, &a)));
            (a, pa) = (new_a, a);
            let new_b = self.add(pb, self.neg(self.mul_ref(q, &b)));
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = self.factor_fav_assoc(x);
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(self.is_unit(unit.clone()));
        (
            ass_x,
            self.div_rref(pa, &unit).unwrap(),
            self.div(pb, unit).unwrap(),
        )
    }
}

pub trait Field: IntegralDomain {
    //promise that a/b always works, except unless b=0.
    //in other words, a/b must not return not divisible
}

impl<F: Field> FavoriteAssociate for F {
    fn factor_fav_assoc(&self, elem: Self::ElemT) -> (Self::ElemT, Self::ElemT) {
        (elem, self.one())
    }
}

impl<F: Field> EuclideanDomain for F {
    fn norm(&self, elem: &Self::ElemT) -> Option<Natural> {
        if elem == &self.zero() {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(&self, a: Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if b == self.zero() {
            None
        } else {
            Some((self.div(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if b == self.zero() {
            None
        } else {
            Some((self.div_lref(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if b == &self.zero() {
            None
        } else {
            Some((self.div_rref(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if b == &self.zero() {
            None
        } else {
            Some((self.div_refs(a, b).unwrap(), self.zero()))
        }
    }

    fn quo(&self, a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if b == self.zero() {
            None
        } else {
            Some(self.div(a, b).unwrap())
        }
    }

    fn quo_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if b == self.zero() {
            None
        } else {
            Some(self.div_lref(a, b).unwrap())
        }
    }

    fn quo_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if b == &self.zero() {
            None
        } else {
            Some(self.div_rref(a, b).unwrap())
        }
    }

    fn quo_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if b == &self.zero() {
            None
        } else {
            Some(self.div_refs(a, b).unwrap())
        }
    }

    fn rem(&self, _a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if b == self.zero() {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_lref(&self, _a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if b == self.zero() {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_rref(&self, _a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if b == &self.zero() {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_refs(&self, _a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if b == &self.zero() {
            None
        } else {
            Some(self.zero())
        }
    }
}

impl<F: Field + Hash> UniqueFactorizationDomain for F {}

impl<F: Field + Hash> UniquelyFactorable for F
where
    Self::ElemT: Hash,
{
    fn factor(&self, elem: &Self::ElemT) -> Option<Factored<Self::ElemT>> {
        if elem == &self.zero() {
            None
        } else {
            Some(Factored::new_unchecked(
                elem.clone(),
                elem.clone(),
                HashMap::new(),
            ))
        }
    }
}

pub trait FieldOfFractions: Field {
    type R: IntegralDomain;
}

#[cfg(test)]
mod tests {
    use super::super::nzq::*;
    use super::*;

    // #[test]
    // fn factorization_invariants() {
    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(-12),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     f.check_invariants().unwrap();

    //     let f =
    //         Factored::<ZZ>::new_unchecked(Integer::from(1), Integer::from(1), HashMap::from([]));
    //     f.check_invariants().unwrap();

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(-12),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //             (Integer::from(5), Natural::from(0u8)),
    //         ]),
    //     );
    //     assert_eq!(
    //         f.check_invariants().is_ok(),
    //         false,
    //         "can't have a power of zero"
    //     );

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(-13),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     assert_eq!(f.check_invariants().is_ok(), false, "product is incorrect");

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(12),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     assert_eq!(f.check_invariants().is_ok(), false, "unit is wrong");

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(12),
    //         Integer::from(3),
    //         HashMap::from([(Integer::from(2), Natural::from(2u8))]),
    //     );
    //     assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(0),
    //         Integer::from(1),
    //         HashMap::from([
    //             (Integer::from(0), Natural::from(1u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     assert_eq!(
    //         f.check_invariants().is_ok(),
    //         false,
    //         "prime factors must not be zero"
    //     );

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(-12),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(4), Natural::from(1u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     assert_eq!(
    //         f.check_invariants().is_ok(),
    //         false,
    //         "prime factors must be prime"
    //     );

    //     let f = Factored::<ZZ>::new_unchecked(
    //         Integer::from(-12),
    //         Integer::from(-1),
    //         HashMap::from([
    //             (Integer::from(-2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ]),
    //     );
    //     assert_eq!(
    //         f.check_invariants().is_ok(),
    //         false,
    //         "prime factors must be fav assoc"
    //     );
    // }

    #[test]
    fn test_xgcd_list() {
        use malachite_q::Rational;
        let a = Rational::from(7);
        let (g, taps) = QQ.xgcd_list(vec![&a]);
        assert_eq!(g, QQ.one());
        assert_eq!(taps.len(), 1);
        assert_eq!(g, &taps[0] * a);
    }

    // #[test]
    // fn test_divisors() {
    //     for a in 1u8..25 {
    //         let b = Integer::from(a);
    //         let fs = ZZ::factor(&b).unwrap();
    //         assert_eq!(
    //             fs.count_divisors().unwrap(),
    //             Natural::from(fs.divisors().collect::<Vec<Integer>>().len())
    //         );
    //     }
    // }
}
