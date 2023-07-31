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

pub trait ComRing: Sized + Clone + PartialEq + Eq + Debug + ToString {
    fn zero() -> Self;
    fn one() -> Self;
    fn neg_mut(&mut self);
    fn neg_ref(&self) -> Self {
        self.clone().neg()
    }
    fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    fn add_mut(&mut self, x: &Self);
    fn add(mut a: Self, b: Self) -> Self {
        a.add_mut(&b);
        a
    }
    fn add_ref(mut a: Self, b: &Self) -> Self {
        a.add_mut(b);
        a
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        new_a.add_mut(b);
        new_a
    }

    fn mul_mut(&mut self, x: &Self);
    fn mul(mut a: Self, b: Self) -> Self {
        a.mul_mut(&b);
        a
    }
    fn mul_ref(mut a: Self, b: &Self) -> Self {
        a.mul_mut(b);
        a
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        new_a.mul_mut(b);
        new_a
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError>;
    fn div_lref(a: &Self, b: Self) -> Result<Self, RingDivisionError> {
        Self::div(a.clone(), b)
    }
    fn div_rref(a: Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::div(a, b.clone())
    }
    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::div(a.clone(), b.clone())
    }

    fn divisible(a: Self, b: Self) -> bool {
        match Self::div(a, b) {
            Ok(_q) => true,
            Err(RingDivisionError::NotDivisible) => false,
            Err(RingDivisionError::DivideByZero) => false,
        }
    }

    fn sum(elems: Vec<Self>) -> Self {
        let mut ans = Self::zero();
        for elem in elems {
            ans = Self::add(ans, elem);
        }
        ans
    }

    fn product(elems: Vec<Self>) -> Self {
        let mut ans = Self::one();
        for elem in elems {
            ans = Self::mul(ans, elem);
        }
        ans
    }

    fn pow(&self, n: &Natural) -> Self {
        if *n == 0 {
            Self::one()
        } else if *n == 1 {
            self.clone()
        } else {
            debug_assert!(*n >= 2);
            let (q, r) = n.div_rem(Natural::from(2u8));
            Self::mul(self.pow(&q), self.pow(&(&q + r)))
        }
    }

    fn from_int(x: &Integer) -> Self {
        if *x < 0 {
            Self::from_int(&-x).neg()
        } else if *x == 0 {
            Self::zero()
        } else if *x == 1 {
            Self::one()
        } else {
            let two = Self::add(Self::one(), Self::one());
            debug_assert!(*x >= 2);
            let bits: Vec<bool> = x.unsigned_abs().bits().collect();
            let mut ans = Self::zero();
            let mut v = Self::one();
            for i in 0..bits.len() {
                if bits[i] {
                    ans.add_mut(&v);
                }
                v.mul_mut(&two);
            }
            ans
        }
    }

    fn is_unit(self) -> bool {
        match Self::div(Self::one(), self) {
            Ok(_inv) => true,
            Err(RingDivisionError::DivideByZero) => false,
            Err(RingDivisionError::NotDivisible) => false,
            // Err(_) => panic!(),
        }
    }

    fn inv(self) -> Result<Self, RingDivisionError> {
        Self::div(Self::one(), self)
    }

    fn inv_ref(a: &Self) -> Result<Self, RingDivisionError> {
        Self::div_rref(Self::one(), a)
    }
}

pub trait InfiniteRing: ComRing {
    //generate an infinite sequence of distinct elements
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>>;
}

pub trait CharacteristicZero: ComRing {
    //promise that the integers are distinct in the ring
}

impl<R: CharacteristicZero> InfiniteRing for R {
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>> {
        struct IntegerIterator {
            next: Integer,
        }

        impl Iterator for IntegerIterator {
            type Item = Integer;

            fn next(&mut self) -> Option<Self::Item> {
                let next = self.next.clone();
                if 0 < next {
                    self.next.neg_mut()
                } else {
                    self.next.neg_mut();
                    self.next.add_mut(&Integer::one());
                }
                Some(next)
            }
        }
        Box::new(
            IntegerIterator {
                next: Integer::from(0),
            }
            .map(|n| R::from_int(&n)),
        )
    }
}

pub trait FiniteUnits: ComRing {
    //a commutative ring with finitely many units
    fn all_units() -> Vec<Self>;
}

pub trait IntegralDomain: ComRing {
    //promise that mul(a, b) == 0 implies a == 0 or b == 0
}

pub trait FavoriteAssociate: IntegralDomain {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)

    //it happens that usually the product of favorite associates is another favorite associate. Should this be a requirement?

    fn factor_fav_assoc(self) -> (Self, Self);
    fn factor_fav_assoc_ref(&self) -> (Self, Self) {
        self.clone().factor_fav_assoc()
    }
    fn is_fav_assoc(&self) -> bool {
        let (_u, a) = self.clone().factor_fav_assoc();
        self == &a
    }
}

pub trait GreatestCommonDivisorDomain: FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(x: Self, y: Self) -> Self;
    fn gcd_list(elems: Vec<&Self>) -> Self {
        let mut ans = Self::zero();
        for x in elems {
            ans = Self::gcd(ans, x.clone());
        }
        ans
    }
}

pub trait UniqueFactorizationDomain: FavoriteAssociate + Hash {
    //unique factorizations exist
}

pub trait UniquelyFactorable: UniqueFactorizationDomain {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self) -> Option<Factored<Self>>;
    fn is_irreducible(&self) -> Option<bool> {
        match self.factor() {
            None => None,
            Some(factored) => Some(factored.is_irreducible()),
        }
    }
}

#[derive(Debug)]
pub struct Factored<R: UniqueFactorizationDomain> {
    elem: R,
    unit: R,
    factors: HashMap<R, Natural>,
}

impl<R: UniqueFactorizationDomain> PartialEq for Factored<R> {
    fn eq(&self, other: &Self) -> bool {
        self.unit == other.unit && self.factors == other.factors
    }
}

impl<R: UniqueFactorizationDomain> Eq for Factored<R> {}

impl<R: UniquelyFactorable> Factored<R> {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !self.unit.clone().is_unit() {
            return Err("unit must be a unit");
        }
        let mut prod = self.unit.clone();
        for (p, k) in &self.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if p == &R::zero() {
                return Err("prime factor must not be zero");
            }
            if !p.is_fav_assoc() {
                return Err("prime factor must be their fav assoc");
            }
            if !p.is_irreducible().unwrap() {
                return Err("prime factor must not be reducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                prod.mul_mut(p);
                i += Natural::from(1u8);
            }
        }
        if self.elem != prod {
            return Err("product is incorrect");
        }
        Ok(())
    }
}

impl<R: UniqueFactorizationDomain> ToString for Factored<R> {
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

impl<R: UniqueFactorizationDomain> Factored<R> {
    pub fn one() -> Self {
        Self {
            elem: R::one(),
            unit: R::one(),
            factors: HashMap::new(),
        }
    }

    pub fn new_unchecked(elem: R, unit: R, factors: HashMap<R, Natural>) -> Self {
        Self {
            elem,
            unit,
            factors,
        }
    }

    pub fn new_unit_unchecked(unit: R) -> Self {
        Self {
            elem: unit.clone(),
            unit,
            factors: HashMap::new(),
        }
    }

    pub fn new_irreducible_unchecked(elem: R) -> Self {
        let (unit, assoc) = elem.factor_fav_assoc();
        Self {
            elem: assoc.clone(),
            unit,
            factors: HashMap::from([(assoc, Natural::from(1u8))]),
        }
    }

    pub fn unit(&self) -> &R {
        &self.unit
    }

    pub fn factors(&self) -> &HashMap<R, Natural> {
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

    pub fn mul(a: Self, b: Self) -> Self {
        let mut mul_factors = a.factors;
        for (p, k) in b.factors {
            *mul_factors.entry(p.clone()).or_insert(Natural::from(0u8)) += k;
        }
        Self::new_unchecked(R::mul(a.elem, b.elem), R::mul(a.unit, b.unit), mul_factors)
    }

    pub fn divisors<'a>(&self) -> Box<dyn Iterator<Item = R> + 'a>
    where
        R: 'a,
    {
        if self.factors.len() == 0 {
            Box::new(vec![R::one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in &self.factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = R::one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = R::mul_ref(p_pow, &p);
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(|prime_power_factors| R::product(prime_power_factors).clone()),
            )
        }
    }

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
    fn xgcd(x: Self, y: Self) -> (Self, Self, Self);
    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        match elems.len() {
            0 => (Self::zero(), vec![]),
            1 => {
                let (unit, assoc) = elems[0].factor_fav_assoc_ref();
                (assoc, vec![unit.inv().unwrap()])
            }
            2 => {
                let (g, x, y) = Self::xgcd(elems[0].clone(), elems[1].clone());
                (g, vec![x, y])
            }
            n => {
                let k = n / 2;
                let (g1, coeffs1) = Self::xgcd_list((0..k).map(|i| elems[i]).collect());
                let (g2, coeffs2) = Self::xgcd_list((k..n).map(|i| elems[i]).collect());
                let (g, x, y) = Self::xgcd(g1, g2);
                let mut coeffs = vec![];
                for c in coeffs1 {
                    coeffs.push(Self::mul_refs(&x, &c));
                }
                for c in coeffs2 {
                    coeffs.push(Self::mul_refs(&y, &c));
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
    fn norm(&self) -> Option<Natural>;
    fn quorem(a: Self, b: Self) -> Option<(Self, Self)>;
    fn quorem_lref(a: &Self, b: Self) -> Option<(Self, Self)> {
        Self::quorem(a.clone(), b)
    }
    fn quorem_rref(a: Self, b: &Self) -> Option<(Self, Self)> {
        Self::quorem(a, b.clone())
    }
    fn quorem_refs(a: &Self, b: &Self) -> Option<(Self, Self)> {
        Self::quorem(a.clone(), b.clone())
    }

    fn quo(a: Self, b: Self) -> Option<Self> {
        match Self::quorem(a, b) {
            Some((q, _r)) => Some(q),
            None => None,
        }
    }
    fn quo_lref(a: &Self, b: Self) -> Option<Self> {
        Self::quo(a.clone(), b)
    }
    fn quo_rref(a: Self, b: &Self) -> Option<Self> {
        Self::quo(a, b.clone())
    }
    fn quo_refs(a: &Self, b: &Self) -> Option<Self> {
        Self::quo(a.clone(), b.clone())
    }

    fn rem(a: Self, b: Self) -> Option<Self> {
        match Self::quorem(a, b) {
            Some((_q, r)) => Some(r),
            None => None,
        }
    }
    fn rem_lref(a: &Self, b: Self) -> Option<Self> {
        Self::rem(a.clone(), b)
    }
    fn rem_rref(a: Self, b: &Self) -> Option<Self> {
        Self::rem(a, b.clone())
    }
    fn rem_refs(a: &Self, b: &Self) -> Option<Self> {
        Self::rem(a.clone(), b.clone())
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> GreatestCommonDivisorDomain for R {
    fn gcd(mut x: Self, mut y: Self) -> Self {
        //Euclidean algorithm
        while y != Self::zero() {
            let r = Self::rem_rref(x, &y).unwrap();
            (x, y) = (y, r)
        }
        let (_unit, assoc) = x.factor_fav_assoc();
        assoc
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> PrincipalIdealDomain for R {
    fn xgcd(mut x: Self, mut y: Self) -> (Self, Self, Self) {
        let mut pa = Self::one();
        let mut a = Self::zero();
        let mut pb = Self::zero();
        let mut b = Self::one();

        while y != Self::zero() {
            let (q, r) = Self::quorem_rref(x, &y).unwrap();
            let new_a = Self::add(pa, Self::mul_refs(&q, &a).neg());
            (a, pa) = (new_a, a);
            let new_b = Self::add(pb, Self::mul_ref(q, &b).neg());
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = x.factor_fav_assoc();
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(unit.clone().is_unit());
        (
            ass_x,
            R::div_rref(pa, &unit).unwrap(),
            R::div(pb, unit).unwrap(),
        )
    }
}

pub trait Field: IntegralDomain {
    //promise that a/b always works, except unless b=0.
    //in other words, a/b must not return not divisible
}

impl<F: Field> FavoriteAssociate for F {
    fn factor_fav_assoc(self) -> (Self, Self) {
        (self, Self::one())
    }
}

impl<F: Field> EuclideanDomain for F {
    fn norm(&self) -> Option<Natural> {
        if self == &Self::zero() {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(a: Self, b: Self) -> Option<(Self, Self)> {
        if b == Self::zero() {
            None
        } else {
            Some((Self::div(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_lref(a: &Self, b: Self) -> Option<(Self, Self)> {
        if b == Self::zero() {
            None
        } else {
            Some((Self::div_lref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_rref(a: Self, b: &Self) -> Option<(Self, Self)> {
        if b == &Self::zero() {
            None
        } else {
            Some((Self::div_rref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_refs(a: &Self, b: &Self) -> Option<(Self, Self)> {
        if b == &Self::zero() {
            None
        } else {
            Some((Self::div_refs(a, b).unwrap(), Self::zero()))
        }
    }

    fn quo(a: Self, b: Self) -> Option<Self> {
        if b == Self::zero() {
            None
        } else {
            Some(Self::div(a, b).unwrap())
        }
    }

    fn quo_lref(a: &Self, b: Self) -> Option<Self> {
        if b == Self::zero() {
            None
        } else {
            Some(Self::div_lref(a, b).unwrap())
        }
    }

    fn quo_rref(a: Self, b: &Self) -> Option<Self> {
        if b == &Self::zero() {
            None
        } else {
            Some(Self::div_rref(a, b).unwrap())
        }
    }

    fn quo_refs(a: &Self, b: &Self) -> Option<Self> {
        if b == &Self::zero() {
            None
        } else {
            Some(Self::div_refs(a, b).unwrap())
        }
    }

    fn rem(_a: Self, b: Self) -> Option<Self> {
        if b == Self::zero() {
            None
        } else {
            Some(Self::zero())
        }
    }

    fn rem_lref(_a: &Self, b: Self) -> Option<Self> {
        if b == Self::zero() {
            None
        } else {
            Some(Self::zero())
        }
    }

    fn rem_rref(_a: Self, b: &Self) -> Option<Self> {
        if b == &Self::zero() {
            None
        } else {
            Some(Self::zero())
        }
    }

    fn rem_refs(_a: &Self, b: &Self) -> Option<Self> {
        if b == &Self::zero() {
            None
        } else {
            Some(Self::zero())
        }
    }
}

impl<F: Field + Hash> UniqueFactorizationDomain for F {}

impl<F: Field + Hash> UniquelyFactorable for F {
    fn factor(&self) -> Option<Factored<Self>> {
        if self == &Self::zero() {
            None
        } else {
            Some(Factored::new_unchecked(
                self.clone(),
                self.clone(),
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
    use super::*;

    #[test]
    fn factorization_invariants() {
        let f = Factored::new_unchecked(
            Integer::from(-12),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(Integer::from(1), Integer::from(1), HashMap::from([]));
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(
            Integer::from(-12),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ]),
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "can't have a power of zero"
        );

        let f = Factored::new_unchecked(
            Integer::from(-13),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        assert_eq!(f.check_invariants().is_ok(), false, "product is incorrect");

        let f = Factored::new_unchecked(
            Integer::from(12),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit is wrong");

        let f = Factored::new_unchecked(
            Integer::from(12),
            Integer::from(3),
            HashMap::from([(Integer::from(2), Natural::from(2u8))]),
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

        let f = Factored::new_unchecked(
            Integer::from(0),
            Integer::from(1),
            HashMap::from([
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must not be zero"
        );

        let f = Factored::new_unchecked(
            Integer::from(-12),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be prime"
        );

        let f = Factored::new_unchecked(
            Integer::from(-12),
            Integer::from(-1),
            HashMap::from([
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be fav assoc"
        );
    }

    #[test]
    fn test_xgcd_list() {
        use malachite_q::Rational;
        let a = Rational::from(7);
        let (g, taps) = Rational::xgcd_list(vec![&a]);
        assert_eq!(g, Rational::one());
        assert_eq!(taps.len(), 1);
        assert_eq!(g, &taps[0] * a);
    }

    #[test]
    fn test_divisors() {
        for a in 1u8..25 {
            let b = Integer::from(a);
            let fs = b.factor().unwrap();
            assert_eq!(
                fs.count_divisors().unwrap(),
                Natural::from(fs.divisors().collect::<Vec<Integer>>().len())
            );
        }
    }
}
