#![allow(dead_code)]

use std::{borrow::Borrow, collections::HashMap, fmt::Debug, hash::Hash};

use malachite_base::num::{
    arithmetic::traits::{DivRem, UnsignedAbs},
    logic::traits::BitIterable,
};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

#[derive(Debug)]
pub enum RingDivisionError {
    DivideByZero,
    NotDivisible,
}

pub trait ComRing: Sized + Clone + Debug + PartialEq + Eq {
    //todo: remove sized here
    type ElemT: Sized + Clone + Debug;

    fn to_string(&self, elem: &Self::ElemT) -> String;

    fn equal(&self, a: &Self::ElemT, b: &Self::ElemT) -> bool;

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

    fn are_associate(&self, a: &Self::ElemT, b: &Self::ElemT) -> bool {
        if self.equal(a, &self.zero()) && self.equal(b, &self.zero()) {
            true
        } else {
            self.div_refs(a, b).is_ok() && self.div_refs(b, a).is_ok()
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
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![elem.clone()];
            while pows.len() < bits.len() {
                pows.push(self.mul_refs(&pows.last().unwrap(), &pows.last().unwrap()));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = self.one();
            for i in 0..count {
                if bits[i] {
                    self.mul_mut(&mut ans, &pows[i]);
                }
            }
            ans
        }
    }

    fn int_pow(&self, elem: &Self::ElemT, n: &Integer) -> Option<Self::ElemT> {
        println!("{:?} {:?}", elem, n);
        if *n == 0 {
            Some(self.one())
        } else if self.equal(elem, &self.zero()) {
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

    fn from_rat(&self, x: &Rational) -> Result<Self::ElemT, RingDivisionError> {
        self.div(
            self.from_int(&super::nzq::QQ.numerator(x)),
            self.from_int(&super::nzq::QQ.denominator(x)),
        )
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
    fn generate_distinct_elements<'a>(&'a self) -> Box<dyn Iterator<Item = Self::ElemT> + 'a>;
}

pub trait CharacteristicZero: ComRing {
    //promise that the integers are distinct in the ring
}

impl<R: CharacteristicZero> InfiniteRing for R {
    fn generate_distinct_elements<'a>(&'a self) -> Box<dyn Iterator<Item = Self::ElemT> + 'a> {
        struct IntegerIterator {
            next: Integer,
        }

        impl Iterator for IntegerIterator {
            type Item = Integer;

            fn next(&mut self) -> Option<Self::Item> {
                let next = self.next.clone();
                if 0 < next {
                    self.next = -self.next.clone();
                } else {
                    self.next = Integer::from(1) - self.next.clone();
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
    fn fav_assoc(&self, elem: Self::ElemT) -> Self::ElemT {
        self.factor_fav_assoc(elem).1
    }
    fn factor_fav_assoc_ref(&self, elem: &Self::ElemT) -> (Self::ElemT, Self::ElemT) {
        self.factor_fav_assoc(elem.clone())
    }
    fn is_fav_assoc(&self, elem: &Self::ElemT) -> bool {
        let (_u, a) = self.factor_fav_assoc(elem.clone());
        self.equal(elem, &a)
    }
}

pub trait GreatestCommonDivisorDomain: FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(&self, x: Self::ElemT, y: Self::ElemT) -> Self::ElemT;
    fn gcd_list<BorrowElemT: Borrow<Self::ElemT>>(&self, elems: Vec<BorrowElemT>) -> Self::ElemT {
        let mut ans = self.zero();
        for x in elems {
            ans = self.gcd(ans, x.borrow().clone());
        }
        ans
    }
    fn lcm(&self, x: Self::ElemT, y: Self::ElemT) -> Self::ElemT {
        if self.equal(&x, &self.zero()) && self.equal(&y, &self.zero()) {
            self.zero()
        } else {
            let g = self.gcd(x.clone(), y.clone());
            self.div(self.mul(x, y), g).unwrap()
        }
    }

    fn lcm_list<BorrowElemT: Borrow<Self::ElemT>>(&self, elems: Vec<BorrowElemT>) -> Self::ElemT {
        let mut ans = self.one();
        for x in elems {
            ans = self.lcm(ans, x.borrow().clone());
        }
        ans
    }
}

#[derive(Debug)]
pub struct Factored<ElemT: Clone> {
    unit: ElemT,
    factors: Vec<(ElemT, Natural)>,
}

impl<ElemT: Clone + Debug> Factored<ElemT> {
    fn check_invariants<Ring: UniqueFactorizationDomain<ElemT = ElemT>>(
        &self,
        ring: &Ring,
    ) -> Result<(), &'static str> {
        if !ring.is_unit(self.unit.clone()) {
            return Err("unit must be a unit");
        }
        for (p, k) in &self.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if ring.equal(p, &ring.zero()) {
                return Err("prime factor must not be zero");
            }
            if !ring.is_fav_assoc(p) {
                return Err("prime factor must be their fav assoc");
            }
            if !ring.is_irreducible(p).unwrap() {
                return Err("prime factor must not be reducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                i += Natural::from(1u8);
            }
        }
        Ok(())
    }

    pub fn expand<Ring: ComRing<ElemT = ElemT>>(&self, ring: &Ring) -> ElemT {
        let mut ans = self.unit.clone();
        for (p, k) in &self.factors {
            ring.mul_mut(&mut ans, &ring.nat_pow(p, k));
        }
        ans
    }

    pub fn new_unchecked(unit: ElemT, factors: Vec<(ElemT, Natural)>) -> Self {
        Self { unit, factors }
    }

    pub fn equal<Ring: ComRing<ElemT = ElemT>>(
        ring: &Ring,
        a: &Factored<ElemT>,
        b: &Factored<ElemT>,
    ) -> bool {
        ring.equal(&a.expand(ring), &b.expand(ring))
    }

    pub fn unit(&self) -> &ElemT {
        &self.unit
    }

    pub fn factors(&self) -> &Vec<(ElemT, Natural)> {
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

    fn mul_by_unchecked<Ring: ComRing<ElemT = ElemT>>(
        &mut self,
        ring: &Ring,
        p: ElemT,
        k: Natural,
    ) {
        for (q, t) in &mut self.factors {
            match (ring.div_refs(&p, q), ring.div_refs(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if ring.is_unit(u) && ring.is_unit(v.clone()) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        ring.mul_mut(&mut self.unit, &ring.nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        self.factors.push((p, k));
    }

    pub fn mul<Ring: ComRing<ElemT = ElemT>>(
        ring: &Ring,
        mut a: Factored<ElemT>,
        b: Factored<ElemT>,
    ) -> Factored<ElemT> {
        ring.mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            a.mul_by_unchecked(ring, p, k)
        }
        a
    }

    pub fn factored_one<Ring: ComRing<ElemT = ElemT>>(ring: &Ring) -> Self {
        Factored {
            unit: ring.one(),
            factors: vec![],
        }
    }

    pub fn factored_irreducible_unchecked<Ring: FavoriteAssociate<ElemT = ElemT>>(
        ring: &Ring,
        elem: ElemT,
    ) -> Self {
        let (unit, assoc) = ring.factor_fav_assoc(elem);
        Factored {
            unit,
            factors: vec![(assoc, Natural::from(1u8))],
        }
    }

    pub fn factored_unit_unchecked<Ring: ComRing<ElemT = ElemT>>(
        _ring: &Ring,
        unit: ElemT,
    ) -> Self {
        Factored {
            unit,
            factors: vec![],
        }
    }

    pub fn to_string<Ring: ComRing<ElemT = ElemT>>(&self, ring: &Ring) -> String {
        let mut s = ring.to_string(&self.unit);
        for (f, k) in &self.factors {
            s += " * (";
            s += &ring.to_string(f);
            s += ")";
            if k != &Natural::from(1u8) {
                s += "^";
                s += &k.to_string();
            }
        }
        s
    }
}

// impl<ElemT: Clone + PartialEq> PartialEq for Factored<ElemT> {
//     fn eq(&self, other: &Self) -> bool {
//         self.unit == other.unit && self.factors == other.factors
//     }
// }

// impl<ElemT: Clone + Eq> Eq for Factored<ElemT> {}

pub trait UniqueFactorizationDomain: FavoriteAssociate {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self, elem: &Self::ElemT) -> Option<Factored<Self::ElemT>>;

    fn is_irreducible(&self, elem: &Self::ElemT) -> Option<bool> {
        match self.factor(elem) {
            None => None,
            Some(factored) => Some(factored.is_irreducible()),
        }
    }

    fn divisors<'a>(
        &'a self,
        factored: &Factored<Self::ElemT>,
    ) -> Box<dyn Iterator<Item = Self::ElemT> + 'a> {
        if factored.factors.len() == 0 {
            Box::new(vec![self.one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in &factored.factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = self.one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = self.mul_ref(p_pow, &p);
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(|prime_power_factors| self.product(prime_power_factors).clone()),
            )
        }
    }

    fn count_divisors(&self, factored: &Factored<Self::ElemT>) -> Option<Natural> {
        let mut count = Natural::from(1u8);
        for (_p, k) in &factored.factors {
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
        while !self.equal(&y, &self.zero()) {
            let r = self.rem_rref(x, &y).unwrap();
            (x, y) = (y, r)
        }
        let (_unit, assoc) = self.factor_fav_assoc(x);
        assoc
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> PrincipalIdealDomain for R {
    fn xgcd(
        &self,
        mut x: Self::ElemT,
        mut y: Self::ElemT,
    ) -> (Self::ElemT, Self::ElemT, Self::ElemT) {
        let mut pa = self.one();
        let mut a = self.zero();
        let mut pb = self.zero();
        let mut b = self.one();

        while !self.equal(&y, &self.zero()) {
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

trait QuotientRing {
    type Ring: ComRing;
}

#[derive(Debug, Clone)]
pub struct EuclideanQuotient<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain>
{
    ed: ED,
    n: ED::ElemT,
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> PartialEq
    for EuclideanQuotient<IS_FIELD, ED>
{
    fn eq(&self, other: &Self) -> bool {
        if self.ed == other.ed {
            let ans = self.ed.equal(&self.n, &other.n);
            debug_assert_eq!(ans, other.ed.equal(&self.n, &other.n));
            ans
        } else {
            false
        }
    }
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> Eq
    for EuclideanQuotient<IS_FIELD, ED>
{
}

impl<ED: EuclideanDomain + UniqueFactorizationDomain> EuclideanQuotient<false, ED> {
    pub fn new_ring(ed: ED, n: ED::ElemT) -> Self {
        Self { ed, n }
    }
}

impl<ED: EuclideanDomain + UniqueFactorizationDomain> EuclideanQuotient<true, ED> {
    pub fn new_field(ed: ED, n: ED::ElemT) -> Self {
        Self { ed, n }
    }
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain>
    EuclideanQuotient<IS_FIELD, ED>
{
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if self.ed.equal(&self.n, &self.ed.zero()) {
            return Err("cant quotient by zero");
        }

        if IS_FIELD {
            if !self.ed.is_irreducible(&self.n).unwrap() {
                return Err("marked as field but element is not irreducible");
            }
        }

        Ok(())
    }
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> ComRing
    for EuclideanQuotient<IS_FIELD, ED>
{
    type ElemT = ED::ElemT;

    fn to_string(&self, elem: &Self::ElemT) -> String {
        self.ed.to_string(elem)
    }

    fn equal(&self, a: &Self::ElemT, b: &Self::ElemT) -> bool {
        self.ed.equal(
            &self
                .ed
                .rem_rref(self.ed.add_ref(self.ed.neg_ref(a), b), &self.n)
                .unwrap(),
            &self.zero(),
        )
    }

    fn zero(&self) -> Self::ElemT {
        self.ed.zero()
    }

    fn one(&self) -> Self::ElemT {
        self.ed.one()
    }

    fn neg_mut(&self, elem: &mut Self::ElemT) {
        *elem = self.ed.rem_rref(self.ed.neg_ref(elem), &self.n).unwrap();
    }

    fn add_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT) {
        *elem = self
            .ed
            .rem_rref(self.ed.add_refs(elem, offset), &self.n)
            .unwrap();
    }

    fn mul_mut(&self, elem: &mut Self::ElemT, mul: &Self::ElemT) {
        *elem = self
            .ed
            .rem_rref(self.ed.mul_refs(elem, mul), &self.n)
            .unwrap();
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        todo!();
    }
}

//IdealQuotient
//PrimeQuotient
//MaximalQuotient

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
        if self.equal(elem, &self.zero()) {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(&self, a: Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some((self.div(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some((self.div_lref(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some((self.div_rref(a, b).unwrap(), self.zero()))
        }
    }

    fn quorem_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some((self.div_refs(a, b).unwrap(), self.zero()))
        }
    }

    fn quo(&self, a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some(self.div(a, b).unwrap())
        }
    }

    fn quo_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some(self.div_lref(a, b).unwrap())
        }
    }

    fn quo_rref(&self, a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some(self.div_rref(a, b).unwrap())
        }
    }

    fn quo_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some(self.div_refs(a, b).unwrap())
        }
    }

    fn rem(&self, _a: Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_lref(&self, _a: &Self::ElemT, b: Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(&b, &self.zero()) {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_rref(&self, _a: Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some(self.zero())
        }
    }

    fn rem_refs(&self, _a: &Self::ElemT, b: &Self::ElemT) -> Option<Self::ElemT> {
        if self.equal(b, &self.zero()) {
            None
        } else {
            Some(self.zero())
        }
    }
}

impl<F: Field + Hash> UniqueFactorizationDomain for F
where
    Self::ElemT: Hash,
{
    fn factor(&self, elem: &Self::ElemT) -> Option<Factored<Self::ElemT>> {
        if self.equal(elem, &self.zero()) {
            None
        } else {
            Some(Factored::new_unchecked(elem.clone(), vec![]))
        }
    }
}

pub trait FieldOfFractions: Field {
    type R: IntegralDomain;

    fn base_ring(&self) -> &Self::R;
    fn from_base_ring(&self, elem: <Self::R as ComRing>::ElemT) -> Self::ElemT;
    fn numerator(&self, elem: &Self::ElemT) -> <Self::R as ComRing>::ElemT;
    fn denominator(&self, elem: &Self::ElemT) -> <Self::R as ComRing>::ElemT;
    fn as_base_ring(&self, elem: Self::ElemT) -> Option<<Self::R as ComRing>::ElemT> {
        if self
            .base_ring()
            .equal(&self.denominator(&elem), &self.base_ring().one())
        {
            Some(self.numerator(&elem))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::nzq::*;
    use super::*;

    #[test]
    fn factorization_invariants() {
        let f = Factored::new_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        f.check_invariants(&ZZ).unwrap();

        let f = Factored::new_unchecked(Integer::from(1), vec![]);
        f.check_invariants(&ZZ).unwrap();

        let f = Factored::new_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ],
        );
        assert_eq!(
            f.check_invariants(&ZZ).is_ok(),
            false,
            "can't have a power of zero"
        );

        let f = Factored::new_unchecked(
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert_eq!(
            f.check_invariants(&ZZ).is_ok(),
            false,
            "unit should be a unit"
        );

        let f = Factored::new_unchecked(
            Integer::from(1),
            vec![
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants(&ZZ).is_ok(),
            false,
            "prime factors must not be zero"
        );

        let f = Factored::new_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants(&ZZ).is_ok(),
            false,
            "prime factors must be prime"
        );

        let f = Factored::new_unchecked(
            Integer::from(-1),
            vec![
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants(&ZZ).is_ok(),
            false,
            "prime factors must be fav assoc"
        );
    }

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
