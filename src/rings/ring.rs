use std::{collections::HashMap, fmt::Debug, hash::Hash};

use malachite_base::num::{
    arithmetic::traits::{DivRem, UnsignedAbs},
    logic::traits::BitIterable,
};
use malachite_nz::{integer::Integer, natural::Natural};

#[derive(Debug)]
pub enum RingOppErr {
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

    fn div(a: Self, b: Self) -> Result<Self, RingOppErr>;
    fn div_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::div(a.clone(), b)
    }
    fn div_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::div(a, b.clone())
    }
    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::div(a.clone(), b.clone())
    }

    fn divisible(a: Self, b: Self) -> bool {
        match Self::div(a, b) {
            Ok(_q) => true,
            Err(RingOppErr::NotDivisible) => false,
            Err(RingOppErr::DivideByZero) => false,
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
            Err(RingOppErr::DivideByZero) => false,
            Err(RingOppErr::NotDivisible) => false,
            // Err(_) => panic!(),
        }
    }

    fn inv(self) -> Result<Self, RingOppErr> {
        Self::div(Self::one(), self)
    }

    fn inv_ref(a: &Self) -> Result<Self, RingOppErr> {
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

    fn factor_fav_assoc(self) -> (Self, Self);
    fn factor_fav_assoc_ref(&self) -> (Self, Self) {
        self.clone().factor_fav_assoc()
    }
    fn is_fav_assoc(&self) -> bool {
        let (_u, a) = self.clone().factor_fav_assoc();
        self == &a
    }
}

pub trait UniqueFactorizationDomain: IntegralDomain {
    //unique factorizations exist
}

pub trait UniquelyFactorable: UniqueFactorizationDomain + FavoriteAssociate + Hash {
    //UFD with an explicit factorizer
    fn make_factorizer() -> Box<dyn UniqueFactorizer<Self>>;
    fn factor(&self) -> Option<UniqueFactorization<Self>> {
        Self::make_factorizer().factor(self)
    }
    fn is_irreducible(&self) -> Option<bool> {
        Self::make_factorizer().is_irreducible(self)
    }
}

#[derive(Debug)]
pub struct UniqueFactorization<R: UniqueFactorizationDomain + FavoriteAssociate + Hash> {
    elem: R,
    unit: R,
    factors: HashMap<R, Natural>,
}

impl<R: UniqueFactorizationDomain + FavoriteAssociate + Hash> UniqueFactorization<R> {
    pub fn check_invariants(
        &self,
        factorizer: &mut impl UniqueFactorizer<R>,
    ) -> Result<(), &'static str> {
        if !self.unit.clone().is_unit() {
            return Err("unit must be a unit");
        }
        let mut prod = self.unit.clone();
        for (p, k) in &self.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if !p.is_fav_assoc() {
                return Err("prime factor must be their fav assoc");
            }
            match factorizer.is_irreducible(p) {
                Some(is_irr) => {
                    if !is_irr {
                        return Err("prime factor must not be irreducible");
                    }
                    let mut i = Natural::from(0u8);
                    while &i < k {
                        prod.mul_mut(p);
                        i += Natural::from(1u8);
                    }
                }
                None => {
                    return Err("prime factor must not be zero");
                }
            }
        }
        if self.elem != prod {
            return Err("product is incorrect");
        }
        Ok(())
    }

    pub fn new_unchecked(elem: R, unit: R, factors: HashMap<R, Natural>) -> Self {
        Self {
            elem,
            unit,
            factors,
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
}

pub trait UniqueFactorizer<R: UniqueFactorizationDomain + FavoriteAssociate + Hash> {
    //factor the non-zero elements
    fn factor(&mut self, a: &R) -> Option<UniqueFactorization<R>>;
    fn is_irreducible(&mut self, a: &R) -> Option<bool> {
        match self.factor(a) {
            Some(f) => Some(f.is_irreducible()),
            None => None,
        }
    }
}

pub trait PrincipalIdealDomain: ComRing + FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(x: Self, y: Self) -> Self;
    fn gcd_list(elems: Vec<&Self>) -> Self {
        let mut ans = Self::zero();
        for x in elems {
            ans = Self::gcd(ans, x.clone());
        }
        ans
    }
    fn xgcd(x: Self, y: Self) -> (Self, Self, Self);
    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        match elems.len() {
            0 => (Self::zero(), vec![]),
            1 => {
                let (unit, assoc) = elems[0].factor_fav_assoc_ref();
                (assoc, vec![unit])
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
    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr>;
    fn quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a.clone(), b)
    }
    fn quorem_rref(a: Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a, b.clone())
    }
    fn quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a.clone(), b.clone())
    }

    fn quo(a: Self, b: Self) -> Result<Self, RingOppErr> {
        match Self::quorem(a, b) {
            Ok((q, _r)) => Ok(q),
            Err(e) => Err(e),
        }
    }
    fn quo_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::quo(a.clone(), b)
    }
    fn quo_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::quo(a, b.clone())
    }
    fn quo_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::quo(a.clone(), b.clone())
    }

    fn rem(a: Self, b: Self) -> Result<Self, RingOppErr> {
        match Self::quorem(a, b) {
            Ok((_q, r)) => Ok(r),
            Err(e) => Err(e),
        }
    }
    fn rem_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::rem(a.clone(), b)
    }
    fn rem_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::rem(a, b.clone())
    }
    fn rem_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::rem(a.clone(), b.clone())
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> PrincipalIdealDomain for R {
    fn gcd(mut x: Self, mut y: Self) -> Self {
        //Euclidean algorithm
        while y != Self::zero() {
            let r = Self::rem_rref(x, &y).unwrap();
            (x, y) = (y, r)
        }
        let (_unit, assoc) = x.factor_fav_assoc();
        assoc
    }

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

// impl<R: EuclideanDomain> PrincipalIdealDomain for R {

// }

pub trait Field: IntegralDomain + UniqueFactorizationDomain {
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

    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_lref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_rref(a: Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_rref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_refs(a, b).unwrap(), Self::zero()))
        }
    }

    fn quo(a: Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div(a, b).unwrap())
        }
    }

    fn quo_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_lref(a, b).unwrap())
        }
    }

    fn quo_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_rref(a, b).unwrap())
        }
    }

    fn quo_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_refs(a, b).unwrap())
        }
    }

    fn rem(_a: Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_lref(_a: &Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_rref(_a: Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_refs(_a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
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

    #[test]
    fn factorization_invariants() {
        let mut naive_integer_factorizer = NaiveIntegerFactorizer();

        let f = UniqueFactorization {
            elem: Integer::from(-12),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            true
        );

        let f = UniqueFactorization {
            elem: Integer::from(1),
            unit: Integer::from(1),
            factors: HashMap::from([]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            true
        );

        let f = UniqueFactorization {
            elem: Integer::from(-12),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "can't have a power of zero"
        );

        let f = UniqueFactorization {
            elem: Integer::from(-13),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "product is incorrect"
        );

        let f = UniqueFactorization {
            elem: Integer::from(12),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "unit is wrong"
        );

        let f = UniqueFactorization {
            elem: Integer::from(12),
            unit: Integer::from(3),
            factors: HashMap::from([(Integer::from(2), Natural::from(2u8))]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "unit should be a unit"
        );

        let f = UniqueFactorization {
            elem: Integer::from(0),
            unit: Integer::from(1),
            factors: HashMap::from([
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "prime factors must not be zero"
        );

        let f = UniqueFactorization {
            elem: Integer::from(-12),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "prime factors must be prime"
        );

        let f = UniqueFactorization {
            elem: Integer::from(-12),
            unit: Integer::from(-1),
            factors: HashMap::from([
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ]),
        };
        assert_eq!(
            f.check_invariants(&mut naive_integer_factorizer).is_ok(),
            false,
            "prime factors must be fav assoc"
        );
    }
}
