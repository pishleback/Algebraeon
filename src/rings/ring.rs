use std::{
    borrow::Borrow,
    collections::HashMap,
    fmt::{Debug, Display},
    hash::Hash,
};

use malachite_base::num::{arithmetic::traits::UnsignedAbs, logic::traits::BitIterable};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

#[derive(Debug)]
pub enum RingDivisionError {
    DivideByZero,
    NotDivisible,
}

pub trait ComRing: Clone + Debug + PartialEq + Eq {
    //todo: remove sized here
    // type ElemT: Sized + Clone + Debug;

    // fn to_string(&self) -> String; //TODO: remove this, replace with std::fmt::Display triat

    // fn equal(&self, a: &Self::ElemT, b: &Self::ElemT) -> bool;

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

    fn add_mut(&mut self, offset: &Self);
    fn add(mut a: Self, b: Self) -> Self {
        Self::add_mut(&mut a, &b);
        a
    }
    fn add_ref(mut a: Self, b: &Self) -> Self {
        Self::add_mut(&mut a, b);
        a
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        Self::add_mut(&mut new_a, b);
        new_a
    }

    fn mul_mut(&mut self, mul: &Self);
    fn mul(mut a: Self, b: Self) -> Self {
        Self::mul_mut(&mut a, &b);
        a
    }
    fn mul_ref(mut a: Self, b: &Self) -> Self {
        Self::mul_mut(&mut a, b);
        a
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        Self::mul_mut(&mut new_a, b);
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

    fn are_associate(a: &Self, b: &Self) -> bool {
        if a == &Self::zero() && b == &Self::zero() {
            true
        } else {
            Self::div_refs(a, b).is_ok() && Self::div_refs(b, a).is_ok()
        }
    }

    fn sum(elems: Vec<&Self>) -> Self {
        let mut ans = Self::zero();
        for elem in elems {
            ans = Self::add_ref(ans, elem);
        }
        ans
    }

    fn product(elems: Vec<&Self>) -> Self {
        let mut ans = Self::one();
        for elem in elems {
            ans = Self::mul_ref(ans, elem);
        }
        ans
    }

    fn nat_pow(&self, n: impl Borrow<Natural>) -> Self {
        let n = n.borrow();
        if *n == 0 {
            Self::one()
        } else if *n == 1 {
            self.clone()
        } else {
            debug_assert!(*n >= 2);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![self.clone()];
            while pows.len() < bits.len() {
                pows.push(Self::mul_refs(&pows.last().unwrap(), &pows.last().unwrap()));
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = Self::one();
            for i in 0..count {
                if bits[i] {
                    Self::mul_mut(&mut ans, &pows[i]);
                }
            }
            ans
        }
    }

    fn int_pow(&self, n: &Integer) -> Option<Self> {
        // println!("{:?} {:?}", elem, n);
        if *n == 0 {
            Some(Self::one())
        } else if self == &Self::zero() {
            Some(Self::zero())
        } else if *n > 0 {
            Some(Self::nat_pow(self, &n.unsigned_abs()))
        } else {
            match Self::inv_ref(self) {
                Ok(self_inv) => Some(Self::nat_pow(&self_inv, &(-n).unsigned_abs())),
                Err(RingDivisionError::NotDivisible) => None,
                Err(RingDivisionError::DivideByZero) => panic!(),
            }
        }
    }

    fn from_int(x: &Integer) -> Self {
        if *x < 0 {
            Self::neg(Self::from_int(&-x))
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
                    Self::add_mut(&mut ans, &v);
                }
                Self::mul_mut(&mut v, &two);
            }
            ans
        }
    }

    fn from_rat(x: &Rational) -> Result<Self, RingDivisionError> {
        Self::div(
            Self::from_int(&Rational::numerator(x)),
            Self::from_int(&Rational::denominator(x)),
        )
    }

    fn is_unit(elem: Self) -> bool {
        match Self::div(Self::one(), elem) {
            Ok(_inv) => true,
            Err(RingDivisionError::DivideByZero) => false,
            Err(RingDivisionError::NotDivisible) => false,
            // Err(_) => panic!(),
        }
    }

    fn inv(self) -> Result<Self, RingDivisionError> {
        Self::div(Self::one(), self)
    }

    fn inv_ref(&self) -> Result<Self, RingDivisionError> {
        Self::div_rref(Self::one(), self)
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
            .map(|n| Self::from_int(&n)),
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
    fn fav_assoc(self) -> Self {
        self.factor_fav_assoc().1
    }
    fn factor_fav_assoc_ref(&self) -> (Self, Self) {
        self.clone().factor_fav_assoc()
    }
    fn is_fav_assoc(&self) -> bool {
        let (_u, a) = Self::factor_fav_assoc_ref(self);
        self == &a
    }
}

pub trait GreatestCommonDivisorDomain: FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(x: Self, y: Self) -> Self;
    fn gcd_list<BorrowSelfT: Borrow<Self>>(elems: Vec<BorrowSelfT>) -> Self {
        let mut ans = Self::zero();
        for x in elems {
            ans = Self::gcd(ans, x.borrow().clone());
        }
        ans
    }
    fn lcm(x: Self, y: Self) -> Self {
        if x == Self::zero() && y == Self::zero() {
            Self::zero()
        } else {
            let g = Self::gcd(x.clone(), y.clone());
            Self::div(Self::mul(x, y), g).unwrap()
        }
    }

    fn lcm_list<BorrowSelfT: Borrow<Self>>(elems: Vec<BorrowSelfT>) -> Self {
        let mut ans = Self::one();
        for x in elems {
            ans = Self::lcm(ans, x.borrow().clone());
        }
        ans
    }
}

#[derive(Debug)]
pub struct Factored<Ring: UniqueFactorizationDomain> {
    unit: Ring,
    factors: Vec<(Ring, Natural)>,
}

impl<Ring: UniqueFactorizationDomain + Display> Display for Factored<Ring> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.unit);
        for (factor, k) in &self.factors {
            write!(f, " * (");
            write!(f, "{}", factor);
            write!(f, ")");
            if k != &Natural::from(1u8) {
                write!(f, "^");
                write!(f, "{}", k);
            }
        }
        Ok(())
    }
}

impl<Ring: UniqueFactorizationDomain> Factored<Ring> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if !Ring::is_unit(self.unit.clone()) {
            return Err("unit must be a unit");
        }
        for (p, k) in &self.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if p == &Ring::zero() {
                return Err("prime factor must not be zero");
            }
            if !Ring::is_fav_assoc(p) {
                return Err("prime factor must be their fav assoc");
            }
            if !Ring::is_irreducible(p).unwrap() {
                return Err("prime factor must not be reducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                i += Natural::from(1u8);
            }
        }
        Ok(())
    }

    pub fn expand(&self) -> Ring {
        let mut ans = self.unit.clone();
        for (p, k) in &self.factors {
            Ring::mul_mut(&mut ans, &Ring::nat_pow(p, k));
        }
        ans
    }

    pub fn new_unchecked(unit: Ring, factors: Vec<(Ring, Natural)>) -> Self {
        Self { unit, factors }
    }

    pub fn equal(a: &Self, b: &Self) -> bool {
        a.unit() == b.unit()
            && a.factors()
                .iter()
                .map(|(x, k)| (x, k))
                .collect::<HashMap<_, _>>()
                == b.factors()
                    .iter()
                    .map(|(x, k)| (x, k))
                    .collect::<HashMap<_, _>>()
    }

    pub fn unit(&self) -> &Ring {
        &self.unit
    }

    pub fn factors(&self) -> &Vec<(Ring, Natural)> {
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

    fn mul_by_unchecked(&mut self, p: Ring, k: Natural) {
        for (q, t) in &mut self.factors {
            match (Ring::div_refs(&p, q), Ring::div_refs(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if Ring::is_unit(u) && Ring::is_unit(v.clone()) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        Ring::mul_mut(&mut self.unit, &Ring::nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        self.factors.push((p, k));
    }

    pub fn mul(mut a: Self, b: Self) -> Self {
        Ring::mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            a.mul_by_unchecked(p, k)
        }
        a
    }

    pub fn factored_one() -> Self {
        Factored {
            unit: Ring::one(),
            factors: vec![],
        }
    }

    pub fn factored_irreducible_unchecked(elem: Ring) -> Self {
        let (unit, assoc) = Ring::factor_fav_assoc(elem);
        Factored {
            unit,
            factors: vec![(assoc, Natural::from(1u8))],
        }
    }

    pub fn factored_unit_unchecked(unit: Ring) -> Self {
        Factored {
            unit,
            factors: vec![],
        }
    }
}

pub fn full_factor_using_partial_factor<R: UniqueFactorizationDomain>(
    elem: R,
    partial_factor: &impl Fn(R) -> Option<(R, R)>,
) -> Factored<R> {
    pub fn full_factor_using_partial_factor_no_unit<R: UniqueFactorizationDomain>(
        elem: R,
        partial_factor: &impl Fn(R) -> Option<(R, R)>,
    ) -> Factored<R> {
        debug_assert_ne!(elem, R::zero());
        debug_assert!(!R::is_unit(elem.clone()));
        match partial_factor(elem.clone()) {
            Some((g, h)) => Factored::mul(
                full_factor_using_partial_factor(g, partial_factor),
                full_factor_using_partial_factor(h, partial_factor),
            ),
            None => {
                //f is irreducible
                Factored::factored_irreducible_unchecked(elem)
            }
        }
    }

    if R::is_unit(elem.clone()) {
        Factored::factored_unit_unchecked(elem)
    } else {
        full_factor_using_partial_factor_no_unit(elem, partial_factor)
    }
}

impl<Ring: UniqueFactorizationDomain> PartialEq for Factored<Ring> {
    fn eq(&self, other: &Self) -> bool {
        self.unit == other.unit && self.factors == other.factors
    }
}

impl<Ring: UniqueFactorizationDomain> Eq for Factored<Ring> {}

pub trait UniqueFactorizationDomain: FavoriteAssociate + Hash {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self) -> Option<Factored<Self>>;

    fn is_irreducible(&self) -> Option<bool> {
        match self.factor() {
            None => None,
            Some(factored) => Some(factored.is_irreducible()),
        }
    }

    fn divisors<'a>(factored: &'a Factored<Self>) -> Box<dyn Iterator<Item = Self> + 'a> {
        if factored.factors.len() == 0 {
            Box::new(vec![Self::one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in &factored.factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = Self::one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = Self::mul_ref(p_pow, &p);
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(|prime_power_factors| {
                    Self::product(prime_power_factors.iter().collect()).clone()
                }),
            )
        }
    }

    fn count_divisors(factored: &Factored<Self>) -> Option<Natural> {
        let mut count = Natural::from(1u8);
        for (_p, k) in &factored.factors {
            count *= k + Natural::from(1u8);
        }
        Some(count)
    }
}

pub trait PrincipalIdealDomain: GreatestCommonDivisorDomain {
    //any gcds should be the standard associate representative
    fn xgcd(a: Self, b: Self) -> (Self, Self, Self); //(g, x, y) s.t. g = ax + by
    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        match elems.len() {
            0 => (Self::zero(), vec![]),
            1 => {
                let (unit, assoc) = Self::factor_fav_assoc_ref(elems[0]);
                (assoc, vec![Self::inv(unit).unwrap()])
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
    }
}

pub trait EuclideanDomain: IntegralDomain {
    //should return None for 0, and Some(norm) for everything else
    fn norm(elem: &Self) -> Option<Natural>;
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
        let (_unit, assoc) = Self::factor_fav_assoc(x);
        assoc
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> PrincipalIdealDomain for R {
    fn xgcd(mut x: Self, mut y: Self) -> (Self, Self, Self) //return (g, a, b)
    {
        let mut pa = Self::one();
        let mut a = Self::zero();
        let mut pb = Self::zero();
        let mut b = Self::one();

        while y != Self::zero() {
            let (q, r) = Self::quorem_rref(x, &y).unwrap();
            let new_a = Self::add(pa, Self::neg(Self::mul_refs(&q, &a)));
            (a, pa) = (new_a, a);
            let new_b = Self::add(pb, Self::neg(Self::mul_ref(q, &b)));
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = Self::factor_fav_assoc(x);
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(Self::is_unit(unit.clone()));
        (
            ass_x,
            Self::div_rref(pa, &unit).unwrap(),
            Self::div(pb, unit).unwrap(),
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
    fn norm(elem: &Self) -> Option<Natural> {
        if elem == &Self::zero() {
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

// impl<F: Field + Hash> UniqueFactorizationDomain for F
// where
//     Self::ElemT: Hash,
// {
//     fn factor(&self, elem: &Self::ElemT) -> Option<Factored<Self::ElemT>> {
//         if self.equal(elem, &self.zero()) {
//             None
//         } else {
//             Some(Factored::new_unchecked(elem.clone(), vec![]))
//         }
//     }
// }

pub trait FieldOfFractions: Field {
    type R: IntegralDomain;

    fn from_base_ring(elem: Self::R) -> Self;
    fn numerator(elem: &Self) -> Self::R;
    fn denominator(elem: &Self) -> Self::R;
    fn as_base_ring(elem: Self) -> Option<Self::R> {
        if Self::denominator(&elem) == Self::R::one() {
            Some(Self::numerator(&elem))
        } else {
            None
        }
    }
}

#[derive(Debug, Clone)]
pub struct EuclideanQuotient<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain>
{
    a: ED,
    n: Option<ED>,
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> PartialEq
    for EuclideanQuotient<IS_FIELD, ED>
{
    fn eq(&self, other: &Self) -> bool {
        let n_opt = match (&self.n, &other.n) {
            (None, None) => None,
            (None, Some(other_n)) => Some(other_n),
            (Some(self_n), None) => Some(self_n),
            (Some(self_n), Some(other_n)) => {
                assert_eq!(self_n, other_n);
                Some(self_n)
            }
        };
        match n_opt {
            Some(n) => {
                ED::rem_refs(&ED::add_ref(other.a.neg_ref(), &self.a), n).unwrap() == ED::zero()
            }
            None => self.a == other.a,
        }
    }
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> Eq
    for EuclideanQuotient<IS_FIELD, ED>
{
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain>
    EuclideanQuotient<IS_FIELD, ED>
{
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match &self.n {
            Some(n) => {
                if n == &ED::zero() {
                    return Err("Can't quotient by zero");
                }
                if IS_FIELD {
                    if !n.is_irreducible().unwrap() {
                        return Err("Marked as field but modulus is not irreducible");
                    }
                }
            }
            None => {}
        }
        Ok(())
    }

    pub fn new(a: ED, n: ED) -> Self {
        debug_assert!(n != ED::zero());
        if IS_FIELD {
            debug_assert!(n.is_irreducible().unwrap());
        }
        Self { a, n: Some(n) }
    }

    pub fn lift(self) -> ED {
        self.a
    }

    pub fn enforce_modulus(&mut self, n: ED) {
        match &self.n {
            Some(self_n) => {
                assert_eq!(self_n, &n);
            }
            None => {
                self.n = Some(n);
            }
        }
    }
}

impl<const IS_FIELD: bool, ED: EuclideanDomain + UniqueFactorizationDomain> ComRing
    for EuclideanQuotient<IS_FIELD, ED>
{
    fn zero() -> Self {
        Self {
            a: ED::zero(),
            n: None,
        }
    }

    fn one() -> Self {
        Self {
            a: ED::one(),
            n: None,
        }
    }

    fn neg_mut(&mut self) {
        self.a.neg_mut()
    }

    fn add_mut(&mut self, other: &Self) {
        self.a.add_mut(&other.a);
        match (&self.n, &other.n) {
            (None, None) => {}
            (None, Some(other_n)) => {
                self.a = ED::rem_refs(&self.a, other_n).unwrap();
            }
            (Some(self_n), None) => {
                self.a = ED::rem_refs(&self.a, self_n).unwrap();
            }
            (Some(self_n), Some(other_n)) => {
                assert_eq!(self_n, other_n);
                self.a = ED::rem_refs(&self.a, self_n).unwrap();
            }
        }
    }

    fn mul_mut(&mut self, other: &Self) {
        self.a.mul_mut(&other.a);
        match (&self.n, &other.n) {
            (None, None) => {}
            (None, Some(other_n)) => {
                self.a = ED::rem_refs(&self.a, other_n).unwrap();
            }
            (Some(self_n), None) => {
                self.a = ED::rem_refs(&self.a, self_n).unwrap();
            }
            (Some(self_n), Some(other_n)) => {
                assert_eq!(self_n, other_n);
                self.a = ED::rem_refs(&self.a, self_n).unwrap();
            }
        }
    }

    fn div(top: Self, bot: Self) -> Result<Self, RingDivisionError> {
        let n = {
            match (&top.n, &bot.n) {
                (None, None) => {
                    panic!("At least one modulus must be present for division in quotient of Euclidean domain");
                }
                (None, Some(bot_n)) => bot_n,
                (Some(top_n), None) => top_n,
                (Some(top_n), Some(bot_n)) => {
                    assert_eq!(top_n, bot_n);
                    top_n
                }
            }
        };

        if bot == Self::zero() {
            Err(RingDivisionError::DivideByZero)
        } else {
            let (g, _x, y) = ED::xgcd(n.clone(), bot.a);
            //g = xn + yb  so   g = yb mod n
            //if z = a/g works then
            //zyb = a mod n
            match ED::div(top.a, g) {
                Ok(z) => Ok(Self {
                    a: ED::mul(z, y),
                    n: Some(n.clone()),
                }),
                Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
                Err(RingDivisionError::DivideByZero) => panic!(),
            }
        }
    }
}

impl<ED: EuclideanDomain + UniqueFactorizationDomain> IntegralDomain
    for EuclideanQuotient<true, ED>
{
}

impl<ED: EuclideanDomain + UniqueFactorizationDomain> Field for EuclideanQuotient<true, ED> {}

//IdealQuotient
//PrimeQuotient
//MaximalQuotient

pub trait Real: ComRing {
    fn as_f64(x: Self) -> f64;
    fn as_f32(x: Self) -> f32 {
        Self::as_f64(x) as f32
    }
}

impl<F: FieldOfFractions> Real for F
where
    F::R: Real,
{
    fn as_f64(x: Self) -> f64 {
        F::R::as_f64(Self::numerator(&x)) / F::R::as_f64(Self::denominator(&x))
    }
}

pub trait DenseReal: Real {
    fn from_f64_approx(x: f64) -> Self;
    fn from_f32_approx(x: f32) -> Self {
        Self::from_f64_approx(x as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::super::nzq::*;
    use super::*;

    // #[test]
    // fn factorization_invariants() {
    //     let f = Factored::new_unchecked(
    //         Integer::from(-1),
    //         vec![
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ],
    //     );
    //     f.check_invariants(&ZZ).unwrap();

    //     let f = Factored::new_unchecked(Integer::from(1), vec![]);
    //     f.check_invariants(&ZZ).unwrap();

    //     let f = Factored::new_unchecked(
    //         Integer::from(-1),
    //         vec![
    //             (Integer::from(2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //             (Integer::from(5), Natural::from(0u8)),
    //         ],
    //     );
    //     assert_eq!(
    //         f.check_invariants(&ZZ).is_ok(),
    //         false,
    //         "can't have a power of zero"
    //     );

    //     let f = Factored::new_unchecked(
    //         Integer::from(3),
    //         vec![(Integer::from(2), Natural::from(2u8))],
    //     );
    //     assert_eq!(
    //         f.check_invariants(&ZZ).is_ok(),
    //         false,
    //         "unit should be a unit"
    //     );

    //     let f = Factored::new_unchecked(
    //         Integer::from(1),
    //         vec![
    //             (Integer::from(0), Natural::from(1u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ],
    //     );
    //     assert_eq!(
    //         f.check_invariants(&ZZ).is_ok(),
    //         false,
    //         "prime factors must not be zero"
    //     );

    //     let f = Factored::new_unchecked(
    //         Integer::from(-1),
    //         vec![
    //             (Integer::from(4), Natural::from(1u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ],
    //     );
    //     assert_eq!(
    //         f.check_invariants(&ZZ).is_ok(),
    //         false,
    //         "prime factors must be prime"
    //     );

    //     let f = Factored::new_unchecked(
    //         Integer::from(-1),
    //         vec![
    //             (Integer::from(-2), Natural::from(2u8)),
    //             (Integer::from(3), Natural::from(1u8)),
    //         ],
    //     );
    //     assert_eq!(
    //         f.check_invariants(&ZZ).is_ok(),
    //         false,
    //         "prime factors must be fav assoc"
    //     );
    // }

    #[test]
    fn test_xgcd_list() {
        use malachite_q::Rational;
        let a = Rational::from(7);
        let (g, taps) = Rational::xgcd_list(vec![&a]);
        assert_eq!(g, Rational::one());
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
