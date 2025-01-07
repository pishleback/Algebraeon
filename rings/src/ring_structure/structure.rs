use std::{borrow::Borrow, fmt::Debug, rc::Rc};

use malachite_base::num::arithmetic::traits::UnsignedAbs;
use malachite_base::num::logic::traits::BitIterable;
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::number::natural::nat_to_usize;
use crate::polynomial::polynomial::Polynomial;
use crate::polynomial::polynomial::PolynomialStructure;

use algebraeon_sets::structure::*;

use super::factorization::*;

#[derive(Debug, PartialEq, Eq)]
pub enum RingDivisionError {
    DivideByZero,
    NotDivisible,
}

pub trait RingStructure: EqStructure {
    fn is_zero(&self, a: &Self::Set) -> bool {
        self.equal(a, &self.zero())
    }

    fn zero(&self) -> Self::Set;
    fn one(&self) -> Self::Set;

    fn neg(&self, a: &Self::Set) -> Self::Set;
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn add_mut(&self, a: &mut Self::Set, b: &Self::Set) {
        *a = self.add(a, b);
    }
    fn sum(&self, vals: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut sum = self.zero();
        for val in vals {
            self.add_mut(&mut sum, val.borrow());
        }
        sum
    }
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn mul_mut(&self, a: &mut Self::Set, b: &Self::Set) {
        *a = self.mul(a, b);
    }
    fn product(&self, vals: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut prod = self.one();
        for val in vals {
            self.mul_mut(&mut prod, val.borrow());
        }
        prod
    }

    fn nat_pow(&self, a: &Self::Set, n: &Natural) -> Self::Set {
        if *n == 0 {
            self.one()
        } else if *n == 1 {
            a.clone()
        } else {
            debug_assert!(*n >= 2);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(self.mul(&pows.last().unwrap(), &pows.last().unwrap()));
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
    fn from_int(&self, x: &Integer) -> Self::Set {
        if *x < 0 {
            self.neg(&self.from_int(&-x))
        } else if *x == 0 {
            self.zero()
        } else if *x == 1 {
            self.one()
        } else {
            let two = self.add(&self.one(), &self.one());
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
}

pub trait MetaRing: MetaType
where
    Self::Structure: RingStructure,
{
    fn zero() -> Self {
        Self::structure().zero()
    }
    fn one() -> Self {
        Self::structure().one()
    }

    fn neg(&self) -> Self {
        Self::structure().neg(self)
    }
    fn add(a: &Self, b: &Self) -> Self {
        Self::structure().add(a, b)
    }
    fn sum(vals: Vec<impl Borrow<Self>>) -> Self {
        Self::structure().sum(vals)
    }
    fn mul(a: &Self, b: &Self) -> Self {
        Self::structure().mul(a, b)
    }
    fn product(vals: Vec<impl Borrow<Self>>) -> Self {
        Self::structure().product(vals)
    }

    fn nat_pow(&self, n: &Natural) -> Self {
        Self::structure().nat_pow(self, n)
    }
    fn from_int(x: &Integer) -> Self {
        Self::structure().from_int(x)
    }
}
impl<R: MetaType> MetaRing for R where Self::Structure: RingStructure {}

pub trait MetaRingEq: MetaType
where
    Self::Structure: RingStructure + EqStructure,
{
    fn is_zero(&self) -> bool {
        Self::structure().is_zero(self)
    }
}
impl<R: MetaType> MetaRingEq for R where Self::Structure: RingStructure + EqStructure {}

// pub trait DisplayableRingStructure: RingStructure {
//     fn fmt_elem(&self, elem: &Self::Set, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result;
// }

pub trait IntegralDomainStructure: RingStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError>;
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }

    fn from_rat(&self, x: &Rational) -> Option<Self::Set> {
        match self.div(
            &self.from_int(&x.numerator()),
            &self.from_int(&x.denominator()),
        ) {
            Ok(d) => Some(d),
            Err(RingDivisionError::NotDivisible) => None,
            Err(RingDivisionError::DivideByZero) => panic!(),
        }
    }

    fn int_pow(&self, a: &Self::Set, n: &Integer) -> Option<Self::Set> {
        // println!("{:?} {:?}", elem, n);
        if *n == 0 {
            Some(self.one())
        } else if self.is_zero(a) {
            Some(self.zero())
        } else if *n > 0 {
            Some(self.nat_pow(a, &n.unsigned_abs()))
        } else {
            match self.inv(a) {
                Ok(self_inv) => Some(self.nat_pow(&self_inv, &(-n).unsigned_abs())),
                Err(RingDivisionError::NotDivisible) => None,
                Err(RingDivisionError::DivideByZero) => panic!(),
            }
        }
    }

    fn divisible(&self, a: &Self::Set, b: &Self::Set) -> bool {
        match self.div(a, b) {
            Ok(_q) => true,
            Err(RingDivisionError::NotDivisible) => false,
            Err(RingDivisionError::DivideByZero) => false,
        }
    }
    fn are_associate(&self, a: &Self::Set, b: &Self::Set) -> bool {
        if self.equal(a, &self.zero()) && self.equal(b, &self.zero()) {
            true
        } else {
            self.div(a, b).is_ok() && self.div(b, a).is_ok()
        }
    }
    fn is_unit(&self, a: &Self::Set) -> bool {
        match self.div(&self.one(), &a) {
            Ok(_inv) => true,
            Err(RingDivisionError::DivideByZero) => false,
            Err(RingDivisionError::NotDivisible) => false,
        }
    }
}
pub trait MetaIntegralDomain: MetaRing
where
    Self::Structure: IntegralDomainStructure,
{
    fn div(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::structure().div(a, b)
    }
    fn inv(&self) -> Result<Self, RingDivisionError> {
        Self::structure().inv(self)
    }

    fn from_rat(x: &Rational) -> Option<Self> {
        Self::structure().from_rat(x)
    }

    fn int_pow(&self, n: &Integer) -> Option<Self> {
        Self::structure().int_pow(self, n)
    }

    fn divisible(a: &Self, b: &Self) -> bool {
        Self::structure().divisible(a, b)
    }
    fn are_associate(a: &Self, b: &Self) -> bool {
        Self::structure().are_associate(a, b)
    }
    fn is_unit(&self) -> bool {
        Self::structure().is_unit(self)
    }
}
impl<R: MetaRing> MetaIntegralDomain for R where Self::Structure: IntegralDomainStructure<Set = R> {}

pub trait OrderedRingStructure: IntegralDomainStructure {
    // <= satisfying translation invariance and multiplication by positive scalar
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering;
    fn abs(&self, a: &Self::Set) -> Self::Set {
        match self.ring_cmp(a, &self.zero()) {
            std::cmp::Ordering::Less => self.neg(a),
            std::cmp::Ordering::Equal => self.zero(),
            std::cmp::Ordering::Greater => a.clone(),
        }
    }
}
pub trait MetaOrderedRing: MetaType
where
    Self::Structure: OrderedRingStructure,
{
    fn ring_cmp(a: &Self, b: &Self) -> std::cmp::Ordering {
        Self::structure().ring_cmp(a, b)
    }

    fn abs(a: &Self) -> Self {
        Self::structure().abs(a)
    }
}
impl<R: MetaType> MetaOrderedRing for R where Self::Structure: OrderedRingStructure<Set = R> {}

pub trait FiniteUnitsStructure: RingStructure {
    fn all_units(&self) -> Vec<Self::Set>;
}

pub trait MetaFiniteUnits: MetaRing
where
    Self::Structure: FiniteUnitsStructure,
{
    fn all_units(&self) -> Vec<Self> {
        Self::structure().all_units()
    }
}
impl<R: MetaRing> MetaFiniteUnits for R where Self::Structure: FiniteUnitsStructure<Set = R> {}

pub trait FavoriteAssociateStructure: IntegralDomainStructure {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)
    //every unit u is required to return (u, 1)

    //it seems to happen that the product of favorite associates is another favorite associate. Should this be a requirement?

    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set);
    fn fav_assoc(&self, a: &Self::Set) -> Self::Set {
        self.factor_fav_assoc(a).1
    }
    fn is_fav_assoc(&self, a: &Self::Set) -> bool {
        let (_u, b) = self.factor_fav_assoc(a);
        self.equal(a, &b)
    }
}
pub trait MetaFavoriteAssociate: MetaIntegralDomain
where
    Self::Structure: FavoriteAssociateStructure,
{
    fn factor_fav_assoc(&self) -> (Self, Self) {
        Self::structure().factor_fav_assoc(self)
    }

    fn fav_assoc(&self) -> Self {
        Self::structure().fav_assoc(self)
    }
}
impl<R: MetaRing> MetaFavoriteAssociate for R where
    Self::Structure: FavoriteAssociateStructure<Set = R>
{
}

pub trait GreatestCommonDivisorStructure: FavoriteAssociateStructure {
    //any gcds should be the standard associate representative
    //euclidean_gcd can be used to implement this
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
    fn gcd_list(&self, elems: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut gcd = self.zero();
        for x in elems {
            gcd = self.gcd(&gcd, x.borrow());
        }
        gcd
    }
    fn lcm(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        if self.is_zero(x) && self.is_zero(y) {
            self.zero()
        } else {
            self.div(&self.mul(x, y), &self.gcd(x, y)).unwrap()
        }
    }
    fn lcm_list(&self, elems: Vec<impl Borrow<Self::Set>>) -> Self::Set {
        let mut lcm = self.one();
        for x in elems {
            lcm = self.lcm(&lcm, x.borrow());
        }
        lcm
    }
}
pub trait MetaGreatestCommonDivisor: MetaFavoriteAssociate
where
    Self::Structure: GreatestCommonDivisorStructure,
{
    fn gcd(x: &Self, y: &Self) -> Self {
        Self::structure().gcd(x, y)
    }

    fn gcd_list(elems: Vec<impl Borrow<Self>>) -> Self {
        Self::structure().gcd_list(elems)
    }

    fn lcm(x: &Self, y: &Self) -> Self {
        Self::structure().lcm(x, y)
    }

    fn lcm_list(elems: Vec<impl Borrow<Self>>) -> Self {
        Self::structure().lcm_list(elems)
    }
}
impl<R: MetaRing> MetaGreatestCommonDivisor for R where
    Self::Structure: GreatestCommonDivisorStructure<Set = R>
{
}

pub trait UniqueFactorizationStructure: FavoriteAssociateStructure {
    //a UFD with an explicit algorithm to compute unique factorizations
    fn factor(&self, a: &Self::Set) -> Option<super::factorization::Factored<Self>>;

    fn is_irreducible(&self, a: &Self::Set) -> bool {
        match self.factor(a) {
            None => false, //zero is not irreducible
            Some(factored) => factored.is_irreducible(),
        }
    }
}
pub trait MetaUniqueFactorization: MetaFavoriteAssociate
where
    Self::Structure: UniqueFactorizationStructure,
{
    fn factor(&self) -> Option<Factored<Self::Structure>> {
        Self::structure().factor(self)
    }

    fn is_irreducible(&self) -> bool {
        Self::structure().is_irreducible(self)
    }
}
impl<R: MetaRing> MetaUniqueFactorization for R where
    Self::Structure: UniqueFactorizationStructure<Set = R>
{
}

pub trait BezoutDomainStructure: GreatestCommonDivisorStructure {
    //any gcds should be the standard associate representative
    fn xgcd(&self, a: &Self::Set, b: &Self::Set) -> (Self::Set, Self::Set, Self::Set); //(g, x, y) s.t. g = ax + by
    fn xgcd_list(&self, elems: Vec<&Self::Set>) -> (Self::Set, Vec<Self::Set>) {
        // println!("{:?}", elems);
        match elems.len() {
            0 => (self.zero(), vec![]),
            1 => {
                let (unit, assoc) = self.factor_fav_assoc(elems[0]);
                (assoc, vec![self.inv(&unit).unwrap()])
            }
            2 => {
                let (g, x, y) = self.xgcd(elems[0], elems[1]);
                (g, vec![x, y])
            }
            n => {
                let k = n / 2;
                let (g1, coeffs1) = self.xgcd_list((0..k).map(|i| elems[i]).collect());
                let (g2, coeffs2) = self.xgcd_list((k..n).map(|i| elems[i]).collect());
                let (g, x, y) = self.xgcd(&g1, &g2);
                let mut coeffs = vec![];
                for c in coeffs1 {
                    coeffs.push(self.mul(&x, &c));
                }
                for c in coeffs2 {
                    coeffs.push(self.mul(&y, &c));
                }
                (g, coeffs)
            }
        }
    }
}

pub trait MetaBezoutDomain: MetaGreatestCommonDivisor
where
    Self::Structure: BezoutDomainStructure,
{
    fn xgcd(x: &Self, y: &Self) -> (Self, Self, Self) {
        Self::structure().xgcd(x, y)
    }

    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        Self::structure().xgcd_list(elems)
    }
}
impl<R: MetaRing> MetaBezoutDomain for R where Self::Structure: BezoutDomainStructure<Set = R> {}

pub trait EuclideanDivisionStructure: IntegralDomainStructure {
    //should return None for 0, and Some(norm) for everything else
    fn norm(&self, elem: &Self::Set) -> Option<Natural>;

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)>;

    fn quo(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        match self.quorem(a, b) {
            Some((q, _r)) => Some(q),
            None => None,
        }
    }

    fn rem(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        if !self.is_zero(&b) {
            let (_q, r) = self.quorem(a, b).unwrap();
            r
        } else {
            a.clone()
        }
    }

    fn euclidean_gcd(&self, mut x: Self::Set, mut y: Self::Set) -> Self::Set
    where
        Self: FavoriteAssociateStructure,
    {
        //Euclidean algorithm
        while !self.is_zero(&y) {
            let r = self.rem(&x, &y);
            (x, y) = (y, r)
        }
        let (_unit, assoc) = self.factor_fav_assoc(&x);
        assoc
    }

    fn euclidean_xgcd(
        &self,
        mut x: Self::Set,
        mut y: Self::Set,
    ) -> (Self::Set, Self::Set, Self::Set)
    where
        Self: FavoriteAssociateStructure,
    {
        let orig_x = x.clone();
        let orig_y = y.clone();

        let mut pa = self.one();
        let mut a = self.zero();
        let mut pb = self.zero();
        let mut b = self.one();

        while !self.is_zero(&y) {
            let (q, r) = self.quorem(&x, &y).unwrap();
            let new_a = self.add(&pa, &self.neg(&self.mul(&q, &a)));
            (a, pa) = (new_a, a);
            let new_b = self.add(&pb, &self.neg(&self.mul(&q, &b)));
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = self.factor_fav_assoc(&x);
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(self.is_unit(&unit));
        let (g, a, b) = (
            ass_x,
            self.div(&pa, &unit).unwrap(),
            self.div(&pb, &unit).unwrap(),
        );
        // println!("{:?} = {:?} * {:?} + {:?} * {:?}", g, a, orig_x, b, orig_y);
        debug_assert!(self.equal(
            &self.add(&self.mul(&a, &orig_x), &self.mul(&b, &orig_y)),
            &g
        ));
        debug_assert!(self.equal(&g, &self.euclidean_gcd(orig_x, orig_y)));
        (g, a, b)
    }
}

pub trait MetaEuclideanDivision: MetaIntegralDomain
where
    Self::Structure: EuclideanDivisionStructure,
{
    fn norm(&self) -> Option<malachite_nz::natural::Natural> {
        Self::structure().norm(self)
    }

    fn quorem(a: &Self, b: &Self) -> Option<(Self, Self)> {
        Self::structure().quorem(a, b)
    }

    fn quo(a: &Self, b: &Self) -> Option<Self> {
        Self::structure().quo(a, b)
    }

    fn rem(a: &Self, b: &Self) -> Self {
        Self::structure().rem(a, b)
    }

    fn euclidean_gcd(x: Self, y: Self) -> Self
    where
        Self::Structure: FavoriteAssociateStructure,
    {
        Self::structure().euclidean_gcd(x, y)
    }

    fn euclidean_xgcd(x: Self, y: Self) -> (Self, Self, Self)
    where
        Self::Structure: GreatestCommonDivisorStructure,
    {
        Self::structure().euclidean_xgcd(x, y)
    }
}
impl<R: MetaRing> MetaEuclideanDivision for R where
    Self::Structure: EuclideanDivisionStructure<Set = R>
{
}

pub trait InfiniteStructure: Structure {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::Set>>;
}
pub trait Infinite: MetaType {
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>>;
}
impl<T: MetaType> Infinite for T
where
    CannonicalStructure<T>: InfiniteStructure,
{
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>> {
        todo!()
    }
}

pub trait CharZeroStructure: RingStructure {}
pub trait MetaCharZero: MetaRing
where
    Self::Structure: CharZeroStructure,
{
}
impl<R: MetaType> MetaCharZero for R where Self::Structure: CharZeroStructure<Set = R> {}

impl<RS: CharZeroStructure + 'static> InfiniteStructure for RS {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = <Self as Structure>::Set>> {
        struct IntegerIterator<RS: CharZeroStructure> {
            ring: RS,
            next: Integer,
        }

        impl<RS: CharZeroStructure> Iterator for IntegerIterator<RS> {
            type Item = RS::Set;

            fn next(&mut self) -> Option<Self::Item> {
                let next = self.next.clone();
                if 0 < next {
                    self.next = -self.next.clone();
                } else {
                    self.next = Integer::from(1) - self.next.clone();
                }
                Some(self.ring.from_int(&next))
            }
        }

        Box::new(IntegerIterator {
            ring: self.clone(),
            next: Integer::from(0),
        })
    }
}

pub trait FieldStructure: IntegralDomainStructure {}

impl<FS: FieldStructure> FavoriteAssociateStructure for FS {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        if self.is_zero(a) {
            (self.one(), self.zero())
        } else {
            (a.clone(), self.one())
        }
    }
}

impl<FS: FieldStructure> EuclideanDivisionStructure for FS {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if self.is_zero(elem) {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        if self.is_zero(b) {
            None
        } else {
            Some((self.div(a, b).unwrap(), self.zero()))
        }
    }
}

impl<FS: FieldStructure> GreatestCommonDivisorStructure for FS {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.euclidean_gcd(x.clone(), y.clone())
    }
}

impl<FS: FieldStructure> BezoutDomainStructure for FS {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        self.euclidean_xgcd(x.clone(), y.clone())
    }
}

impl<FS: FieldStructure> UniqueFactorizationStructure for FS {
    fn factor(&self, a: &Self::Set) -> Option<super::factorization::Factored<Self>> {
        if self.is_zero(a) {
            None
        } else {
            Some(Factored::new_unchecked(
                Rc::new(self.clone()),
                a.clone(),
                vec![],
            ))
        }
    }
}

pub trait FiniteFieldStructure: FieldStructure + FiniteUnitsStructure {
    //return (p, k) where p is a prime and |F| = p^k
    fn characteristic_and_power(&self) -> (Natural, Natural);

    fn all_elements(&self) -> Vec<Self::Set> {
        let mut elems = vec![self.zero()];
        elems.append(&mut self.all_units());
        elems
    }
}

pub trait FieldOfFractionsStructure: FieldStructure {
    type RS: IntegralDomainStructure;

    fn base_ring_structure(&self) -> Rc<Self::RS>;
    fn from_base_ring(&self, elem: <Self::RS as Structure>::Set) -> Self::Set;
    fn numerator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set;
    fn denominator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set;
    fn as_base_ring(&self, elem: Self::Set) -> Option<<Self::RS as Structure>::Set> {
        let base_ring = self.base_ring_structure();
        if base_ring.equal(&self.denominator(&elem), &base_ring.one()) {
            Some(self.numerator(&elem))
        } else {
            None
        }
    }
}
pub trait MetaFieldOfFractions: MetaRing
where
    Self::Structure: FieldOfFractionsStructure,
{
    fn base_ring_structure() -> Rc<<Self::Structure as FieldOfFractionsStructure>::RS> {
        Self::structure().base_ring_structure()
    }
    fn from_base_ring(
        a: <<Self::Structure as FieldOfFractionsStructure>::RS as Structure>::Set,
    ) -> Self {
        Self::structure().from_base_ring(a)
    }
    fn numerator(&self) -> <<Self::Structure as FieldOfFractionsStructure>::RS as Structure>::Set {
        Self::structure().numerator(self)
    }
    fn denominator(
        &self,
    ) -> <<Self::Structure as FieldOfFractionsStructure>::RS as Structure>::Set {
        Self::structure().denominator(self)
    }
    fn as_base_ring(
        self,
    ) -> Option<<<Self::Structure as FieldOfFractionsStructure>::RS as Structure>::Set> {
        Self::structure().as_base_ring(self)
    }
}
impl<R: MetaRing> MetaFieldOfFractions for R where
    Self::Structure: FieldOfFractionsStructure<Set = R>
{
}

impl<FS: FieldOfFractionsStructure> CharZeroStructure for FS where FS::RS: CharZeroStructure {}

impl<FS: FieldOfFractionsStructure> ComplexSubsetStructure for FS where
    FS::RS: ComplexSubsetStructure
{
}

impl<FS: FieldOfFractionsStructure> RealSubsetStructure for FS where FS::RS: RealSubsetStructure {}

impl<FS: FieldOfFractionsStructure> RealToFloatStructure for FS
where
    FS::RS: RealToFloatStructure,
{
    fn as_f64(&self, x: &Self::Set) -> f64 {
        let base_ring = self.base_ring_structure();
        RealToFloatStructure::as_f64(base_ring.as_ref(), &self.numerator(x))
            / RealToFloatStructure::as_f64(base_ring.as_ref(), &self.denominator(x))
    }
}

//is a subset of the complex numbers
pub trait ComplexSubsetStructure: Structure {}

//is a subset of the real numbers
pub trait RealSubsetStructure: ComplexSubsetStructure {}

// pub trait OrderedFieldStructure: RealSubsetStructure {
//     fn cmp(&self, x: &Self::Set, y: &Self::Set) -> std::cmp::Ordering;
// }

pub trait RealRoundingStructure: RealSubsetStructure {
    fn floor(&self, x: &Self::Set) -> Integer; //round down
    fn ceil(&self, x: &Self::Set) -> Integer; //round up
    fn round(&self, x: &Self::Set) -> Integer; //round closets, either direction is fine if mid way

    /*
    fn ceil(&self, x: &Self::Set) -> Integer {
        -self.floor(&self.neg(x))
    }
    fn round(&self, x: &Self::Set) -> Integer {
        self.floor(&self.add(
            x,
            &self.from_rat(&Rational::from_integers(Integer::ONE, Integer::TWO)).unwrap(),
        ))
    }
    */
}

pub trait MetaRealRounding: MetaType
where
    Self::Structure: RealRoundingStructure,
{
    fn floor(&self) -> Integer {
        Self::structure().floor(self)
    }

    fn ceil(&self) -> Integer {
        Self::structure().ceil(self)
    }

    fn round(&self) -> Integer {
        Self::structure().round(self)
    }
}
impl<R: MetaType> MetaRealRounding for R where Self::Structure: RealRoundingStructure<Set = R> {}

pub trait RealToFloatStructure: RealSubsetStructure {
    fn as_f64(&self, x: &Self::Set) -> f64;
    fn as_f32(&self, x: &Self::Set) -> f32 {
        RealToFloatStructure::as_f64(self, x) as f32
    }
}

pub trait MetaRealToFloat: MetaType
where
    Self::Structure: RealToFloatStructure,
{
    fn as_f64(&self) -> f64 {
        Self::structure().as_f64(self)
    }

    fn as_f32(&self) -> f32 {
        Self::structure().as_f32(self)
    }
}
impl<R: MetaType> MetaRealToFloat for R where Self::Structure: RealToFloatStructure {}

pub trait RealFromFloatStructure: RealSubsetStructure {
    fn from_f64_approx(&self, x: f64) -> Self::Set;
    fn from_f32_approx(&self, x: f32) -> Self::Set {
        self.from_f64_approx(x as f64)
    }
}

pub trait MetaRealFromFloat: MetaType
where
    Self::Structure: RealFromFloatStructure,
{
    fn from_f64_approx(x: f64) -> Self {
        Self::structure().from_f64_approx(x)
    }

    fn from_f32_approx(x: f32) -> Self {
        Self::structure().from_f32_approx(x)
    }
}
impl<R: MetaType> MetaRealFromFloat for R where Self::Structure: RealFromFloatStructure<Set = R> {}

pub trait ComplexConjugateStructure: Structure {
    fn conjugate(&self, x: &Self::Set) -> Self::Set;
}
pub trait MetaComplexConjugate: MetaType
where
    Self::Structure: ComplexConjugateStructure,
{
    fn conjugate(&self) -> Self {
        Self::structure().conjugate(self)
    }
}
impl<R: MetaType> MetaComplexConjugate for R where
    Self::Structure: ComplexConjugateStructure<Set = R>
{
}

impl<RS: RealSubsetStructure> ComplexConjugateStructure for RS {
    fn conjugate(&self, x: &Self::Set) -> Self::Set {
        x.clone()
    }
}

pub trait PositiveRealNthRootStructure: ComplexSubsetStructure {
    //if x is a non-negative real number, return the nth root of x
    //may also return Ok for other well-defined values such as for 1st root of any x and 0th root of any non-zero x, but is not required to
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()>;
}

pub trait MetaPositiveRealNthRoot: MetaType
where
    Self::Structure: PositiveRealNthRootStructure,
{
    fn nth_root(&self, n: usize) -> Result<Self, ()> {
        Self::structure().nth_root(self, n)
    }
}
impl<R: MetaType> MetaPositiveRealNthRoot for R where
    Self::Structure: PositiveRealNthRootStructure<Set = R>
{
}

pub trait AlgebraicClosureStructure: FieldStructure
where
    PolynomialStructure<Self::BFS>:
        UniqueFactorizationStructure + Structure<Set = Polynomial<<Self::BFS as Structure>::Set>>,
{
    type BFS: FieldStructure; //base field structure

    fn base_field(&self) -> Rc<Self::BFS>;

    fn base_field_inclusion(&self, x: &<Self::BFS as Structure>::Set) -> Self::Set;

    //return None for the zero polynomial
    fn all_roots_list(
        &self,
        poly: &Polynomial<<Self::BFS as Structure>::Set>,
    ) -> Option<Vec<Self::Set>>;
    fn all_roots_unique(
        &self,
        poly: &Polynomial<<Self::BFS as Structure>::Set>,
    ) -> Option<Vec<Self::Set>> {
        self.all_roots_list(
            &PolynomialStructure::new(self.base_field())
                .factor(poly)
                .unwrap()
                .expand_squarefree(),
        )
    }
    fn all_roots_powers(
        &self,
        poly: &Polynomial<<Self::BFS as Structure>::Set>,
    ) -> Option<Vec<(Self::Set, usize)>> {
        let mut root_powers = vec![];
        for (factor, k) in PolynomialStructure::new(self.base_field())
            .factor(poly)?
            .unit_and_factors()
            .1
        {
            for root in self.all_roots_list(&factor).unwrap() {
                root_powers.push((root, nat_to_usize(&k).unwrap()))
            }
        }
        Some(root_powers)
    }
}
