use super::*;
use crate::polynomial::*;
use algebraeon_nzq::{Integer, Natural, Rational, traits::*};
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, fmt::Debug};

#[derive(Debug, PartialEq, Eq)]
pub enum RingDivisionError {
    DivideByZero,
    NotDivisible,
}

pub trait AdditiveMonoidSignature: SetSignature {
    fn zero(&self) -> Self::Set;

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
}

pub trait AdditiveMonoidEqSignature: AdditiveMonoidSignature + EqSignature {
    fn is_zero(&self, a: &Self::Set) -> bool {
        self.equal(a, &self.zero())
    }
}
impl<R: AdditiveMonoidSignature + EqSignature> AdditiveMonoidEqSignature for R {}

pub trait MetaAdditiveMonoidEq: MetaType
where
    Self::Signature: RingSignature + EqSignature,
{
    fn is_zero(&self) -> bool {
        Self::structure().is_zero(self)
    }
}
impl<R: MetaType> MetaAdditiveMonoidEq for R where Self::Signature: RingSignature + EqSignature {}

pub trait SemiRingSignature: AdditiveMonoidSignature {
    fn one(&self) -> Self::Set;
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
        if *n == Natural::ZERO {
            self.one()
        } else if *n == Natural::ONE {
            a.clone()
        } else {
            debug_assert!(*n >= Natural::TWO);
            let bits: Vec<_> = n.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(self.mul(pows.last().unwrap(), pows.last().unwrap()));
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

    fn from_nat(&self, x: impl Into<Natural>) -> Self::Set {
        let x = x.into();
        if x == Natural::ZERO {
            self.zero()
        } else if x == Natural::ONE {
            self.one()
        } else {
            let two = self.add(&self.one(), &self.one());
            debug_assert!(x >= Natural::TWO);
            let bits: Vec<bool> = x.bits().collect();
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

pub trait MetaSemiRing: MetaType
where
    Self::Signature: SemiRingSignature,
{
    fn zero() -> Self {
        Self::structure().zero()
    }
    fn one() -> Self {
        Self::structure().one()
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

    fn from_nat(x: impl Into<Natural>) -> Self {
        Self::structure().from_nat(x)
    }
}
impl<R: MetaType> MetaSemiRing for R where Self::Signature: SemiRingSignature {}

pub trait SemiRingEqSignature: SemiRingSignature + EqSignature {}
impl<R: SemiRingSignature + EqSignature> SemiRingEqSignature for R {}

pub trait CharacteristicSignature: SemiRingSignature {
    fn characteristic(&self) -> Natural;
}
pub trait MetaCharacteristic: MetaType
where
    Self::Signature: CharacteristicSignature,
{
    fn characteristic() -> Natural {
        Self::structure().characteristic()
    }
}
impl<R: MetaType> MetaCharacteristic for R where Self::Signature: CharacteristicSignature {}

pub trait AdditiveGroupSignature: AdditiveMonoidSignature {
    fn neg(&self, a: &Self::Set) -> Self::Set;

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.add(a, &self.neg(b))
    }
}

pub trait RingSignature: SemiRingSignature + AdditiveGroupSignature {
    /// Determine whether the ring is reduced.
    ///
    /// Returns `Ok(true)` if the ring is reduced, `Ok(false)` if it is not,
    /// and `Err` when the implementation cannot decide.
    fn is_reduced(&self) -> Result<bool, String> {
        Err("unable to decide whether the ring is reduced".to_string())
    }

    fn bracket(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.sub(&self.mul(a, b), &self.mul(b, a))
    }

    fn from_int(&self, x: impl Into<Integer>) -> Self::Set {
        let x = x.into();
        if x < Integer::ZERO {
            self.neg(&self.from_int(-x))
        } else {
            self.from_nat(x.abs())
        }
    }

    fn principal_subring_inclusion<'a>(&'a self) -> PrincipalSubringInclusion<Self, &'a Self> {
        PrincipalSubringInclusion::new(self)
    }

    fn into_principal_subring_inclusion(self) -> PrincipalSubringInclusion<Self, Self> {
        PrincipalSubringInclusion::new(self)
    }
}

pub trait MetaRing: MetaType
where
    Self::Signature: RingSignature,
{
    fn ring_is_reduced() -> Result<bool, String> {
        Self::structure().is_reduced()
    }

    fn neg(&self) -> Self {
        Self::structure().neg(self)
    }

    fn from_int(x: impl Into<Integer>) -> Self {
        Self::structure().from_int(x)
    }
}
impl<R: MetaType> MetaRing for R where Self::Signature: RingSignature {}

pub trait RingEqSignature: RingSignature + EqSignature {}
impl<R: RingSignature + EqSignature> RingEqSignature for R {}

pub trait SemiRingUnitsSignature: SemiRingSignature {
    /// b such that a*b=1 and b*a=1
    /// Err(DivideByZero) if b is zero
    /// Err(NotDivisible) if no such b exists
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError>;

    fn is_unit(&self, a: &Self::Set) -> bool {
        match self.inv(a) {
            Ok(_inv) => true,
            Err(RingDivisionError::DivideByZero | RingDivisionError::NotDivisible) => false,
        }
    }
}
pub trait MetaSemiRingUnitsSignature: MetaType
where
    Self::Signature: SemiRingUnitsSignature,
{
    fn inv(&self) -> Result<Self, RingDivisionError> {
        Self::structure().inv(self)
    }
    fn is_unit(&self) -> bool {
        Self::structure().is_unit(self)
    }
}
impl<R: MetaType> MetaSemiRingUnitsSignature for R where
    Self::Signature: SemiRingUnitsSignature<Set = R>
{
}

pub trait RingUnitsSignature: RingSignature + SemiRingUnitsSignature {}
impl<Ring: RingSignature + SemiRingUnitsSignature> RingUnitsSignature for Ring {}

pub trait IntegralDomainSignature: RingUnitsSignature + EqSignature {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError>;

    fn from_rat(&self, x: &Rational) -> Option<Self::Set> {
        match self.div(
            &self.from_int(Fraction::numerator(x)),
            &self.from_nat(Fraction::denominator(x)),
        ) {
            Ok(d) => Some(d),
            Err(RingDivisionError::NotDivisible) => None,
            Err(RingDivisionError::DivideByZero) => panic!(),
        }
    }

    fn int_pow(&self, a: &Self::Set, n: &Integer) -> Option<Self::Set> {
        if *n == Integer::ZERO {
            Some(self.one())
        } else if self.is_zero(a) {
            Some(self.zero())
        } else if *n > Integer::ZERO {
            Some(self.nat_pow(a, &n.abs()))
        } else {
            match self.inv(a) {
                Ok(self_inv) => Some(self.nat_pow(&self_inv, &n.abs())),
                Err(RingDivisionError::NotDivisible) => None,
                Err(RingDivisionError::DivideByZero) => panic!(),
            }
        }
    }

    /// return true iff a is divisible by b
    fn divisible(&self, a: &Self::Set, b: &Self::Set) -> bool {
        match self.div(a, b) {
            Ok(_q) => true,
            Err(RingDivisionError::NotDivisible | RingDivisionError::DivideByZero) => false,
        }
    }
    fn are_associate(&self, a: &Self::Set, b: &Self::Set) -> bool {
        if self.equal(a, &self.zero()) && self.equal(b, &self.zero()) {
            true
        } else {
            self.div(a, b).is_ok() && self.div(b, a).is_ok()
        }
    }
}
pub trait MetaIntegralDomain: MetaRing
where
    Self::Signature: IntegralDomainSignature,
{
    fn div(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::structure().div(a, b)
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
}
impl<R: MetaRing> MetaIntegralDomain for R where Self::Signature: IntegralDomainSignature<Set = R> {}

pub trait OrderedRingSignature: IntegralDomainSignature {
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
    Self::Signature: OrderedRingSignature,
{
    fn ring_cmp(a: &Self, b: &Self) -> std::cmp::Ordering {
        Self::structure().ring_cmp(a, b)
    }

    fn abs(a: &Self) -> Self {
        Self::structure().abs(a)
    }
}
impl<R: MetaType> MetaOrderedRing for R where Self::Signature: OrderedRingSignature<Set = R> {}

pub trait FiniteUnitsSignature: RingSignature {
    fn all_units(&self) -> Vec<Self::Set>;

    fn all_units_and_zero(&self) -> Vec<Self::Set> {
        let mut elems = vec![self.zero()];
        elems.append(&mut self.all_units());
        elems
    }
}

pub trait MetaFiniteUnits: MetaRing
where
    Self::Signature: FiniteUnitsSignature,
{
    fn all_units(&self) -> Vec<Self> {
        Self::structure().all_units()
    }
}
impl<R: MetaRing> MetaFiniteUnits for R where Self::Signature: FiniteUnitsSignature<Set = R> {}

pub trait FavoriteAssociateSignature: IntegralDomainSignature {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)
    //every unit u is required to return (u, 1) i.e. 1 is the favorite associate of every unit

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
    Self::Signature: FavoriteAssociateSignature,
{
    fn factor_fav_assoc(&self) -> (Self, Self) {
        Self::structure().factor_fav_assoc(self)
    }
    fn fav_assoc(&self) -> Self {
        Self::structure().fav_assoc(self)
    }
    fn is_fav_assoc(&self) -> bool {
        Self::structure().is_fav_assoc(self)
    }
}
impl<R: MetaRing> MetaFavoriteAssociate for R where
    Self::Signature: FavoriteAssociateSignature<Set = R>
{
}

pub trait GreatestCommonDivisorSignature: FavoriteAssociateSignature {
    //any gcds should be the standard associate representative
    //euclidean_gcd can be used to implement this
    fn gcd<'a>(&'a self, x: &Self::Set, y: &Self::Set) -> Self::Set;
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
    Self::Signature: GreatestCommonDivisorSignature,
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
    Self::Signature: GreatestCommonDivisorSignature<Set = R>
{
}

pub trait BezoutDomainSignature: GreatestCommonDivisorSignature {
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
    Self::Signature: BezoutDomainSignature,
{
    fn xgcd(x: &Self, y: &Self) -> (Self, Self, Self) {
        Self::structure().xgcd(x, y)
    }

    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        Self::structure().xgcd_list(elems)
    }
}
impl<R: MetaRing> MetaBezoutDomain for R where Self::Signature: BezoutDomainSignature<Set = R> {}

pub trait EuclideanDivisionSignature: SemiRingEqSignature {
    /// None for 0 and Some(norm) for everything else
    fn norm(&self, elem: &Self::Set) -> Option<Natural>;

    /// None if b is 0 and Some((q, r)) such that a=bq+r and r=0 or norm(r) < norm(b)
    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)>;

    fn quo(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        self.quorem(a, b).map(|(q, _r)| q)
    }

    fn rem(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        if self.is_zero(b) {
            a.clone()
        } else {
            let (_q, r) = self.quorem(a, b).unwrap();
            r
        }
    }

    fn euclidean_gcd(&self, mut x: Self::Set, mut y: Self::Set) -> Self::Set
    where
        Self: FavoriteAssociateSignature,
    {
        //Euclidean algorithm
        while !self.is_zero(&y) {
            let r = self.rem(&x, &y);
            (x, y) = (y, r);
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
        Self: FavoriteAssociateSignature,
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

pub trait MetaEuclideanDivision: MetaType
where
    Self::Signature: EuclideanDivisionSignature,
{
    fn norm(&self) -> Option<Natural> {
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
        Self::Signature: FavoriteAssociateSignature,
    {
        Self::structure().euclidean_gcd(x, y)
    }

    fn euclidean_xgcd(x: Self, y: Self) -> (Self, Self, Self)
    where
        Self::Signature: GreatestCommonDivisorSignature,
    {
        Self::structure().euclidean_xgcd(x, y)
    }
}
impl<R: MetaType> MetaEuclideanDivision for R where
    Self::Signature: EuclideanDivisionSignature<Set = R>
{
}

pub trait EuclideanDomainSignature: EuclideanDivisionSignature + IntegralDomainSignature {}
impl<Ring: EuclideanDivisionSignature + IntegralDomainSignature> EuclideanDomainSignature for Ring {}

pub trait InfiniteSignature: SetSignature {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::Set>>;
}
pub trait Infinite: MetaType {
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>>;
}
impl<T: MetaType> Infinite for T
where
    T::Signature: InfiniteSignature,
{
    fn generate_distinct_elements() -> Box<dyn Iterator<Item = Self>> {
        todo!()
    }
}

pub trait FieldSignature: IntegralDomainSignature {}

impl<FS: FieldSignature> FavoriteAssociateSignature for FS {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        if self.is_zero(a) {
            (self.one(), self.zero())
        } else {
            (a.clone(), self.one())
        }
    }
}

impl<FS: FieldSignature> EuclideanDivisionSignature for FS {
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

impl<FS: FieldSignature> GreatestCommonDivisorSignature for FS {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.euclidean_gcd(x.clone(), y.clone())
    }
}

impl<FS: FieldSignature> BezoutDomainSignature for FS {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        self.euclidean_xgcd(x.clone(), y.clone())
    }
}

// A trait to indicate that characteristic() always returns 0
pub trait CharZeroRingSignature: RingSignature + CharacteristicSignature {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer>;
}
pub trait MetaCharZeroRing: MetaRing
where
    Self::Signature: CharZeroRingSignature,
{
    fn try_to_int(&self) -> Option<Integer> {
        Self::structure().try_to_int(self)
    }
}
impl<R: MetaType> MetaCharZeroRing for R where Self::Signature: CharZeroRingSignature<Set = R> {}

impl<RS: CharZeroRingSignature + 'static> InfiniteSignature for RS {
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = <Self as SetSignature>::Set>> {
        struct IntegerIterator<RS: CharZeroRingSignature> {
            ring: RS,
            next: Integer,
        }

        impl<RS: CharZeroRingSignature> Iterator for IntegerIterator<RS> {
            type Item = RS::Set;

            fn next(&mut self) -> Option<Self::Item> {
                let next = self.next.clone();
                if Integer::ZERO < next {
                    self.next = -self.next.clone();
                } else {
                    self.next = Integer::from(1) - self.next.clone();
                }
                Some(self.ring.from_int(next))
            }
        }

        Box::new(IntegerIterator {
            ring: self.clone(),
            next: Integer::from(0),
        })
    }
}

pub trait CharZeroFieldSignature: FieldSignature + CharZeroRingSignature {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational>;

    fn rational_extension<'a>(&'a self) -> PrincipalRationalSubfieldInclusion<Self, &'a Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }
    fn into_rational_extension(self) -> PrincipalRationalSubfieldInclusion<Self, Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }
}
pub trait MetaCharZeroField: MetaRing
where
    Self::Signature: CharZeroFieldSignature,
{
    fn try_to_rat(&self) -> Option<Rational> {
        Self::structure().try_to_rat(self)
    }
}
impl<R: MetaType> MetaCharZeroField for R where Self::Signature: CharZeroFieldSignature<Set = R> {}

pub trait FiniteFieldSignature: FieldSignature + FiniteUnitsSignature + FiniteSetSignature {
    // Return (p, k) where p is a prime and |F| = p^k
    fn characteristic_and_power(&self) -> (Natural, Natural);
}

//is a subset of the complex numbers
pub trait ComplexSubsetSignature: SetSignature {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32);
    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64);
}
pub trait MetaComplexSubset: MetaType
where
    Self::Signature: ComplexSubsetSignature,
{
    fn as_f32_real_and_imaginary_parts(&self) -> (f32, f32) {
        Self::structure().as_f32_real_and_imaginary_parts(self)
    }
    fn as_f64_real_and_imaginary_parts(&self) -> (f64, f64) {
        Self::structure().as_f64_real_and_imaginary_parts(self)
    }
}
impl<R: MetaType> MetaComplexSubset for R where Self::Signature: ComplexSubsetSignature<Set = R> {}

//is a subset of the real numbers
pub trait RealSubsetSignature: ComplexSubsetSignature {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        let (r, i) = self.as_f64_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }
    fn as_f32(&self, x: &Self::Set) -> f32 {
        let (r, i) = self.as_f32_real_and_imaginary_parts(x);
        debug_assert_eq!(i, 0.0);
        r
    }
}
pub trait MetaRealSubset: MetaType
where
    Self::Signature: RealSubsetSignature,
{
    fn as_f64(&self) -> f64 {
        Self::structure().as_f64(self)
    }

    fn as_f32(&self) -> f32 {
        Self::structure().as_f32(self)
    }
}
impl<R: MetaType> MetaRealSubset for R where Self::Signature: RealSubsetSignature {}

pub trait RealRoundingSignature: RealSubsetSignature {
    fn floor(&self, x: &Self::Set) -> Integer; //round down
    fn ceil(&self, x: &Self::Set) -> Integer; //round up
    fn round(&self, x: &Self::Set) -> Integer; //round closets, either direction is fine if mid way
}

pub trait MetaRealRounding: MetaType
where
    Self::Signature: RealRoundingSignature,
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
impl<R: MetaType> MetaRealRounding for R where Self::Signature: RealRoundingSignature<Set = R> {}

#[allow(clippy::wrong_self_convention)]
pub trait RealFromFloatSignature: RealSubsetSignature {
    fn from_f64_approx(&self, x: f64) -> Self::Set;
    fn from_f32_approx(&self, x: f32) -> Self::Set {
        self.from_f64_approx(f64::from(x))
    }
}

pub trait MetaRealFromFloat: MetaType
where
    Self::Signature: RealFromFloatSignature,
{
    fn from_f64_approx(x: f64) -> Self {
        Self::structure().from_f64_approx(x)
    }

    fn from_f32_approx(x: f32) -> Self {
        Self::structure().from_f32_approx(x)
    }
}
impl<R: MetaType> MetaRealFromFloat for R where Self::Signature: RealFromFloatSignature<Set = R> {}

pub trait ComplexConjugateSignature: SetSignature {
    fn conjugate(&self, x: &Self::Set) -> Self::Set;
}
pub trait MetaComplexConjugate: MetaType
where
    Self::Signature: ComplexConjugateSignature,
{
    fn conjugate(&self) -> Self {
        Self::structure().conjugate(self)
    }
}
impl<R: MetaType> MetaComplexConjugate for R where
    Self::Signature: ComplexConjugateSignature<Set = R>
{
}

impl<RS: RealSubsetSignature> ComplexConjugateSignature for RS {
    fn conjugate(&self, x: &Self::Set) -> Self::Set {
        x.clone()
    }
}

pub trait PositiveRealNthRootSignature: ComplexSubsetSignature {
    //if x is a non-negative real number, return the nth root of x
    //may also return Ok for other well-defined values such as for 1st root of any x and 0th root of any non-zero x, but is not required to
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()>;
}

pub trait MetaPositiveRealNthRoot: MetaType
where
    Self::Signature: PositiveRealNthRootSignature,
{
    fn nth_root(&self, n: usize) -> Result<Self, ()> {
        Self::structure().nth_root(self, n)
    }
}
impl<R: MetaType> MetaPositiveRealNthRoot for R where
    Self::Signature: PositiveRealNthRootSignature<Set = R>
{
}

//TODO: Move this sort of struture to the field inclusion homomorphism
pub trait AlgebraicClosureSignature: FieldSignature
where
    //TODO: can this allow polynomial structures taking a reference to the base field rather than an instance?
    PolynomialStructure<Self::BFS, Self::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<Self::BFS as SetSignature>::Set>>,
{
    type BFS: FieldSignature; //base field structure

    fn base_field(&self) -> Self::BFS;

    fn base_field_inclusion(&self, x: &<Self::BFS as SetSignature>::Set) -> Self::Set;

    //return None for the zero polynomial
    fn all_roots_list(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<Self::Set>>;

    fn all_roots_unique(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<Self::Set>> {
        let base_field_poly = self.base_field().into_polynomial_ring();
        self.all_roots_list(
            &base_field_poly
                .factorizations()
                .expanded_squarefree(&base_field_poly.factor(poly).unwrap()),
        )
    }

    fn all_roots_powers(
        &self,
        poly: &Polynomial<<Self::BFS as SetSignature>::Set>,
    ) -> Option<Vec<(Self::Set, usize)>> {
        let mut root_powers = vec![];
        let base_field_poly = self.base_field().into_polynomial_ring();
        for (factor, k) in base_field_poly
            .factorizations()
            .into_powers(base_field_poly.factor(poly)?)
        {
            for root in self.all_roots_list(&factor).unwrap() {
                root_powers.push((root, (&k).try_into().unwrap()));
            }
        }
        Some(root_powers)
    }
}

/// The free ring of rank 0 is the integers
/// The free ring of rank 1 is the polynomial ring over the integers
/// The free ring of rank n is the multipolynomial ring over the integers
pub trait FreeRingSignature: RingSignature {
    type Generator: Clone + Debug + PartialEq + Eq + std::hash::Hash + Send + Sync;

    fn free_generators(&self) -> std::collections::HashSet<Self::Generator>;
    fn free_rank(&self) -> usize {
        self.free_generators().len()
    }
}
