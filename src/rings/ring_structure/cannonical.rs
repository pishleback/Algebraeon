use std::fmt::Display;
use std::rc::Rc;
use std::{borrow::Borrow, fmt::Debug, marker::PhantomData};

use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use super::super::structure::*;
use super::elements::*;
use super::factorization::*;
use super::structure::*;

pub trait Ring: StructuredType
where
    Self::Structure: RingStructure<Set = Self>,
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
    fn sum(vals: Vec<&Self>) -> Self {
        Self::structure().sum(vals)
    }
    fn mul(a: &Self, b: &Self) -> Self {
        Self::structure().mul(a, b)
    }
    fn product(vals: Vec<&Self>) -> Self {
        Self::structure().product(vals)
    }
    fn div(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::structure().div(a, b)
    }
    fn inv(&self) -> Result<Self, RingDivisionError> {
        Self::structure().inv(self)
    }

    fn nat_pow(a: &Self, n: &Natural) -> Self {
        Self::structure().nat_pow(a, n)
    }
    fn int_pow(a: &Self, n: &Integer) -> Option<Self> {
        Self::structure().int_pow(a, n)
    }
    fn from_int(x: &Integer) -> Self {
        Self::structure().from_int(x)
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

//TODO: autogenerate these cannonical versions of the ring structures and their implementations in a proc macro?
//In a proc macro how would I do the ones which go from opp(a : &Self::Set) to opp(&self)?

impl<R: StructuredType> Ring for R where Self::Structure: RingStructure<Set = R> {}

pub trait IntegralDomain: Ring
where
    Self::Structure: RingStructure<Set = Self>,
{
}

impl<R: Ring> IntegralDomain for R where Self::Structure: IntegralDomainStructure<Set = R> {}

pub trait FiniteUnitsDomain: Ring
where
    Self::Structure: FiniteUnitsStructure<Set = Self>,
{
    fn all_units(&self) -> Vec<Self> {
        Self::structure().all_units()
    }
}

impl<R: Ring> FiniteUnitsDomain for R where Self::Structure: FiniteUnitsStructure<Set = R> {}

pub trait FavoriteAssociateDomain: IntegralDomain
where
    Self::Structure: FavoriteAssociateStructure<Set = Self>,
{
    fn factor_fav_assoc(&self) -> (Self, Self) {
        Self::structure().factor_fav_assoc(self)
    }

    fn fav_assoc(&self) -> Self {
        Self::structure().fav_assoc(self)
    }
}

impl<R: Ring> FavoriteAssociateDomain for R where
    Self::Structure: FavoriteAssociateStructure<Set = R>
{
}

pub trait GreatestCommonDivisorDomain: FavoriteAssociateDomain
where
    Self::Structure: GreatestCommonDivisorStructure<Set = Self>,
{
    fn gcd(x: &Self, y: &Self) -> Self {
        Self::structure().gcd(x, y)
    }

    fn gcd_list(elems: Vec<&Self>) -> Self {
        Self::structure().gcd_list(elems)
    }

    fn lcm(x: &Self, y: &Self) -> Self {
        Self::structure().lcm(x, y)
    }

    fn lcm_list(elems: Vec<&Self>) -> Self {
        Self::structure().lcm_list(elems)
    }
}

impl<R: Ring> GreatestCommonDivisorDomain for R where
    Self::Structure: GreatestCommonDivisorStructure<Set = R>
{
}

pub trait UniqueFactorizationDomain: FavoriteAssociateDomain
where
    Self::Structure: UniqueFactorizationStructure<Set = Self>,
{
    fn factor(&self) -> Option<Factored<Self::Structure>> {
        Self::structure().factor(self)
    }

    fn is_irreducible(&self) -> bool {
        Self::structure().is_irreducible(self)
    }
}

impl<R: Ring> UniqueFactorizationDomain for R where
    Self::Structure: UniqueFactorizationStructure<Set = R>
{
}

pub trait BezoutDomain: GreatestCommonDivisorDomain
where
    Self::Structure: BezoutDomainStructure<Set = Self>,
{
    fn xgcd(x: &Self, y: &Self) -> (Self, Self, Self) {
        Self::structure().xgcd(x, y)
    }

    fn xgcd_list(elems: Vec<&Self>) -> (Self, Vec<Self>) {
        Self::structure().xgcd_list(elems)
    }
}

impl<R: Ring> BezoutDomain for R where Self::Structure: BezoutDomainStructure<Set = R> {}

pub trait EuclideanDivisionDomain: IntegralDomain
where
    Self::Structure: EuclideanDivisionStructure<Set = Self>,
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

    fn euclidean_gcd(mut x: Self, mut y: Self) -> Self
    where
        Self::Structure: FavoriteAssociateStructure,
    {
        Self::structure().euclidean_gcd(x, y)
    }

    fn euclidean_xgcd(mut x: Self, mut y: Self) -> (Self, Self, Self)
    where
        Self::Structure: GreatestCommonDivisorStructure,
    {
        Self::structure().euclidean_xgcd(x, y)
    }
}

impl<R: Ring> EuclideanDivisionDomain for R where
    Self::Structure: EuclideanDivisionStructure<Set = R>
{
}

pub trait CharZero: Ring
where
    Self::Structure: CharZeroStructure<Set = Self>,
{
}

impl<R: Ring> CharZero for R where Self::Structure: CharZeroStructure<Set = R> {}

pub trait RealDomain: CharZero
where
    Self::Structure: RealStructure<Set = Self>,
{
    fn as_f64(&self) -> f64 {
        Self::structure().as_f64(self)
    }

    fn as_f32(&self) -> f32 {
        Self::structure().as_f32(self)
    }
}

impl<R: CharZero> RealDomain for R where Self::Structure: RealStructure<Set = R> {}

pub trait DenseRealDomain: RealDomain
where
    Self::Structure: DenseRealStructure<Set = Self>,
{
    fn from_f64_approx(x: f64) -> Self {
        Self::structure().from_f64_approx(x)
    }

    fn from_f32_approx(x: f32) -> Self {
        Self::structure().from_f32_approx(x)
    }
}

impl<R: RealDomain> DenseRealDomain for R where Self::Structure: DenseRealStructure<Set = R> {}

// impl<R: Ring> RingStructure for CannonicalStructure<R> {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         a == b
//     }

//     fn zero(&self) -> Self::Set {
//         R::zero()
//     }

//     fn one(&self) -> Self::Set {
//         R::one()
//     }

//     fn neg(&self, a: &Self::Set) -> Self::Set {
//         R::neg(a)
//     }

//     fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
//         R::add(a, b)
//     }

//     fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
//         R::mul(a, b)
//     }

//     fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
//         R::div(a, b)
//     }
// }

// impl<R: Ring + Display> DisplayableStructure for CannonicalStructure<R> {
//     fn fmt_elem(&self, elem: &Self::Set, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{}", elem)
//     }
// }
