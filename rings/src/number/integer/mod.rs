use factor::factor;
use itertools::Itertools;
use malachite_base::num::arithmetic::traits::DivMod;
use malachite_base::num::arithmetic::traits::UnsignedAbs;
use malachite_base::num::basic::traits::{One, Two, Zero};

use std::{
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    str::FromStr,
};

use crate::number::natural::*;
use crate::structure::factorization::*;
use crate::structure::structure::*;

use algebraeon_sets::structure::*;

pub mod berlekamp_zassenhaus;
pub mod modulo;
pub mod polynomial;
pub mod zimmermann_polys;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Integer(malachite_nz::integer::Integer);

impl Integer {
    pub(crate) fn from_malachite(value: malachite_nz::integer::Integer) -> Self {
        Self(value)
    }
    pub(crate) fn to_malachite(self) -> malachite_nz::integer::Integer {
        self.0
    }
    pub(crate) fn to_malachite_ref(&self) -> &malachite_nz::integer::Integer {
        &self.0
    }
}

impl std::fmt::Display for Integer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl From<u8> for Integer {
    fn from(value: u8) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<u16> for Integer {
    fn from(value: u16) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<u32> for Integer {
    fn from(value: u32) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<u64> for Integer {
    fn from(value: u64) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<u128> for Integer {
    fn from(value: u128) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<usize> for Integer {
    fn from(value: usize) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<i8> for Integer {
    fn from(value: i8) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<i16> for Integer {
    fn from(value: i16) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<i32> for Integer {
    fn from(value: i32) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<i64> for Integer {
    fn from(value: i64) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<i128> for Integer {
    fn from(value: i128) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<isize> for Integer {
    fn from(value: isize) -> Self {
        Self(malachite_nz::integer::Integer::from(value))
    }
}
impl From<Natural> for Integer {
    fn from(value: Natural) -> Self {
        Self(malachite_nz::integer::Integer::from(value.to_malachite()))
    }
}
impl From<&Natural> for Integer {
    fn from(value: &Natural) -> Self {
        Self(malachite_nz::integer::Integer::from(
            value.to_malachite_ref(),
        ))
    }
}

impl FromStr for Integer {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(malachite_nz::integer::Integer::from_str(s)?))
    }
}

impl Integer {
    pub const ZERO: Self = Self(malachite_nz::integer::Integer::ZERO);
    pub const ONE: Self = Self(malachite_nz::integer::Integer::ONE);
    pub const TWO: Self = Self(malachite_nz::integer::Integer::TWO);
}

impl AddAssign<Integer> for Integer {
    fn add_assign(&mut self, rhs: Integer) {
        self.0.add_assign(rhs.0)
    }
}
impl AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, rhs: &Integer) {
        self.0.add_assign(&rhs.0)
    }
}

impl SubAssign<Integer> for Integer {
    fn sub_assign(&mut self, rhs: Integer) {
        self.0.sub_assign(rhs.0)
    }
}
impl SubAssign<&Integer> for Integer {
    fn sub_assign(&mut self, rhs: &Integer) {
        self.0.sub_assign(&rhs.0)
    }
}

impl MulAssign<Integer> for Integer {
    fn mul_assign(&mut self, rhs: Integer) {
        self.0.mul_assign(rhs.0)
    }
}
impl MulAssign<&Integer> for Integer {
    fn mul_assign(&mut self, rhs: &Integer) {
        self.0.mul_assign(&rhs.0)
    }
}

impl Neg for Integer {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer(self.0.neg())
    }
}
impl Neg for &Integer {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer((&self.0).neg())
    }
}

impl Add<Integer> for Integer {
    type Output = Integer;

    fn add(self, rhs: Integer) -> Self::Output {
        Integer(self.0.add(rhs.0))
    }
}
impl Add<&Integer> for Integer {
    type Output = Integer;

    fn add(self, rhs: &Integer) -> Self::Output {
        Integer(self.0.add(&rhs.0))
    }
}
impl Add<Integer> for &Integer {
    type Output = Integer;

    fn add(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).add(rhs.0))
    }
}
impl Add<&Integer> for &Integer {
    type Output = Integer;

    fn add(self, rhs: &Integer) -> Self::Output {
        Integer((&self.0).add(&rhs.0))
    }
}

impl Sub<Integer> for Integer {
    type Output = Integer;

    fn sub(self, rhs: Integer) -> Self::Output {
        Integer(self.0.sub(rhs.0))
    }
}
impl Sub<&Integer> for Integer {
    type Output = Integer;

    fn sub(self, rhs: &Integer) -> Self::Output {
        Integer(self.0.sub(&rhs.0))
    }
}
impl Sub<Integer> for &Integer {
    type Output = Integer;

    fn sub(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).sub(rhs.0))
    }
}
impl Sub<&Integer> for &Integer {
    type Output = Integer;

    fn sub(self, rhs: &Integer) -> Self::Output {
        Integer((&self.0).sub(&rhs.0))
    }
}

impl Mul<Integer> for Integer {
    type Output = Integer;

    fn mul(self, rhs: Integer) -> Self::Output {
        Integer(self.0.mul(rhs.0))
    }
}
impl Mul<&Integer> for Integer {
    type Output = Integer;

    fn mul(self, rhs: &Integer) -> Self::Output {
        Integer(self.0.mul(&rhs.0))
    }
}
impl Mul<Integer> for &Integer {
    type Output = Integer;

    fn mul(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).mul(rhs.0))
    }
}
impl Mul<&Integer> for &Integer {
    type Output = Integer;

    fn mul(self, rhs: &Integer) -> Self::Output {
        Integer((&self.0).mul(&rhs.0))
    }
}

impl Rem<Integer> for Integer {
    type Output = Integer;

    fn rem(self, rhs: Integer) -> Self::Output {
        Integer(self.0.rem(rhs.0))
    }
}
impl Rem<&Integer> for Integer {
    type Output = Integer;

    fn rem(self, rhs: &Integer) -> Self::Output {
        Integer(self.0.rem(&rhs.0))
    }
}
impl Rem<Integer> for &Integer {
    type Output = Integer;

    fn rem(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).rem(rhs.0))
    }
}
impl Rem<&Integer> for &Integer {
    type Output = Integer;

    fn rem(self, rhs: &Integer) -> Self::Output {
        Integer((&self.0).rem(&rhs.0))
    }
}

impl Div<Integer> for Integer {
    type Output = Integer;

    fn div(self, rhs: Integer) -> Self::Output {
        Integer(self.0.div(rhs.0))
    }
}
impl Div<&Integer> for Integer {
    type Output = Integer;

    fn div(self, rhs: &Integer) -> Self::Output {
        Integer(self.0.div(&rhs.0))
    }
}
impl Div<Integer> for &Integer {
    type Output = Integer;

    fn div(self, rhs: Integer) -> Self::Output {
        Integer((&self.0).div(rhs.0))
    }
}
impl Div<&Integer> for &Integer {
    type Output = Integer;

    fn div(self, rhs: &Integer) -> Self::Output {
        Integer((&self.0).div(&rhs.0))
    }
}

impl DivMod<Integer> for Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: Integer) -> (Self::DivOutput, Self::ModOutput) {
        let (q, r) = self.0.div_mod(other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<&Integer> for Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: &Integer) -> (Self::DivOutput, Self::ModOutput) {
        let (q, r) = self.0.div_mod(&other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<Integer> for &Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: Integer) -> (Self::DivOutput, Self::ModOutput) {
        let (q, r) = (&self.0).div_mod(other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<&Integer> for &Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: &Integer) -> (Self::DivOutput, Self::ModOutput) {
        let (q, r) = (&self.0).div_mod(&other.0);
        (Integer(q), Integer(r))
    }
}

impl Integer{
    pub fn unsigned_abs(self) -> Natural {
        Natural::from_malachite(self.0.unsigned_abs())
    }

    pub fn unsigned_abs_ref(&self) -> Natural {
        Natural::from_malachite(self.0.unsigned_abs())
    }
}

impl MetaType for Integer {
    type Structure = CannonicalStructure<Integer>;

    fn structure() -> std::rc::Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl SemiRingStructure for CannonicalStructure<Integer> {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }

    fn one(&self) -> Self::Set {
        Integer::ONE
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl RingStructure for CannonicalStructure<Integer> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }
}

impl IntegralDomainStructure for CannonicalStructure<Integer> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.quorem(a, b) {
            Some((q, r)) => {
                if r == self.zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            None => Err(RingDivisionError::DivideByZero),
        }
    }
}

impl OrderedRingStructure for CannonicalStructure<Integer> {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FiniteUnitsStructure for CannonicalStructure<Integer> {
    fn all_units(&self) -> Vec<Self::Set> {
        vec![Integer::ONE, -Integer::ONE]
    }
}

impl FavoriteAssociateStructure for CannonicalStructure<Integer> {
    fn factor_fav_assoc(&self, a: &Self::Set) -> (Self::Set, Self::Set) {
        if a == &Integer::ZERO {
            (Integer::ONE, Integer::ZERO)
        } else if a < &Integer::ZERO {
            (-Integer::ONE, -a)
        } else {
            (Integer::ONE, a.clone())
        }
    }
}

impl UniqueFactorizationStructure for CannonicalStructure<Integer> {
    fn factor(&self, a: &Self::Set) -> Option<Factored<Self>> {
        if a == &Integer::ZERO {
            None
        } else {
            let unit;
            if a < &Integer::ZERO {
                unit = Integer::from(-1);
            } else {
                unit = Integer::from(1);
            }
            let f = factor(a.unsigned_abs_ref()).unwrap();
            Some(Factored::new_unchecked(
                self.clone().into(),
                unit,
                f.into_powers()
                    .into_iter()
                    .map(|(p, k)| (Integer::from(p), k))
                    .collect(),
            ))
        }
    }
}

impl EuclideanDivisionStructure for CannonicalStructure<Integer> {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        if elem == &Integer::ZERO {
            None
        } else {
            Some(elem.unsigned_abs_ref())
        }
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        if b == &Integer::ZERO {
            None
        } else {
            Some(a.div_mod(b.clone()))
        }
    }
}

impl GreatestCommonDivisorStructure for CannonicalStructure<Integer> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        Integer::structure().euclidean_gcd(x.clone(), y.clone())
    }
}

impl BezoutDomainStructure for CannonicalStructure<Integer> {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        Integer::euclidean_xgcd(x.clone(), y.clone())
    }
}

impl CharZeroStructure for CannonicalStructure<Integer> {}

impl ComplexSubsetStructure for CannonicalStructure<Integer> {}

impl RealSubsetStructure for CannonicalStructure<Integer> {}

impl RealToFloatStructure for CannonicalStructure<Integer> {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        if x < &Integer::ZERO {
            -self.as_f64(&-x)
        } else {
            let limbs = x.clone().to_malachite().into_twos_complement_limbs_asc();
            let mut flt = 0.0;
            for (i, k) in limbs.into_iter().enumerate() {
                flt += (k as f64) * (2.0 as f64).powf(i as f64 * 64.0);
            }
            flt
        }
    }
}
