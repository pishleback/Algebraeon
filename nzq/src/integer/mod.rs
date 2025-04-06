//! The Integer type and operations.

use crate::Rational;
use crate::natural::*;
use crate::traits::*;
use algebraeon_sets::structure::*;
use malachite_base::num::basic::traits::{One, Two, Zero};
use std::{
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    str::FromStr,
};

/// Represent an integer {..., -2, -1, 0, 1, 2, ...}
#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

impl TryFrom<Rational> for Integer {
    type Error = ();
    fn try_from(value: Rational) -> Result<Self, Self::Error> {
        if let Ok(value) = malachite_nz::integer::Integer::try_from(value.to_malachite()) {
            Ok(Integer::from_malachite(value))
        } else {
            Err(())
        }
    }
}
impl TryFrom<&Rational> for Integer {
    type Error = ();
    fn try_from(value: &Rational) -> Result<Self, Self::Error> {
        if let Ok(value) = malachite_nz::integer::Integer::try_from(value.to_malachite_ref()) {
            Ok(Integer::from_malachite(value))
        } else {
            Err(())
        }
    }
}

impl Into<f64> for Integer {
    fn into(self) -> f64 {
        if self < Integer::ZERO {
            -<Self as Into<f64>>::into(-self)
        } else {
            let limbs = self.to_malachite().into_twos_complement_limbs_asc();
            let mut flt = 0.0;
            for (i, k) in limbs.into_iter().enumerate() {
                flt += (k as f64) * (2.0 as f64).powf(i as f64 * 64.0);
            }
            flt
        }
    }
}
impl Into<f64> for &Integer {
    fn into(self) -> f64 {
        self.clone().into()
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

impl PartialEq<Natural> for Integer {
    fn eq(&self, other: &Natural) -> bool {
        self.0.eq(other.to_malachite_ref())
    }
}
impl PartialEq<&Natural> for Integer {
    fn eq(&self, other: &&Natural) -> bool {
        self.eq(*other)
    }
}
impl PartialOrd<Natural> for Integer {
    fn partial_cmp(&self, other: &Natural) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(other.to_malachite_ref())
    }
}
impl PartialOrd<&Natural> for Integer {
    fn partial_cmp(&self, other: &&Natural) -> Option<std::cmp::Ordering> {
        self.partial_cmp(*other)
    }
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
        use malachite_base::num::arithmetic::traits::Mod;
        Integer(self.0.mod_op(rhs.0))
    }
}
impl Rem<&Integer> for Integer {
    type Output = Integer;

    fn rem(self, rhs: &Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::Mod;
        Integer(self.0.mod_op(&rhs.0))
    }
}
impl Rem<Integer> for &Integer {
    type Output = Integer;

    fn rem(self, rhs: Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::Mod;
        Integer((&self.0).mod_op(rhs.0))
    }
}
impl Rem<&Integer> for &Integer {
    type Output = Integer;

    fn rem(self, rhs: &Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::Mod;
        Integer((&self.0).mod_op(&rhs.0))
    }
}

impl Rem<Natural> for Integer {
    type Output = Natural;

    fn rem(self, rhs: Natural) -> Self::Output {
        self.rem(Integer::from(rhs)).abs()
    }
}
impl Rem<&Natural> for Integer {
    type Output = Natural;

    fn rem(self, rhs: &Natural) -> Self::Output {
        self.rem(Integer::from(rhs)).abs()
    }
}
impl Rem<Natural> for &Integer {
    type Output = Natural;

    fn rem(self, rhs: Natural) -> Self::Output {
        self.rem(Integer::from(rhs)).abs()
    }
}
impl Rem<&Natural> for &Integer {
    type Output = Natural;

    fn rem(self, rhs: &Natural) -> Self::Output {
        self.rem(Integer::from(rhs)).abs()
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
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<&Integer> for Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: &Integer) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(&other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<Integer> for &Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: Integer) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = (&self.0).div_mod(other.0);
        (Integer(q), Integer(r))
    }
}
impl DivMod<&Integer> for &Integer {
    type DivOutput = Integer;

    type ModOutput = Integer;

    fn div_mod(self, other: &Integer) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = (&self.0).div_mod(&other.0);
        (Integer(q), Integer(r))
    }
}

impl Abs for Integer {
    type Output = Natural;
    fn abs(self) -> Self::Output {
        use malachite_base::num::arithmetic::traits::UnsignedAbs;
        Natural::from_malachite(self.0.unsigned_abs())
    }
}

impl Abs for &Integer {
    type Output = Natural;
    fn abs(self) -> Self::Output {
        use malachite_base::num::arithmetic::traits::UnsignedAbs;
        Natural::from_malachite((&self.0).unsigned_abs())
    }
}

impl MetaType for Integer {
    type Structure = CannonicalStructure<Integer>;

    fn structure() -> std::rc::Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}
