//! The Natural type and operations.

use crate::*;
use algebraeon_macros::CanonicalStructure;
use malachite::base::num::basic::traits::{One, Two, Zero};
use malachite::base::num::conversion::traits::ExactFrom;
use std::iter::{Product, Sum};
use std::{
    borrow::Borrow,
    ops::{
        Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Div, Mul,
        MulAssign, Neg, Rem, Shl, Shr, Sub, SubAssign,
    },
    str::FromStr,
};

/// Represents a natural number {0, 1, 2, ...}
#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, CanonicalStructure)]
#[canonical_structure(eq, partial_ord, ord)]
pub struct Natural(malachite::natural::Natural);

#[allow(clippy::wrong_self_convention)]
impl Natural {
    pub(crate) fn from_malachite(value: malachite::natural::Natural) -> Self {
        Self(value)
    }
    pub(crate) fn to_malachite(self) -> malachite::natural::Natural {
        self.0
    }
    pub(crate) fn to_malachite_ref(&self) -> &malachite::natural::Natural {
        &self.0
    }
}

impl ToStringSignature for NaturalCanonicalStructure {
    fn to_string(&self, elem: &Self::Elem) -> String {
        format!("{}", elem)
    }
}

impl std::fmt::Display for Natural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl From<u8> for Natural {
    fn from(value: u8) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<u16> for Natural {
    fn from(value: u16) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<u32> for Natural {
    fn from(value: u32) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<u64> for Natural {
    fn from(value: u64) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<u128> for Natural {
    fn from(value: u128) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<usize> for Natural {
    fn from(value: usize) -> Self {
        Self(malachite::natural::Natural::from(value))
    }
}
impl From<&Natural> for Natural {
    fn from(value: &Natural) -> Self {
        value.clone()
    }
}

impl TryFrom<Natural> for u8 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(u8::MAX) {
            return Err(());
        }
        Ok(u8::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for u16 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(u16::MAX) {
            return Err(());
        }
        Ok(u16::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for u32 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(u32::MAX) {
            return Err(());
        }
        Ok(u32::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for u64 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(u64::MAX) {
            return Err(());
        }
        Ok(u64::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for u128 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(u128::MAX) {
            return Err(());
        }
        Ok(u128::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for usize {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(usize::MAX) {
            return Err(());
        }
        Ok(usize::exact_from(&value.0))
    }
}

impl TryFrom<Natural> for i8 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(i8::MAX as u8) {
            return Err(());
        }
        Ok(i8::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for i16 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(i16::MAX as u16) {
            return Err(());
        }
        Ok(i16::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for i32 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(i32::MAX as u32) {
            return Err(());
        }
        Ok(i32::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for i64 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(i64::MAX as u64) {
            return Err(());
        }
        Ok(i64::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for i128 {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(i128::MAX as u128) {
            return Err(());
        }
        Ok(i128::exact_from(&value.0))
    }
}
impl TryFrom<Natural> for isize {
    type Error = ();
    fn try_from(value: Natural) -> Result<Self, Self::Error> {
        if value > Natural::from(isize::MAX as usize) {
            return Err(());
        }
        Ok(isize::exact_from(&value.0))
    }
}

impl TryFrom<&Natural> for u8 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(u8::MAX) {
            return Err(());
        }
        Ok(u8::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for u16 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(u16::MAX) {
            return Err(());
        }
        Ok(u16::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for u32 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(u32::MAX) {
            return Err(());
        }
        Ok(u32::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for u64 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(u64::MAX) {
            return Err(());
        }
        Ok(u64::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for u128 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(u128::MAX) {
            return Err(());
        }
        Ok(u128::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for usize {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(usize::MAX) {
            return Err(());
        }
        Ok(usize::exact_from(&value.0))
    }
}

impl TryFrom<&Natural> for i8 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(i8::MAX as u8) {
            return Err(());
        }
        Ok(i8::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for i16 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(i16::MAX as u16) {
            return Err(());
        }
        Ok(i16::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for i32 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(i32::MAX as u32) {
            return Err(());
        }
        Ok(i32::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for i64 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(i64::MAX as u64) {
            return Err(());
        }
        Ok(i64::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for i128 {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(i128::MAX as u128) {
            return Err(());
        }
        Ok(i128::exact_from(&value.0))
    }
}
impl TryFrom<&Natural> for isize {
    type Error = ();
    fn try_from(value: &Natural) -> Result<Self, Self::Error> {
        if value > &Natural::from(isize::MAX as usize) {
            return Err(());
        }
        Ok(isize::exact_from(&value.0))
    }
}

impl FromStr for Natural {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(malachite::natural::Natural::from_str(s)?))
    }
}

impl TryFrom<Integer> for Natural {
    type Error = ();

    fn try_from(value: Integer) -> Result<Self, Self::Error> {
        Ok(Self::from_malachite(
            malachite::natural::Natural::try_from(value.to_malachite()).map_err(|_| ())?,
        ))
    }
}

impl TryFrom<&Integer> for Natural {
    type Error = ();

    fn try_from(value: &Integer) -> Result<Self, Self::Error> {
        Ok(Self::from_malachite(
            malachite::natural::Natural::try_from(value.to_malachite_ref()).map_err(|_| ())?,
        ))
    }
}

impl TryFrom<Rational> for Natural {
    type Error = ();

    fn try_from(value: Rational) -> Result<Self, Self::Error> {
        Ok(Self::from_malachite(
            malachite::natural::Natural::try_from(value.to_malachite()).map_err(|_| ())?,
        ))
    }
}

impl TryFrom<&Rational> for Natural {
    type Error = ();

    fn try_from(value: &Rational) -> Result<Self, Self::Error> {
        Ok(Self::from_malachite(
            malachite::natural::Natural::try_from(value.to_malachite_ref()).map_err(|_| ())?,
        ))
    }
}

impl Natural {
    pub const ZERO: Self = Self(malachite::natural::Natural::ZERO);
    pub const ONE: Self = Self(malachite::natural::Natural::ONE);
    pub const TWO: Self = Self(malachite::natural::Natural::TWO);

    pub fn latex(&self) -> String {
        format!("{}", self)
    }

    pub fn typst(&self) -> String {
        format!("{}", self)
    }
}

impl AddAssign<Natural> for Natural {
    fn add_assign(&mut self, rhs: Natural) {
        self.0.add_assign(rhs.0);
    }
}
impl AddAssign<&Natural> for Natural {
    fn add_assign(&mut self, rhs: &Natural) {
        self.0.add_assign(&rhs.0);
    }
}

impl SubAssign<Natural> for Natural {
    fn sub_assign(&mut self, rhs: Natural) {
        self.0.sub_assign(rhs.0);
    }
}
impl SubAssign<&Natural> for Natural {
    fn sub_assign(&mut self, rhs: &Natural) {
        self.0.sub_assign(&rhs.0);
    }
}

impl MulAssign<Natural> for Natural {
    fn mul_assign(&mut self, rhs: Natural) {
        self.0.mul_assign(rhs.0);
    }
}
impl MulAssign<&Natural> for Natural {
    fn mul_assign(&mut self, rhs: &Natural) {
        self.0.mul_assign(&rhs.0);
    }
}

impl BitAndAssign<Natural> for Natural {
    fn bitand_assign(&mut self, rhs: Natural) {
        self.0.bitand_assign(rhs.0);
    }
}
impl BitAndAssign<&Natural> for Natural {
    fn bitand_assign(&mut self, rhs: &Natural) {
        self.0.bitand_assign(&rhs.0);
    }
}

impl BitOrAssign<Natural> for Natural {
    fn bitor_assign(&mut self, rhs: Natural) {
        self.0.bitor_assign(rhs.0);
    }
}
impl BitOrAssign<&Natural> for Natural {
    fn bitor_assign(&mut self, rhs: &Natural) {
        self.0.bitor_assign(&rhs.0);
    }
}

impl BitXorAssign<Natural> for Natural {
    fn bitxor_assign(&mut self, rhs: Natural) {
        self.0.bitxor_assign(rhs.0);
    }
}
impl BitXorAssign<&Natural> for Natural {
    fn bitxor_assign(&mut self, rhs: &Natural) {
        self.0.bitxor_assign(&rhs.0);
    }
}

impl Neg for Natural {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer::from_malachite(self.0.neg())
    }
}
impl Neg for &Natural {
    type Output = Integer;

    fn neg(self) -> Self::Output {
        Integer::from_malachite((&self.0).neg())
    }
}

impl Add<Natural> for Natural {
    type Output = Natural;

    fn add(self, rhs: Natural) -> Self::Output {
        Natural(self.0.add(rhs.0))
    }
}
impl Add<&Natural> for Natural {
    type Output = Natural;

    fn add(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.add(&rhs.0))
    }
}
impl Add<Natural> for &Natural {
    type Output = Natural;

    fn add(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).add(rhs.0))
    }
}
impl Add<&Natural> for &Natural {
    type Output = Natural;

    fn add(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).add(&rhs.0))
    }
}

impl Sub<Natural> for Natural {
    type Output = Natural;

    fn sub(self, rhs: Natural) -> Self::Output {
        Natural(self.0.sub(rhs.0))
    }
}
impl Sub<&Natural> for Natural {
    type Output = Natural;

    fn sub(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.sub(&rhs.0))
    }
}
impl Sub<Natural> for &Natural {
    type Output = Natural;

    fn sub(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).sub(rhs.0))
    }
}
impl Sub<&Natural> for &Natural {
    type Output = Natural;

    fn sub(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).sub(&rhs.0))
    }
}

impl Mul<Natural> for Natural {
    type Output = Natural;

    fn mul(self, rhs: Natural) -> Self::Output {
        Natural(self.0.mul(rhs.0))
    }
}
impl Mul<&Natural> for Natural {
    type Output = Natural;

    fn mul(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.mul(&rhs.0))
    }
}
impl Mul<Natural> for &Natural {
    type Output = Natural;

    fn mul(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).mul(rhs.0))
    }
}
impl Mul<&Natural> for &Natural {
    type Output = Natural;

    fn mul(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).mul(&rhs.0))
    }
}

impl Rem<Natural> for Natural {
    type Output = Natural;

    fn rem(self, rhs: Natural) -> Self::Output {
        Natural(self.0.rem(rhs.0))
    }
}
impl Rem<&Natural> for Natural {
    type Output = Natural;

    fn rem(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.rem(&rhs.0))
    }
}
impl Rem<Natural> for &Natural {
    type Output = Natural;

    fn rem(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).rem(rhs.0))
    }
}
impl Rem<&Natural> for &Natural {
    type Output = Natural;

    fn rem(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).rem(&rhs.0))
    }
}

impl Div<Natural> for Natural {
    type Output = Natural;

    fn div(self, rhs: Natural) -> Self::Output {
        Natural(self.0.div(rhs.0))
    }
}
impl Div<&Natural> for Natural {
    type Output = Natural;

    fn div(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.div(&rhs.0))
    }
}
impl Div<Natural> for &Natural {
    type Output = Natural;

    fn div(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).div(rhs.0))
    }
}
impl Div<&Natural> for &Natural {
    type Output = Natural;

    fn div(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).div(&rhs.0))
    }
}

impl DivMod<Natural> for Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite::base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<&Natural> for Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: &Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite::base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(&other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<Natural> for &Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite::base::num::arithmetic::traits::DivMod;
        let (q, r) = (&self.0).div_mod(other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<&Natural> for &Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: &Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite::base::num::arithmetic::traits::DivMod;
        let (q, r) = (&self.0).div_mod(&other.0);
        (Natural(q), Natural(r))
    }
}

impl BitAnd<Natural> for Natural {
    type Output = Natural;

    fn bitand(self, rhs: Natural) -> Self::Output {
        Natural(self.0.bitand(rhs.0))
    }
}
impl BitAnd<&Natural> for Natural {
    type Output = Natural;

    fn bitand(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.bitand(&rhs.0))
    }
}
impl BitAnd<Natural> for &Natural {
    type Output = Natural;

    fn bitand(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).bitand(rhs.0))
    }
}
impl BitAnd<&Natural> for &Natural {
    type Output = Natural;

    fn bitand(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).bitand(&rhs.0))
    }
}

impl BitOr<Natural> for Natural {
    type Output = Natural;

    fn bitor(self, rhs: Natural) -> Self::Output {
        Natural(self.0.bitor(rhs.0))
    }
}
impl BitOr<&Natural> for Natural {
    type Output = Natural;

    fn bitor(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.bitor(&rhs.0))
    }
}
impl BitOr<Natural> for &Natural {
    type Output = Natural;

    fn bitor(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).bitor(rhs.0))
    }
}
impl BitOr<&Natural> for &Natural {
    type Output = Natural;

    fn bitor(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).bitor(&rhs.0))
    }
}

impl BitXor<Natural> for Natural {
    type Output = Natural;

    fn bitxor(self, rhs: Natural) -> Self::Output {
        Natural(self.0.bitxor(rhs.0))
    }
}
impl BitXor<&Natural> for Natural {
    type Output = Natural;

    fn bitxor(self, rhs: &Natural) -> Self::Output {
        Natural(self.0.bitxor(&rhs.0))
    }
}
impl BitXor<Natural> for &Natural {
    type Output = Natural;

    fn bitxor(self, rhs: Natural) -> Self::Output {
        Natural((&self.0).bitxor(rhs.0))
    }
}
impl BitXor<&Natural> for &Natural {
    type Output = Natural;

    fn bitxor(self, rhs: &Natural) -> Self::Output {
        Natural((&self.0).bitxor(&rhs.0))
    }
}

impl<T> Shl<T> for Natural
where
    malachite::natural::Natural: Shl<T, Output = malachite::natural::Natural>,
{
    type Output = Natural;

    fn shl(self, rhs: T) -> Self::Output {
        Natural(self.0.shl(rhs))
    }
}
impl<T> Shl<T> for &Natural
where
    for<'a> &'a malachite::natural::Natural: Shl<T, Output = malachite::natural::Natural>,
{
    type Output = Natural;

    fn shl(self, rhs: T) -> Self::Output {
        Natural((&self.0).shl(rhs))
    }
}

impl<T> Shr<T> for Natural
where
    malachite::natural::Natural: Shr<T, Output = malachite::natural::Natural>,
{
    type Output = Natural;

    fn shr(self, rhs: T) -> Self::Output {
        Natural(self.0.shr(rhs))
    }
}
impl<T> Shr<T> for &Natural
where
    for<'a> &'a malachite::natural::Natural: Shr<T, Output = malachite::natural::Natural>,
{
    type Output = Natural;

    fn shr(self, rhs: T) -> Self::Output {
        Natural((&self.0).shr(rhs))
    }
}

impl AbsDiff<Natural> for Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: Natural) -> Self::Output {
        use malachite::base::num::arithmetic::traits::AbsDiff;
        Natural(self.0.abs_diff(rhs.0))
    }
}
impl AbsDiff<&Natural> for Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: &Natural) -> Self::Output {
        use malachite::base::num::arithmetic::traits::AbsDiff;
        Natural(self.0.abs_diff(&rhs.0))
    }
}
impl AbsDiff<Natural> for &Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: Natural) -> Self::Output {
        use malachite::base::num::arithmetic::traits::AbsDiff;
        Natural((&self.0).abs_diff(rhs.0))
    }
}
impl AbsDiff<&Natural> for &Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: &Natural) -> Self::Output {
        use malachite::base::num::arithmetic::traits::AbsDiff;
        Natural((&self.0).abs_diff(&rhs.0))
    }
}

impl<Exponent: Borrow<Natural>, Modulus: Borrow<Natural>> ModPow<Exponent, Modulus> for Natural {
    type Output = Natural;
    fn mod_pow(self, exp: Exponent, m: Modulus) -> Self::Output {
        use malachite::base::num::arithmetic::traits::ModPow;
        Natural((self.0 % &m.borrow().0).mod_pow(&exp.borrow().0, &m.borrow().0))
    }
}

impl<Exponent: Borrow<Natural>, Modulus: Borrow<Natural>> ModPow<Exponent, Modulus> for &Natural {
    type Output = Natural;
    fn mod_pow(self, exp: Exponent, m: Modulus) -> Self::Output {
        use malachite::base::num::arithmetic::traits::ModPow;
        Natural((&self.0 % &m.borrow().0).mod_pow(&exp.borrow().0, &m.borrow().0))
    }
}

impl<Modulus: Borrow<Natural>> ModInv<Modulus> for Natural {
    type Output = Option<Natural>;
    fn mod_inv(self, m: Modulus) -> Self::Output {
        use malachite::base::num::arithmetic::traits::ModInverse;
        Some(Natural(
            (self.0 % &m.borrow().0).mod_inverse(&m.borrow().0)?,
        ))
    }
}

impl<Modulus: Borrow<Natural>> ModInv<Modulus> for &Natural {
    type Output = Option<Natural>;
    fn mod_inv(self, m: Modulus) -> Self::Output {
        use malachite::base::num::arithmetic::traits::ModInverse;
        Some(Natural(
            (&self.0 % &m.borrow().0).mod_inverse(&m.borrow().0)?,
        ))
    }
}

impl Sum for Natural {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_malachite(malachite::natural::Natural::sum(
            iter.map(|x| x.to_malachite()),
        ))
    }
}

impl Product for Natural {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        Self::from_malachite(malachite::natural::Natural::product(
            iter.map(|x| x.to_malachite()),
        ))
    }
}

impl Natural {
    /// 2 raised to the power of `pow`.
    /// ```
    /// use algebraeon_structures::Natural;
    /// assert_eq!(
    ///     Natural::from(32u32),
    ///     Natural::power_of_2(5)
    /// );
    /// ```
    pub fn power_of_2(pow: u64) -> Self {
        use malachite::base::num::arithmetic::traits::PowerOf2;
        Self(malachite::natural::Natural::power_of_2(pow))
    }

    /// An iterator over the bits in the binary expansion.
    /// ```
    /// use algebraeon_structures::Natural;
    /// assert_eq!(
    ///     Natural::from(11u32).bits().collect::<Vec<_>>(),
    ///     vec![true, true, false, true],
    /// );
    /// ```
    pub fn bits<'a>(&'a self) -> impl ExactSizeIterator<Item = bool> + DoubleEndedIterator + 'a {
        use malachite::base::num::logic::traits::BitIterable;
        self.0.bits()
    }

    /// Return the number of bits needed to store n i.e. ceil(log2(n)) for all non-zero n
    pub fn bitcount(&self) -> usize {
        self.bits().len()
    }

    /// Return `self-other` if the result is a natural number
    pub fn try_sub(&self, other: &Self) -> Option<Self> {
        use malachite::base::num::arithmetic::traits::CheckedSub;
        Some(Self::from_malachite(
            self.to_malachite_ref()
                .checked_sub(other.to_malachite_ref())?,
        ))
    }
}

impl CountableSetSignature for NaturalCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        use malachite::natural::exhaustive::exhaustive_naturals;
        exhaustive_naturals().map(Natural::from_malachite)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_latex_and_typst() {
        assert_eq!(Natural::ZERO.latex(), "0");
        assert_eq!(Natural::ZERO.typst(), "0");

        let n = Natural::from(42u32);
        assert_eq!(n.latex(), "42");
        assert_eq!(n.typst(), "42");
    }

    #[test]
    fn test_nat_to_usize() {
        assert_eq!(
            <&Natural as TryInto<usize>>::try_into(&Natural::from(0u8)).unwrap(),
            0
        );
        assert_eq!(
            <&Natural as TryInto<usize>>::try_into(&Natural::from(1u8)).unwrap(),
            1
        );
        assert_eq!(
            <&Natural as TryInto<usize>>::try_into(&Natural::from(2u8)).unwrap(),
            2
        );
    }

    #[test]
    fn test_nat_to_uint() {
        // u8
        assert_eq!(
            <&Natural as TryInto<u8>>::try_into(&Natural::from_str("0").unwrap()),
            Ok(0)
        );
        assert_eq!(
            <&Natural as TryInto<u8>>::try_into(&Natural::from_str("255").unwrap()),
            Ok(255)
        );
        assert_eq!(
            <&Natural as TryInto<u8>>::try_into(&Natural::from_str("256").unwrap()),
            Err(())
        );

        // u16
        assert_eq!(
            <&Natural as TryInto<u16>>::try_into(&Natural::from_str("65535").unwrap()),
            Ok(65535)
        );
        assert_eq!(
            <&Natural as TryInto<u16>>::try_into(&Natural::from_str("65536").unwrap()),
            Err(())
        );

        // u32
        assert_eq!(
            <&Natural as TryInto<u32>>::try_into(&Natural::from_str("4294967295").unwrap()),
            Ok(4294967295)
        );
        assert_eq!(
            <&Natural as TryInto<u32>>::try_into(&Natural::from_str("4294967296").unwrap()),
            Err(())
        );

        // u64
        assert_eq!(
            <&Natural as TryInto<u64>>::try_into(
                &Natural::from_str("18446744073709551615").unwrap()
            ),
            Ok(18446744073709551615)
        );
        assert_eq!(
            <&Natural as TryInto<u64>>::try_into(
                &Natural::from_str("18446744073709551616").unwrap()
            ),
            Err(())
        );

        // u128
        assert_eq!(
            <&Natural as TryInto<u128>>::try_into(
                &Natural::from_str("340282366920938463463374607431768211455").unwrap()
            ),
            Ok(340282366920938463463374607431768211455)
        );
        assert_eq!(
            <&Natural as TryInto<u128>>::try_into(
                &Natural::from_str("340282366920938463463374607431768211456").unwrap()
            ),
            Err(())
        );
    }

    #[test]
    fn test_nat_to_int() {
        // i8
        assert_eq!(
            <&Natural as TryInto<i8>>::try_into(&Natural::from_str("0").unwrap()),
            Ok(0)
        );
        assert_eq!(
            <&Natural as TryInto<i8>>::try_into(&Natural::from_str("127").unwrap()),
            Ok(127)
        );
        assert_eq!(
            <&Natural as TryInto<i8>>::try_into(&Natural::from_str("128").unwrap()),
            Err(())
        );

        // i16
        assert_eq!(
            <&Natural as TryInto<i16>>::try_into(&Natural::from_str("32767").unwrap()),
            Ok(32767)
        );
        assert_eq!(
            <&Natural as TryInto<i16>>::try_into(&Natural::from_str("32768").unwrap()),
            Err(())
        );

        // i32
        assert_eq!(
            <&Natural as TryInto<i32>>::try_into(&Natural::from_str("2147483647").unwrap()),
            Ok(2147483647)
        );
        assert_eq!(
            <&Natural as TryInto<i32>>::try_into(&Natural::from_str("2147483648").unwrap()),
            Err(())
        );

        // i64
        assert_eq!(
            <&Natural as TryInto<i64>>::try_into(
                &Natural::from_str("9223372036854775807").unwrap()
            ),
            Ok(9223372036854775807)
        );
        assert_eq!(
            <&Natural as TryInto<i64>>::try_into(
                &Natural::from_str("9223372036854775808").unwrap()
            ),
            Err(())
        );

        // i128
        assert_eq!(
            <&Natural as TryInto<i128>>::try_into(
                &Natural::from_str("170141183460469231731687303715884105727").unwrap()
            ),
            Ok(170141183460469231731687303715884105727)
        );
        assert_eq!(
            <&Natural as TryInto<i128>>::try_into(
                &Natural::from_str("170141183460469231731687303715884105728").unwrap()
            ),
            Err(())
        );
    }

    #[allow(clippy::redundant_closure_for_method_calls)]
    #[test]
    fn natural_countable_list() {
        assert_eq!(
            Natural::structure()
                .generate_all_elements()
                .take(10)
                .collect::<Vec<_>>(),
            (0..10u32).map(|x| x.into()).collect::<Vec<Natural>>()
        );
    }

    #[test]
    fn natural_checked_sub() {
        assert_eq!(
            Natural::from(3u16).try_sub(&Natural::from(2u16)).unwrap(),
            Natural::from(1u16)
        );
        assert!(Natural::from(2u16).try_sub(&Natural::from(3u16)).is_none());
    }
}
