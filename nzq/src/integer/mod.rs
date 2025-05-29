//! The Integer type and operations.

use crate::Rational;
use crate::natural::Natural;
use crate::traits::{Abs, AbsDiff, DivMod};
use algebraeon_sets::structure::{
    CanonicalStructure, CountableSetSignature, EqSignature, MetaType, SetSignature, Signature,
    ToStringSignature,
};
use malachite_base::num::basic::traits::{One, Two, Zero};
use std::{
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub, SubAssign},
    str::FromStr,
};

/// Represent an integer {..., -2, -1, 0, 1, 2, ...}
#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, CanonicalStructure)]
pub struct Integer(malachite_nz::integer::Integer);

#[allow(clippy::wrong_self_convention)]
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

impl ToStringSignature for IntegerCanonicalStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{}", elem)
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
impl From<&Integer> for Integer {
    fn from(value: &Integer) -> Self {
        value.clone()
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

macro_rules! impl_try_into_unsigned {
    ($($t:ty),*) => {
        $(
            impl TryInto<$t> for Integer {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    (&self).try_into()
                }
            }
            impl TryInto<$t> for &Integer {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    if self < &Integer::ZERO {
                        Err(())
                    } else {
                        if self > &Integer::from(<$t>::MAX) {
                            Err(())
                        } else {
                            let limbs = self.to_malachite_ref().to_twos_complement_limbs_asc();
                            match limbs.len() {
                                0 => {
                                    Ok(0)
                                }
                                1 => {
                                    Ok(limbs[0] as $t)
                                }
                                2 | 3 => {
                                    if limbs.len() == 3 {
                                        debug_assert!(limbs[2] == 0); // malachite sometimes adds a 0 on the end for some reason
                                    }
                                    let low = limbs[0] as u128;
                                    let high = limbs[1] as u128;
                                    let value = (high << 64) | low;
                                    Ok(value as $t)
                                }
                                _ => {unreachable!()}
                            }
                        }
                    }
                }
            }
        )*
    };
}

macro_rules! impl_try_into_signed {
    ($($t:ty),*) => {
        $(
            impl TryInto<$t> for Integer {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    (&self).try_into()
                }
            }
            impl TryInto<$t> for &Integer {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    if self > &Integer::from(<$t>::MAX) {
                        Err(())
                    } else if self < &Integer::from(<$t>::MIN) {
                        Err(())
                    } else {
                        let limbs = self.to_malachite_ref().to_twos_complement_limbs_asc();
                        match limbs.len() {
                            0 => {
                                Ok(0)
                            }
                            1 => {
                                Ok(limbs[0] as $t)
                            }
                            2 | 3 => {
                                if limbs.len() == 3 {
                                    debug_assert!(limbs[2] == 0); // malachite sometimes adds a 0 on the end for some reason
                                }
                                let low = limbs[0] as u128;
                                let high = limbs[1] as u128;
                                let value = (high << 64) | low;
                                Ok(value as $t)
                            }
                            _ => {unreachable!()}
                        }
                    }
                }
            }
        )*
    };
}

impl_try_into_unsigned!(u8, u16, u32, u64, u128, usize);
impl_try_into_signed!(i8, i16, i32, i64, i128, isize);

#[allow(clippy::from_over_into, clippy::cast_precision_loss)]
impl Into<f64> for Integer {
    fn into(self) -> f64 {
        if self < Integer::ZERO {
            -<Self as Into<f64>>::into(-self)
        } else {
            let limbs = self.to_malachite().into_twos_complement_limbs_asc();
            let mut flt = 0.0;
            for (i, k) in limbs.into_iter().enumerate() {
                flt += (k as f64) * (2.0_f64).powf(i as f64 * 64.0);
            }
            flt
        }
    }
}
#[allow(clippy::from_over_into)]
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

    pub fn sqrt_if_square(&self) -> Option<Natural> {
        if self < &Integer::ZERO {
            None
        } else {
            self.abs().sqrt_if_square()
        }
    }

    pub fn is_square(&self) -> bool {
        if self < &Integer::ZERO {
            false
        } else {
            self.abs().is_square()
        }
    }
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
        self.0.add_assign(rhs.0);
    }
}
impl AddAssign<&Integer> for Integer {
    fn add_assign(&mut self, rhs: &Integer) {
        self.0.add_assign(&rhs.0);
    }
}

impl SubAssign<Integer> for Integer {
    fn sub_assign(&mut self, rhs: Integer) {
        self.0.sub_assign(rhs.0);
    }
}
impl SubAssign<&Integer> for Integer {
    fn sub_assign(&mut self, rhs: &Integer) {
        self.0.sub_assign(&rhs.0);
    }
}

impl MulAssign<Integer> for Integer {
    fn mul_assign(&mut self, rhs: Integer) {
        self.0.mul_assign(rhs.0);
    }
}
impl MulAssign<&Integer> for Integer {
    fn mul_assign(&mut self, rhs: &Integer) {
        self.0.mul_assign(&rhs.0);
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

impl AbsDiff<Integer> for Integer {
    type Output = Natural;

    fn abs_diff(self, rhs: Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural::from_malachite(self.0.abs_diff(rhs.0))
    }
}
impl AbsDiff<&Integer> for Integer {
    type Output = Natural;

    fn abs_diff(self, rhs: &Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural::from_malachite(self.0.abs_diff(&rhs.0))
    }
}
impl AbsDiff<Integer> for &Integer {
    type Output = Natural;

    fn abs_diff(self, rhs: Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural::from_malachite((&self.0).abs_diff(rhs.0))
    }
}
impl AbsDiff<&Integer> for &Integer {
    type Output = Natural;

    fn abs_diff(self, rhs: &Integer) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural::from_malachite((&self.0).abs_diff(&rhs.0))
    }
}

impl CountableSetSignature for IntegerCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
        use malachite_nz::integer::exhaustive::exhaustive_integers;
        exhaustive_integers().map(Integer::from_malachite)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int_to_uint() {
        // u8
        assert_eq!(
            <&Integer as TryInto<u8>>::try_into(&Integer::from_str("0").unwrap()),
            Ok(0)
        );
        assert_eq!(
            <&Integer as TryInto<u8>>::try_into(&Integer::from_str("255").unwrap()),
            Ok(255)
        );
        assert_eq!(
            <&Integer as TryInto<u8>>::try_into(&Integer::from_str("256").unwrap()),
            Err(())
        );

        // u16
        assert_eq!(
            <&Integer as TryInto<u16>>::try_into(&Integer::from_str("65535").unwrap()),
            Ok(65535)
        );
        assert_eq!(
            <&Integer as TryInto<u16>>::try_into(&Integer::from_str("65536").unwrap()),
            Err(())
        );

        // u32
        assert_eq!(
            <&Integer as TryInto<u32>>::try_into(&Integer::from_str("4294967295").unwrap()),
            Ok(4294967295)
        );
        assert_eq!(
            <&Integer as TryInto<u32>>::try_into(&Integer::from_str("4294967296").unwrap()),
            Err(())
        );

        // u64
        assert_eq!(
            <&Integer as TryInto<u64>>::try_into(
                &Integer::from_str("18446744073709551615").unwrap()
            ),
            Ok(18446744073709551615)
        );
        assert_eq!(
            <&Integer as TryInto<u64>>::try_into(
                &Integer::from_str("18446744073709551616").unwrap()
            ),
            Err(())
        );

        // u128
        assert_eq!(
            <&Integer as TryInto<u128>>::try_into(
                &Integer::from_str("340282366920938463463374607431768211455").unwrap()
            ),
            Ok(340282366920938463463374607431768211455)
        );
        assert_eq!(
            <&Integer as TryInto<u128>>::try_into(
                &Integer::from_str("340282366920938463463374607431768211456").unwrap()
            ),
            Err(())
        );
    }

    #[test]
    fn test_int_to_int() {
        // i8
        assert_eq!(
            <&Integer as TryInto<i8>>::try_into(&Integer::from_str("-129").unwrap()),
            Err(())
        );
        assert_eq!(
            <&Integer as TryInto<i8>>::try_into(&Integer::from_str("-128").unwrap()),
            Ok(-128)
        );
        assert_eq!(
            <&Integer as TryInto<i8>>::try_into(&Integer::from_str("0").unwrap()),
            Ok(0)
        );
        assert_eq!(
            <&Integer as TryInto<i8>>::try_into(&Integer::from_str("127").unwrap()),
            Ok(127)
        );
        assert_eq!(
            <&Integer as TryInto<i8>>::try_into(&Integer::from_str("128").unwrap()),
            Err(())
        );

        // i16
        assert_eq!(
            <&Integer as TryInto<i16>>::try_into(&Integer::from_str("-32769").unwrap()),
            Err(())
        );
        assert_eq!(
            <&Integer as TryInto<i16>>::try_into(&Integer::from_str("-32768").unwrap()),
            Ok(-32768)
        );
        assert_eq!(
            <&Integer as TryInto<i16>>::try_into(&Integer::from_str("32767").unwrap()),
            Ok(32767)
        );
        assert_eq!(
            <&Integer as TryInto<i16>>::try_into(&Integer::from_str("32768").unwrap()),
            Err(())
        );

        // i32
        assert_eq!(
            <&Integer as TryInto<i32>>::try_into(&Integer::from_str("-2147483649").unwrap()),
            Err(())
        );
        assert_eq!(
            <&Integer as TryInto<i32>>::try_into(&Integer::from_str("-2147483648").unwrap()),
            Ok(-2147483648)
        );
        assert_eq!(
            <&Integer as TryInto<i32>>::try_into(&Integer::from_str("2147483647").unwrap()),
            Ok(2147483647)
        );
        assert_eq!(
            <&Integer as TryInto<i32>>::try_into(&Integer::from_str("2147483648").unwrap()),
            Err(())
        );

        // i64
        assert_eq!(
            <&Integer as TryInto<i64>>::try_into(
                &Integer::from_str("-9223372036854775809").unwrap()
            ),
            Err(())
        );
        assert_eq!(
            <&Integer as TryInto<i64>>::try_into(
                &Integer::from_str("-9223372036854775808").unwrap()
            ),
            Ok(-9223372036854775808)
        );
        assert_eq!(
            <&Integer as TryInto<i64>>::try_into(
                &Integer::from_str("9223372036854775807").unwrap()
            ),
            Ok(9223372036854775807)
        );
        assert_eq!(
            <&Integer as TryInto<i64>>::try_into(
                &Integer::from_str("9223372036854775808").unwrap()
            ),
            Err(())
        );

        // i128
        assert_eq!(
            <&Integer as TryInto<i128>>::try_into(
                &Integer::from_str("-170141183460469231731687303715884105729").unwrap()
            ),
            Err(())
        );
        assert_eq!(
            <&Integer as TryInto<i128>>::try_into(
                &Integer::from_str("-170141183460469231731687303715884105728").unwrap()
            ),
            Ok(-170141183460469231731687303715884105728)
        );
        assert_eq!(
            <&Integer as TryInto<i128>>::try_into(
                &Integer::from_str("170141183460469231731687303715884105727").unwrap()
            ),
            Ok(170141183460469231731687303715884105727)
        );
        assert_eq!(
            <&Integer as TryInto<i128>>::try_into(
                &Integer::from_str("170141183460469231731687303715884105728").unwrap()
            ),
            Err(())
        );
    }
}
