//! The Natural type and operations.

use crate::integer::*;
use crate::traits::*;
use algebraeon_sets::structure::*;
use malachite_base::num::{
    arithmetic::traits::PowerOf2,
    basic::traits::{One, Two, Zero},
};
use std::{
    borrow::Borrow,
    ops::{
        Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Div, Mul,
        MulAssign, Neg, Rem, Shl, Shr, Sub, SubAssign,
    },
    str::FromStr,
};

mod functions;
pub use functions::choose;
pub use functions::gcd;
pub use functions::lcm;

/// Represents a natural number {0, 1, 2, ...}
#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, CanonicalStructure)]
pub struct Natural(malachite_nz::natural::Natural);

impl Natural {
    pub(crate) fn from_malachite(value: malachite_nz::natural::Natural) -> Self {
        Self(value)
    }
    pub(crate) fn to_malachite(self) -> malachite_nz::natural::Natural {
        self.0
    }
    pub(crate) fn to_malachite_ref(&self) -> &malachite_nz::natural::Natural {
        &self.0
    }
}

impl ToStringSignature for NaturalCanonicalStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{}", elem)
    }
}

impl OrdSignature for NaturalCanonicalStructure {
    fn cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Natural::cmp(a, b)
    }
}

impl std::fmt::Display for Natural {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl From<u8> for Natural {
    fn from(value: u8) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<u16> for Natural {
    fn from(value: u16) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<u32> for Natural {
    fn from(value: u32) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<u64> for Natural {
    fn from(value: u64) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<u128> for Natural {
    fn from(value: u128) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<usize> for Natural {
    fn from(value: usize) -> Self {
        Self(malachite_nz::natural::Natural::from(value))
    }
}
impl From<&Natural> for Natural {
    fn from(value: &Natural) -> Self {
        value.clone()
    }
}

impl FromStr for Natural {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(malachite_nz::natural::Natural::from_str(s)?))
    }
}

impl Natural {
    pub const ZERO: Self = Self(malachite_nz::natural::Natural::ZERO);
    pub const ONE: Self = Self(malachite_nz::natural::Natural::ONE);
    pub const TWO: Self = Self(malachite_nz::natural::Natural::TWO);
}

impl AddAssign<Natural> for Natural {
    fn add_assign(&mut self, rhs: Natural) {
        self.0.add_assign(rhs.0)
    }
}
impl AddAssign<&Natural> for Natural {
    fn add_assign(&mut self, rhs: &Natural) {
        self.0.add_assign(&rhs.0)
    }
}

impl SubAssign<Natural> for Natural {
    fn sub_assign(&mut self, rhs: Natural) {
        self.0.sub_assign(rhs.0)
    }
}
impl SubAssign<&Natural> for Natural {
    fn sub_assign(&mut self, rhs: &Natural) {
        self.0.sub_assign(&rhs.0)
    }
}

impl MulAssign<Natural> for Natural {
    fn mul_assign(&mut self, rhs: Natural) {
        self.0.mul_assign(rhs.0)
    }
}
impl MulAssign<&Natural> for Natural {
    fn mul_assign(&mut self, rhs: &Natural) {
        self.0.mul_assign(&rhs.0)
    }
}

impl BitAndAssign<Natural> for Natural {
    fn bitand_assign(&mut self, rhs: Natural) {
        self.0.bitand_assign(rhs.0)
    }
}
impl BitAndAssign<&Natural> for Natural {
    fn bitand_assign(&mut self, rhs: &Natural) {
        self.0.bitand_assign(&rhs.0)
    }
}

impl BitOrAssign<Natural> for Natural {
    fn bitor_assign(&mut self, rhs: Natural) {
        self.0.bitor_assign(rhs.0)
    }
}
impl BitOrAssign<&Natural> for Natural {
    fn bitor_assign(&mut self, rhs: &Natural) {
        self.0.bitor_assign(&rhs.0)
    }
}

impl BitXorAssign<Natural> for Natural {
    fn bitxor_assign(&mut self, rhs: Natural) {
        self.0.bitxor_assign(rhs.0)
    }
}
impl BitXorAssign<&Natural> for Natural {
    fn bitxor_assign(&mut self, rhs: &Natural) {
        self.0.bitxor_assign(&rhs.0)
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
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<&Natural> for Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: &Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = self.0.div_mod(&other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<Natural> for &Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
        let (q, r) = (&self.0).div_mod(other.0);
        (Natural(q), Natural(r))
    }
}
impl DivMod<&Natural> for &Natural {
    type DivOutput = Natural;

    type ModOutput = Natural;

    fn div_mod(self, other: &Natural) -> (Self::DivOutput, Self::ModOutput) {
        use malachite_base::num::arithmetic::traits::DivMod;
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
    malachite_nz::natural::Natural: Shl<T, Output = malachite_nz::natural::Natural>,
{
    type Output = Natural;

    fn shl(self, rhs: T) -> Self::Output {
        Natural(self.0.shl(rhs))
    }
}
impl<T> Shl<T> for &Natural
where
    for<'a> &'a malachite_nz::natural::Natural: Shl<T, Output = malachite_nz::natural::Natural>,
{
    type Output = Natural;

    fn shl(self, rhs: T) -> Self::Output {
        Natural((&self.0).shl(rhs))
    }
}

impl<T> Shr<T> for Natural
where
    malachite_nz::natural::Natural: Shr<T, Output = malachite_nz::natural::Natural>,
{
    type Output = Natural;

    fn shr(self, rhs: T) -> Self::Output {
        Natural(self.0.shr(rhs))
    }
}
impl<T> Shr<T> for &Natural
where
    for<'a> &'a malachite_nz::natural::Natural: Shr<T, Output = malachite_nz::natural::Natural>,
{
    type Output = Natural;

    fn shr(self, rhs: T) -> Self::Output {
        Natural((&self.0).shr(rhs))
    }
}

impl AbsDiff<Natural> for Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: Natural) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural(self.0.abs_diff(rhs.0))
    }
}
impl AbsDiff<&Natural> for Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: &Natural) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural(self.0.abs_diff(&rhs.0))
    }
}
impl AbsDiff<Natural> for &Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: Natural) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural((&self.0).abs_diff(rhs.0))
    }
}
impl AbsDiff<&Natural> for &Natural {
    type Output = Natural;

    fn abs_diff(self, rhs: &Natural) -> Self::Output {
        use malachite_base::num::arithmetic::traits::AbsDiff;
        Natural((&self.0).abs_diff(&rhs.0))
    }
}

impl<Exponent: Borrow<Natural>, Modulus: Borrow<Natural>> ModPow<Exponent, Modulus> for Natural {
    type Output = Natural;
    fn mod_pow(self, exp: Exponent, m: Modulus) -> Self::Output {
        use malachite_base::num::arithmetic::traits::ModPow;
        Natural((self.0 % &m.borrow().0).mod_pow(&exp.borrow().0, &m.borrow().0))
    }
}

impl<Exponent: Borrow<Natural>, Modulus: Borrow<Natural>> ModPow<Exponent, Modulus> for &Natural {
    type Output = Natural;
    fn mod_pow(self, exp: Exponent, m: Modulus) -> Self::Output {
        use malachite_base::num::arithmetic::traits::ModPow;
        Natural((&self.0 % &m.borrow().0).mod_pow(&exp.borrow().0, &m.borrow().0))
    }
}

impl<Modulus: Borrow<Natural>> ModInv<Modulus> for Natural {
    type Output = Option<Natural>;
    fn mod_inv(self, m: Modulus) -> Self::Output {
        use malachite_base::num::arithmetic::traits::ModInverse;
        Some(Natural(
            (self.0 % &m.borrow().0).mod_inverse(&m.borrow().0)?,
        ))
    }
}

impl<Modulus: Borrow<Natural>> ModInv<Modulus> for &Natural {
    type Output = Option<Natural>;
    fn mod_inv(self, m: Modulus) -> Self::Output {
        use malachite_base::num::arithmetic::traits::ModInverse;
        Some(Natural(
            (&self.0 % &m.borrow().0).mod_inverse(&m.borrow().0)?,
        ))
    }
}

impl Natural {
    /// 2 raised to the power of `pow`.
    /// ```
    /// use algebraeon_nzq::Natural;
    /// assert_eq!(
    ///     Natural::from(32u32),
    ///     Natural::power_of_2(5)
    /// );
    /// ```
    pub fn power_of_2(pow: u64) -> Self {
        Self(malachite_nz::natural::Natural::power_of_2(pow))
    }

    /// An iterator over the bits in the binary expansion.
    /// ```
    /// use algebraeon_nzq::Natural;
    /// assert_eq!(
    ///     Natural::from(11u32).bits().collect::<Vec<_>>(),
    ///     vec![true, true, false, true],
    /// );
    /// ```
    pub fn bits<'a>(
        &'a self,
    ) -> impl Iterator<Item = bool> + ExactSizeIterator + DoubleEndedIterator + 'a {
        use malachite_base::num::logic::traits::BitIterable;
        self.0.bits()
    }

    /// Return the number of bits needed to store n i.e. ceil(log2(n)) for all non-zero n
    pub fn bitcount(&self) -> usize {
        self.bits().len()
    }
}

macro_rules! impl_try_into_unsigned {
    ($($t:ty),*) => {
        $(
            impl TryInto<$t> for Natural {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    (&self).try_into()
                }
            }
            impl TryInto<$t> for &Natural {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    let limbs = self.0.to_limbs_asc();
                    match limbs.len() {
                        0 => Ok(0),
                        1 => {
                            let n = limbs[0];
                            if Natural::from(n) > Natural::from(<$t>::MAX) {
                                Err(())
                            } else {
                                Ok(n as $t)
                            }
                        },
                        2 => {
                            if std::mem::size_of::<$t>() >= 16 {
                                let low = limbs[0] as u128;
                                let high = limbs[1] as u128;
                                let value = (high << 64) | low;
                                if value > <$t>::MAX as u128 {
                                    Err(())
                                } else {
                                    Ok(value as $t)
                                }
                            } else {
                                Err(())
                            }
                        },
                        _ => Err(()),
                    }
                }
            }
        )*
    };
}

macro_rules! impl_try_into_signed {
    ($($t:ty),*) => {
        $(
            impl TryInto<$t> for Natural {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    (&self).try_into()
                }
            }
            impl TryInto<$t> for &Natural {
                type Error = ();

                fn try_into(self) -> Result<$t, Self::Error> {
                    let limbs = self.0.to_limbs_asc();
                    match limbs.len() {
                        0 => Ok(0),
                        1 => {
                            let n = limbs[0] as i128;
                            if n > <$t>::MAX as i128 {
                                Err(())
                            } else {
                                Ok(n as $t)
                            }
                        },
                        2 => {
                            if std::mem::size_of::<$t>() >= 16 {
                                let low = limbs[0] as u128;
                                let high = limbs[1] as u128;
                                let value = (high << 64) | low;
                                if value > <$t>::MAX as u128 {
                                    Err(())
                                } else {
                                    Ok(value as $t)
                                }
                            } else {
                                Err(())
                            }
                        },
                        _ => Err(()),
                    }
                }
            }
        )*
    };
}

impl_try_into_unsigned!(u8, u16, u32, u64, u128, usize);
impl_try_into_signed!(i8, i16, i32, i64, i128, isize);

impl CountableSetSignature for NaturalCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
        use malachite_nz::natural::exhaustive::exhaustive_naturals;
        exhaustive_naturals()
            .into_iter()
            .map(|n| Natural::from_malachite(n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
