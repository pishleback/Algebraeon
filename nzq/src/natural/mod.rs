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

#[derive(Default, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

impl Natural {
    pub fn power_of_2(pow: u64) -> Self {
        Self(malachite_nz::natural::Natural::power_of_2(pow))
    }

    pub fn mod_pow(self, exp: impl Borrow<Self>, m: impl Borrow<Self>) -> Self {
        use malachite_base::num::arithmetic::traits::ModPow;
        Self(self.0.mod_pow(&exp.borrow().0, &m.borrow().0))
    }

    pub fn mod_pow_ref(&self, exp: impl Borrow<Self>, m: impl Borrow<Self>) -> Self {
        use malachite_base::num::arithmetic::traits::ModPow;
        Self((&self.0).mod_pow(&exp.borrow().0, &m.borrow().0))
    }

    pub fn mod_inv(self, modulus: &Natural) -> Result<Self, ()> {
        use malachite_base::num::arithmetic::traits::ExtendedGcd;
        if modulus.0 == malachite_nz::natural::Natural::ZERO {
            Err(())
        } else {
            let (g, a, _b) = self.0.extended_gcd(&modulus.0);
            if g == malachite_nz::natural::Natural::ONE {
                Ok(Integer::from_malachite(a) % modulus)
            } else {
                Err(())
            }
        }
    }

    pub fn mod_inv_ref(&self, modulus: &Natural) -> Result<Self, ()> {
        use malachite_base::num::arithmetic::traits::ExtendedGcd;
        if modulus.0 == malachite_nz::natural::Natural::ZERO {
            Err(())
        } else {
            let (g, a, _b) = (&self.0).extended_gcd(&modulus.0);
            if g == malachite_nz::natural::Natural::ONE {
                Ok(Integer::from_malachite(a) % modulus)
            } else {
                Err(())
            }
        }
    }

    pub fn bits<'a>(
        &'a self,
    ) -> impl Iterator<Item = bool> + ExactSizeIterator + DoubleEndedIterator + 'a {
        use malachite_base::num::logic::traits::BitIterable;
        self.0.bits()
    }
}

pub fn primes() -> impl Iterator<Item = usize> {
    use malachite_base::num::factorization::traits::Primes;
    usize::primes()
}

impl TryInto<usize> for Natural {
    type Error = ();

    fn try_into(self) -> Result<usize, Self::Error> {
        (&self).try_into()
    }
}
impl TryInto<usize> for &Natural {
    type Error = ();

    fn try_into(self) -> Result<usize, Self::Error> {
        let limbs = self.0.to_limbs_asc();
        if limbs.len() == 0 {
            Ok(0)
        } else if limbs.len() == 1 {
            let n = limbs[0];
            if Natural::from(n) > Natural::from(usize::MAX) {
                Err(())
            } else {
                Ok(n as usize)
            }
        } else {
            Err(())
        }
    }
}

impl MetaType for Natural {
    type Structure = CannonicalStructure<Natural>;

    fn structure() -> std::rc::Rc<Self::Structure> {
        CannonicalStructure::new().into()
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
}
