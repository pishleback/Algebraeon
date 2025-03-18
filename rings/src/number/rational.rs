use std::rc::Rc;
use std::{
    ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};

use super::super::polynomial::polynomial::*;
use super::super::structure::factorization::*;
use super::super::structure::structure::*;

use super::integer::*;
use super::natural::*;

use algebraeon_sets::structure::*;
use malachite_base::num::basic::traits::{One, OneHalf, Two, Zero};
use malachite_q::arithmetic::traits::SimplestRationalInInterval;
use rayon::collections::binary_heap::Iter;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Rational(malachite_q::Rational);

impl Rational {
    pub(crate) fn from_malachite(value: malachite_q::Rational) -> Self {
        Self(value)
    }
    pub(crate) fn to_malachite(self) -> malachite_q::Rational {
        self.0
    }
    pub(crate) fn to_malachite_ref(&self) -> &malachite_q::Rational {
        &self.0
    }
}

impl std::fmt::Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl From<u8> for Rational {
    fn from(value: u8) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<u16> for Rational {
    fn from(value: u16) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<u32> for Rational {
    fn from(value: u32) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<u64> for Rational {
    fn from(value: u64) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<u128> for Rational {
    fn from(value: u128) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<usize> for Rational {
    fn from(value: usize) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<i8> for Rational {
    fn from(value: i8) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<i16> for Rational {
    fn from(value: i16) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<i32> for Rational {
    fn from(value: i32) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<i64> for Rational {
    fn from(value: i64) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<i128> for Rational {
    fn from(value: i128) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<isize> for Rational {
    fn from(value: isize) -> Self {
        Self(malachite_q::Rational::from(value))
    }
}
impl From<Natural> for Rational {
    fn from(value: Natural) -> Self {
        Self(malachite_q::Rational::from(value.to_malachite()))
    }
}
impl From<&Natural> for Rational {
    fn from(value: &Natural) -> Self {
        Self(malachite_q::Rational::from(value.to_malachite_ref()))
    }
}
impl From<Integer> for Rational {
    fn from(value: Integer) -> Self {
        Self(malachite_q::Rational::from(value.to_malachite()))
    }
}
impl From<&Integer> for Rational {
    fn from(value: &Integer) -> Self {
        Self(malachite_q::Rational::from(value.to_malachite_ref()))
    }
}

impl FromStr for Rational {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(malachite_q::Rational::from_str(s)?))
    }
}

impl Rational {
    pub const ZERO: Self = Self(malachite_q::Rational::ZERO);
    pub const ONE: Self = Self(malachite_q::Rational::ONE);
    pub const TWO: Self = Self(malachite_q::Rational::TWO);
    pub const ONE_HALF: Self = Self(malachite_q::Rational::ONE_HALF);
}

impl AddAssign<Rational> for Rational {
    fn add_assign(&mut self, rhs: Rational) {
        self.0.add_assign(rhs.0)
    }
}
impl AddAssign<&Rational> for Rational {
    fn add_assign(&mut self, rhs: &Rational) {
        self.0.add_assign(&rhs.0)
    }
}

impl SubAssign<Rational> for Rational {
    fn sub_assign(&mut self, rhs: Rational) {
        self.0.sub_assign(rhs.0)
    }
}
impl SubAssign<&Rational> for Rational {
    fn sub_assign(&mut self, rhs: &Rational) {
        self.0.sub_assign(&rhs.0)
    }
}

impl MulAssign<Rational> for Rational {
    fn mul_assign(&mut self, rhs: Rational) {
        self.0.mul_assign(rhs.0)
    }
}
impl MulAssign<&Rational> for Rational {
    fn mul_assign(&mut self, rhs: &Rational) {
        self.0.mul_assign(&rhs.0)
    }
}

impl Neg for Rational {
    type Output = Rational;

    fn neg(self) -> Self::Output {
        Rational(self.0.neg())
    }
}
impl Neg for &Rational {
    type Output = Rational;

    fn neg(self) -> Self::Output {
        Rational((&self.0).neg())
    }
}

impl Add<Rational> for Rational {
    type Output = Rational;

    fn add(self, rhs: Rational) -> Self::Output {
        Rational(self.0.add(rhs.0))
    }
}
impl Add<&Rational> for Rational {
    type Output = Rational;

    fn add(self, rhs: &Rational) -> Self::Output {
        Rational(self.0.add(&rhs.0))
    }
}
impl Add<Rational> for &Rational {
    type Output = Rational;

    fn add(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).add(rhs.0))
    }
}
impl Add<&Rational> for &Rational {
    type Output = Rational;

    fn add(self, rhs: &Rational) -> Self::Output {
        Rational((&self.0).add(&rhs.0))
    }
}

impl Sub<Rational> for Rational {
    type Output = Rational;

    fn sub(self, rhs: Rational) -> Self::Output {
        Rational(self.0.sub(rhs.0))
    }
}
impl Sub<&Rational> for Rational {
    type Output = Rational;

    fn sub(self, rhs: &Rational) -> Self::Output {
        Rational(self.0.sub(&rhs.0))
    }
}
impl Sub<Rational> for &Rational {
    type Output = Rational;

    fn sub(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).sub(rhs.0))
    }
}
impl Sub<&Rational> for &Rational {
    type Output = Rational;

    fn sub(self, rhs: &Rational) -> Self::Output {
        Rational((&self.0).sub(&rhs.0))
    }
}

impl Mul<Rational> for Rational {
    type Output = Rational;

    fn mul(self, rhs: Rational) -> Self::Output {
        Rational(self.0.mul(rhs.0))
    }
}
impl Mul<&Rational> for Rational {
    type Output = Rational;

    fn mul(self, rhs: &Rational) -> Self::Output {
        Rational(self.0.mul(&rhs.0))
    }
}
impl Mul<Rational> for &Rational {
    type Output = Rational;

    fn mul(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).mul(rhs.0))
    }
}
impl Mul<&Rational> for &Rational {
    type Output = Rational;

    fn mul(self, rhs: &Rational) -> Self::Output {
        Rational((&self.0).mul(&rhs.0))
    }
}

impl Div<Rational> for Rational {
    type Output = Rational;

    fn div(self, rhs: Rational) -> Self::Output {
        Rational(self.0.div(rhs.0))
    }
}
impl Div<&Rational> for Rational {
    type Output = Rational;

    fn div(self, rhs: &Rational) -> Self::Output {
        Rational(self.0.div(&rhs.0))
    }
}
impl Div<Rational> for &Rational {
    type Output = Rational;

    fn div(self, rhs: Rational) -> Self::Output {
        Rational((&self.0).div(rhs.0))
    }
}
impl Div<&Rational> for &Rational {
    type Output = Rational;

    fn div(self, rhs: &Rational) -> Self::Output {
        Rational((&self.0).div(&rhs.0))
    }
}

impl Rational {
    pub fn numerator(&self) -> Integer {
        //malachite returns a natural for the numerator for some
        if self >= &Rational::ZERO {
            Integer::from(Natural::from_malachite(self.0.numerator_ref().clone()))
        } else {
            -Natural::from_malachite(self.0.numerator_ref().clone())
        }
    }

    pub fn denominator(&self) -> Natural {
        Natural::from_malachite(self.0.denominator_ref().clone())
    }

    pub fn from_integers(n: impl Into<Integer>, d: impl Into<Integer>) -> Self {
        Self(malachite_q::Rational::from_integers(
            n.into().to_malachite(),
            d.into().to_malachite(),
        ))
    }

    pub fn abs(self) -> Self {
        use malachite_base::num::arithmetic::traits::Abs;
        Self(self.0.abs())
    }

    pub fn abs_ref(&self) -> Self {
        use malachite_base::num::arithmetic::traits::Abs;
        Self((&self.0).abs())
    }

    pub fn into_abs_numerator_and_denominator(self) -> (Natural, Natural) {
        let (n, d) = self.0.into_numerator_and_denominator();
        (Natural::from_malachite(n), Natural::from_malachite(d))
    }

    pub fn simplest_rational_in_closed_interval(a: &Rational, b: &Rational) -> Self {
        Self(malachite_q::Rational::simplest_rational_in_closed_interval(
            &a.0, &b.0,
        ))
    }

    pub fn simplest_rational_in_open_interval(a: &Rational, b: &Rational) -> Self {
        Self(malachite_q::Rational::simplest_rational_in_open_interval(
            &a.0, &b.0,
        ))
    }

    pub fn exhaustive_rationals() -> impl Iterator<Item = Rational> {
        malachite_q::exhaustive::exhaustive_rationals().map(|v| Rational(v))
    }
}

impl MetaType for Rational {
    type Structure = CannonicalStructure<Rational>;

    fn structure() -> std::rc::Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl SemiRingStructure for CannonicalStructure<Rational> {
    fn zero(&self) -> Self::Set {
        Rational::ZERO
    }

    fn one(&self) -> Self::Set {
        Rational::ONE
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl RingStructure for CannonicalStructure<Rational> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }
}

impl IntegralDomainStructure for CannonicalStructure<Rational> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if b == &Rational::ZERO {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}

impl OrderedRingStructure for CannonicalStructure<Rational> {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FieldStructure for CannonicalStructure<Rational> {}

impl FieldOfFractionsStructure for CannonicalStructure<Rational> {
    type RS = CannonicalStructure<Integer>;

    fn base_ring_structure(&self) -> Rc<Self::RS> {
        Integer::structure()
    }

    fn from_base_ring(&self, elem: <Self::RS as Structure>::Set) -> Self::Set {
        Rational::from(elem)
    }

    fn numerator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set {
        elem.numerator()
    }

    fn denominator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set {
        Integer::from(elem.denominator())
    }
}

impl RealRoundingStructure for CannonicalStructure<Rational> {
    fn floor(&self, x: &Self::Set) -> Integer {
        Integer::from_malachite(
            <malachite_q::Rational as malachite_base::num::arithmetic::traits::Floor>::floor(
                x.0.clone(),
            ),
        )
    }
    fn ceil(&self, x: &Self::Set) -> Integer {
        -self.floor(&-x)
    }
    fn round(&self, x: &Self::Set) -> Integer {
        self.floor(&(x + Rational::ONE_HALF))
    }
}

impl RealFromFloatStructure for CannonicalStructure<Rational> {
    fn from_f64_approx(&self, x: f64) -> Self::Set {
        Rational(malachite_q::Rational::try_from_float_simplest(x).unwrap())
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Rational>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
        self.factorize_by_factorize_primitive_part(p)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_rational_floor_ceil_round() {
        let rat = |s: &'static str| Rational::from_str(s).unwrap();

        assert_eq!(rat("-2").floor(), Integer::from(-2));
        assert_eq!(rat("-7/4").floor(), Integer::from(-2));
        assert_eq!(rat("-3/2").floor(), Integer::from(-2));
        assert_eq!(rat("-5/4").floor(), Integer::from(-2));
        assert_eq!(rat("-1").floor(), Integer::from(-1));
        assert_eq!(rat("-3/4").floor(), Integer::from(-1));
        assert_eq!(rat("-1/2").floor(), Integer::from(-1));
        assert_eq!(rat("-1/4").floor(), Integer::from(-1));
        assert_eq!(rat("0").floor(), Integer::from(0));
        assert_eq!(rat("1/4").floor(), Integer::from(0));
        assert_eq!(rat("1/2").floor(), Integer::from(0));
        assert_eq!(rat("3/4").floor(), Integer::from(0));
        assert_eq!(rat("1").floor(), Integer::from(1));
        assert_eq!(rat("5/4").floor(), Integer::from(1));
        assert_eq!(rat("3/2").floor(), Integer::from(1));
        assert_eq!(rat("7/4").floor(), Integer::from(1));
        assert_eq!(rat("2").floor(), Integer::from(2));

        assert_eq!(rat("-2").ceil(), Integer::from(-2));
        assert_eq!(rat("-7/4").ceil(), Integer::from(-1));
        assert_eq!(rat("-3/2").ceil(), Integer::from(-1));
        assert_eq!(rat("-5/4").ceil(), Integer::from(-1));
        assert_eq!(rat("-1").ceil(), Integer::from(-1));
        assert_eq!(rat("-3/4").ceil(), Integer::from(0));
        assert_eq!(rat("-1/2").ceil(), Integer::from(0));
        assert_eq!(rat("-1/4").ceil(), Integer::from(0));
        assert_eq!(rat("0").ceil(), Integer::from(0));
        assert_eq!(rat("1/4").ceil(), Integer::from(1));
        assert_eq!(rat("1/2").ceil(), Integer::from(1));
        assert_eq!(rat("3/4").ceil(), Integer::from(1));
        assert_eq!(rat("1").ceil(), Integer::from(1));
        assert_eq!(rat("5/4").ceil(), Integer::from(2));
        assert_eq!(rat("3/2").ceil(), Integer::from(2));
        assert_eq!(rat("7/4").ceil(), Integer::from(2));
        assert_eq!(rat("2").ceil(), Integer::from(2));

        assert_eq!(rat("-2").round(), Integer::from(-2));
        assert_eq!(rat("-7/4").round(), Integer::from(-2));
        assert!(vec![Integer::from(-2), Integer::from(-1)].contains(&rat("-3/2").round()));
        assert_eq!(rat("-5/4").round(), Integer::from(-1));
        assert_eq!(rat("-1").round(), Integer::from(-1));
        assert_eq!(rat("-3/4").round(), Integer::from(-1));
        assert!(vec![Integer::from(-1), Integer::from(0)].contains(&rat("-1/2").round()));
        assert_eq!(rat("-1/4").round(), Integer::from(0));
        assert_eq!(rat("0").round(), Integer::from(0));
        assert_eq!(rat("1/4").round(), Integer::from(0));
        assert!(vec![Integer::from(0), Integer::from(1)].contains(&rat("1/2").round()));
        assert_eq!(rat("3/4").round(), Integer::from(1));
        assert_eq!(rat("1").round(), Integer::from(1));
        assert_eq!(rat("5/4").round(), Integer::from(1));
        assert!(vec![Integer::from(1), Integer::from(2)].contains(&rat("3/2").round()));
        assert_eq!(rat("7/4").round(), Integer::from(2));
        assert_eq!(rat("2").round(), Integer::from(2));
    }
}
