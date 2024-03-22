use std::fmt::Display;

use malachite_nz::natural::Natural;

use super::super::ring::*;

//the finite field of 4 elements
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum QuaternaryField {
    Zero,
    One,
    Alpha,
    Beta,
}

impl Display for QuaternaryField {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            QuaternaryField::Zero => write!(f, "0"),
            QuaternaryField::One => write!(f, "1"),
            QuaternaryField::Alpha => write!(f, "a"),
            QuaternaryField::Beta => write!(f, "b"),
        }
    }
}

impl ComRing for QuaternaryField {
    fn zero() -> Self {
        Self::Zero
    }

    fn one() -> Self {
        Self::One
    }

    fn neg_mut(&mut self) {}

    fn add_mut(&mut self, offset: &Self) {
        match (&self, offset) {
            (_, QuaternaryField::Zero) => {}
            (QuaternaryField::Zero, x) => *self = x.clone(),
            (QuaternaryField::One, QuaternaryField::One) => *self = Self::Zero,
            (QuaternaryField::One, QuaternaryField::Alpha) => *self = Self::Beta,
            (QuaternaryField::One, QuaternaryField::Beta) => *self = Self::Alpha,
            (QuaternaryField::Alpha, QuaternaryField::One) => *self = Self::Beta,
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => *self = Self::Zero,
            (QuaternaryField::Alpha, QuaternaryField::Beta) => *self = Self::One,
            (QuaternaryField::Beta, QuaternaryField::One) => *self = Self::Alpha,
            (QuaternaryField::Beta, QuaternaryField::Alpha) => *self = Self::One,
            (QuaternaryField::Beta, QuaternaryField::Beta) => *self = Self::Zero,
        }
    }

    fn mul_mut(&mut self, mul: &Self) {
        match (&self, mul) {
            (_, QuaternaryField::Zero) => *self = Self::Zero,
            (QuaternaryField::Zero, _) => *self = Self::Zero,
            (_, QuaternaryField::One) => {}
            (QuaternaryField::One, x) => *self = x.clone(),
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => *self = Self::Beta,
            (QuaternaryField::Alpha, QuaternaryField::Beta) => *self = Self::One,
            (QuaternaryField::Beta, QuaternaryField::Alpha) => *self = Self::One,
            (QuaternaryField::Beta, QuaternaryField::Beta) => *self = Self::Alpha,
        }
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
        match (&a, &b) {
            (_, QuaternaryField::Zero) => Err(RingDivisionError::DivideByZero),
            (_, QuaternaryField::One) => Ok(a),
            (QuaternaryField::Zero, _) => Ok(Self::Zero),
            (QuaternaryField::One, QuaternaryField::Alpha) => Ok(Self::Beta),
            (QuaternaryField::One, QuaternaryField::Beta) => Ok(Self::Alpha),
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => Ok(Self::One),
            (QuaternaryField::Alpha, QuaternaryField::Beta) => Ok(Self::Beta),
            (QuaternaryField::Beta, QuaternaryField::Alpha) => Ok(Self::Alpha),
            (QuaternaryField::Beta, QuaternaryField::Beta) => Ok(Self::One),
        }
    }
}

impl IntegralDomain for QuaternaryField {}

impl Field for QuaternaryField {}

impl FiniteUnits for QuaternaryField {
    fn all_units() -> Vec<Self> {
        vec![Self::One, Self::Alpha, Self::Beta]
    }
}

impl FiniteField for QuaternaryField {
    fn characteristic_and_power() -> (Natural, Natural) {
        (Natural::from(2u8), Natural::from(2u8))
    }
}

impl UniqueFactorizationDomain for super::super::polynomial::poly::Polynomial<QuaternaryField> {
    fn factor(&self) -> Option<Factored<Self>> {
        self.clone().factorize_by_berlekamps_algorithm()
    }
}

#[cfg(test)]
mod tests {
    use crate::{rings_old::numbers::small_fields::QuaternaryField, ComRing, RingDivisionError};

    #[test]
    fn test_neg() {
        assert_eq!(QuaternaryField::Zero.neg(), QuaternaryField::Zero);
        assert_eq!(QuaternaryField::One.neg(), QuaternaryField::One);
        assert_eq!(QuaternaryField::Alpha.neg(), QuaternaryField::Alpha);
        assert_eq!(QuaternaryField::Beta.neg(), QuaternaryField::Beta);
    }

    #[test]
    fn test_add() {
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Zero, QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Zero, QuaternaryField::One),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Zero, QuaternaryField::Alpha),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Zero, QuaternaryField::Beta),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::One, QuaternaryField::Zero),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::One, QuaternaryField::One),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::One, QuaternaryField::Alpha),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::One, QuaternaryField::Beta),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Alpha, QuaternaryField::Zero),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Alpha, QuaternaryField::One),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Alpha, QuaternaryField::Alpha),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Alpha, QuaternaryField::Beta),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Beta, QuaternaryField::Zero),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Beta, QuaternaryField::One),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Beta, QuaternaryField::Alpha),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(QuaternaryField::Beta, QuaternaryField::Beta),
            QuaternaryField::Zero
        );
    }

    #[test]
    fn test_mul() {
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Zero, QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Zero, QuaternaryField::One),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Zero, QuaternaryField::Alpha),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Zero, QuaternaryField::Beta),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::One, QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::One, QuaternaryField::One),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::One, QuaternaryField::Alpha),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::One, QuaternaryField::Beta),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Alpha, QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Alpha, QuaternaryField::One),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Alpha, QuaternaryField::Alpha),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Alpha, QuaternaryField::Beta),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Beta, QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Beta, QuaternaryField::One),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Beta, QuaternaryField::Alpha),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(QuaternaryField::Beta, QuaternaryField::Beta),
            QuaternaryField::Alpha
        );
    }

    #[test]
    fn test_div() {
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Zero, QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Zero, QuaternaryField::One),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Zero, QuaternaryField::Alpha),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Zero, QuaternaryField::Beta),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::One, QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::One, QuaternaryField::One),
            Ok(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::One, QuaternaryField::Alpha),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::One, QuaternaryField::Beta),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Alpha, QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Alpha, QuaternaryField::One),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Alpha, QuaternaryField::Alpha),
            Ok(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Alpha, QuaternaryField::Beta),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Beta, QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Beta, QuaternaryField::One),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Beta, QuaternaryField::Alpha),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(QuaternaryField::Beta, QuaternaryField::Beta),
            Ok(QuaternaryField::One)
        );
    }
}
