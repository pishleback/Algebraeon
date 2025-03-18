use std::fmt::Display;
use std::rc::Rc;

use crate::polynomial::polynomial::*;
use crate::structure::structure::*;
use algebraeon_sets::structure::*;

use algebraeon_nzq::natural::*;

//the finite field of 4 elements
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

impl MetaType for QuaternaryField {
    type Structure = CannonicalStructure<QuaternaryField>;

    fn structure() -> Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl SemiRingStructure for CannonicalStructure<QuaternaryField> {
    fn zero(&self) -> Self::Set {
        QuaternaryField::Zero
    }

    fn one(&self) -> Self::Set {
        QuaternaryField::One
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        match (a, b) {
            (_, QuaternaryField::Zero) => a.clone(),
            (QuaternaryField::Zero, _) => b.clone(),
            (QuaternaryField::One, QuaternaryField::One) => QuaternaryField::Zero,
            (QuaternaryField::One, QuaternaryField::Alpha) => QuaternaryField::Beta,
            (QuaternaryField::One, QuaternaryField::Beta) => QuaternaryField::Alpha,
            (QuaternaryField::Alpha, QuaternaryField::One) => QuaternaryField::Beta,
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => QuaternaryField::Zero,
            (QuaternaryField::Alpha, QuaternaryField::Beta) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::One) => QuaternaryField::Alpha,
            (QuaternaryField::Beta, QuaternaryField::Alpha) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::Beta) => QuaternaryField::Zero,
        }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        match (a, b) {
            (_, QuaternaryField::Zero) => QuaternaryField::Zero,
            (QuaternaryField::Zero, _) => QuaternaryField::Zero,
            (_, QuaternaryField::One) => a.clone(),
            (QuaternaryField::One, _) => b.clone(),
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => QuaternaryField::Beta,
            (QuaternaryField::Alpha, QuaternaryField::Beta) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::Alpha) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::Beta) => QuaternaryField::Alpha,
        }
    }
}

impl RingStructure for CannonicalStructure<QuaternaryField> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        a.clone()
    }
}

impl IntegralDomainStructure for CannonicalStructure<QuaternaryField> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match (&a, &b) {
            (_, QuaternaryField::Zero) => Err(RingDivisionError::DivideByZero),
            (_, QuaternaryField::One) => Ok(a.clone()),
            (QuaternaryField::Zero, _) => Ok(QuaternaryField::Zero),
            (QuaternaryField::One, QuaternaryField::Alpha) => Ok(QuaternaryField::Beta),
            (QuaternaryField::One, QuaternaryField::Beta) => Ok(QuaternaryField::Alpha),
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => Ok(QuaternaryField::One),
            (QuaternaryField::Alpha, QuaternaryField::Beta) => Ok(QuaternaryField::Beta),
            (QuaternaryField::Beta, QuaternaryField::Alpha) => Ok(QuaternaryField::Alpha),
            (QuaternaryField::Beta, QuaternaryField::Beta) => Ok(QuaternaryField::One),
        }
    }
}

impl FieldStructure for CannonicalStructure<QuaternaryField> {}

impl FiniteUnitsStructure for CannonicalStructure<QuaternaryField> {
    fn all_units(&self) -> Vec<Self::Set> {
        vec![
            QuaternaryField::One,
            QuaternaryField::Alpha,
            QuaternaryField::Beta,
        ]
    }
}

impl FiniteFieldStructure for CannonicalStructure<QuaternaryField> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (Natural::from(2u8), Natural::from(2u8))
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<QuaternaryField>> {
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::factorization::Factored<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
            QuaternaryField::add(&QuaternaryField::Zero, &QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Zero, &QuaternaryField::One),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Zero, &QuaternaryField::Alpha),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Zero, &QuaternaryField::Beta),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::One, &QuaternaryField::Zero),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::One, &QuaternaryField::One),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::One, &QuaternaryField::Alpha),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::One, &QuaternaryField::Beta),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Alpha, &QuaternaryField::Zero),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Alpha, &QuaternaryField::One),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Alpha, &QuaternaryField::Alpha),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Alpha, &QuaternaryField::Beta),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Beta, &QuaternaryField::Zero),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Beta, &QuaternaryField::One),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Beta, &QuaternaryField::Alpha),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::add(&QuaternaryField::Beta, &QuaternaryField::Beta),
            QuaternaryField::Zero
        );
    }

    #[test]
    fn test_mul() {
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Zero, &QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Zero, &QuaternaryField::One),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Zero, &QuaternaryField::Alpha),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Zero, &QuaternaryField::Beta),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::One, &QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::One, &QuaternaryField::One),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::One, &QuaternaryField::Alpha),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::One, &QuaternaryField::Beta),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Alpha, &QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Alpha, &QuaternaryField::One),
            QuaternaryField::Alpha
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Alpha, &QuaternaryField::Alpha),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Alpha, &QuaternaryField::Beta),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Beta, &QuaternaryField::Zero),
            QuaternaryField::Zero
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Beta, &QuaternaryField::One),
            QuaternaryField::Beta
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Beta, &QuaternaryField::Alpha),
            QuaternaryField::One
        );
        assert_eq!(
            QuaternaryField::mul(&QuaternaryField::Beta, &QuaternaryField::Beta),
            QuaternaryField::Alpha
        );
    }

    #[test]
    fn test_div() {
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Zero, &QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Zero, &QuaternaryField::One),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Zero, &QuaternaryField::Alpha),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Zero, &QuaternaryField::Beta),
            Ok(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::One, &QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::One, &QuaternaryField::One),
            Ok(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::One, &QuaternaryField::Alpha),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::One, &QuaternaryField::Beta),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Alpha, &QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Alpha, &QuaternaryField::One),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Alpha, &QuaternaryField::Alpha),
            Ok(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Alpha, &QuaternaryField::Beta),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Beta, &QuaternaryField::Zero),
            Err(RingDivisionError::DivideByZero)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Beta, &QuaternaryField::One),
            Ok(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Beta, &QuaternaryField::Alpha),
            Ok(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::div(&QuaternaryField::Beta, &QuaternaryField::Beta),
            Ok(QuaternaryField::One)
        );
    }
}
