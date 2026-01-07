use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::fmt::Display;

//the finite field of 4 elements
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuaternaryField {
    Zero,
    One,
    Alpha,
    Beta,
}

impl ToStringSignature for QuaternaryFieldCanonicalStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{}", elem)
    }
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

impl SetWithZeroSignature for QuaternaryFieldCanonicalStructure {
    fn zero(&self) -> Self::Set {
        QuaternaryField::Zero
    }
}

impl AdditiveMonoidSignature for QuaternaryFieldCanonicalStructure {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        #[allow(clippy::match_same_arms)]
        match (a, b) {
            (_, QuaternaryField::Zero) => *a,
            (QuaternaryField::Zero, _) => *b,
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

    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }

    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl MultiplicativeMonoidSignature for QuaternaryFieldCanonicalStructure {
    fn one(&self) -> Self::Set {
        QuaternaryField::One
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        #[allow(clippy::match_same_arms)]
        match (a, b) {
            (_, QuaternaryField::Zero) => QuaternaryField::Zero,
            (QuaternaryField::Zero, _) => QuaternaryField::Zero,
            (_, QuaternaryField::One) => *a,
            (QuaternaryField::One, _) => *b,
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => QuaternaryField::Beta,
            (QuaternaryField::Alpha, QuaternaryField::Beta) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::Alpha) => QuaternaryField::One,
            (QuaternaryField::Beta, QuaternaryField::Beta) => QuaternaryField::Alpha,
        }
    }
}

impl SemiRingSignature for QuaternaryFieldCanonicalStructure {}

impl RingSignature for QuaternaryFieldCanonicalStructure {
    fn is_reduced(&self) -> Result<bool, String> {
        Ok(true)
    }
}

impl CharacteristicSignature for QuaternaryFieldCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::TWO
    }
}

impl AdditiveGroupSignature for QuaternaryFieldCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        *a
    }
}

impl MultiplicativeMonoidUnitsSignature for QuaternaryFieldCanonicalStructure {
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set> {
        self.try_div(&self.one(), a)
    }
}

impl IntegralDomainSignature for QuaternaryFieldCanonicalStructure {
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        #[allow(clippy::match_same_arms)]
        match (&a, &b) {
            (_, QuaternaryField::Zero) => None,
            (_, QuaternaryField::One) => Some(*a),
            (QuaternaryField::Zero, _) => Some(QuaternaryField::Zero),
            (QuaternaryField::One, QuaternaryField::Alpha) => Some(QuaternaryField::Beta),
            (QuaternaryField::One, QuaternaryField::Beta) => Some(QuaternaryField::Alpha),
            (QuaternaryField::Alpha, QuaternaryField::Alpha) => Some(QuaternaryField::One),
            (QuaternaryField::Alpha, QuaternaryField::Beta) => Some(QuaternaryField::Beta),
            (QuaternaryField::Beta, QuaternaryField::Alpha) => Some(QuaternaryField::Alpha),
            (QuaternaryField::Beta, QuaternaryField::Beta) => Some(QuaternaryField::One),
        }
    }
}

impl FieldSignature for QuaternaryFieldCanonicalStructure {}

impl<B: BorrowedStructure<QuaternaryFieldCanonicalStructure>> CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<QuaternaryFieldCanonicalStructure, B>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        self.list_all_elements().into_iter()
    }
}

impl<B: BorrowedStructure<QuaternaryFieldCanonicalStructure>> FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<QuaternaryFieldCanonicalStructure, B>
{
    fn list_all_elements(&self) -> Vec<Self::Set> {
        vec![
            QuaternaryField::One,
            QuaternaryField::Alpha,
            QuaternaryField::Beta,
        ]
    }
}

impl CountableSetSignature for QuaternaryFieldCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        self.all_units_and_zero().into_iter()
    }
}

impl FiniteSetSignature for QuaternaryFieldCanonicalStructure {}

impl FiniteFieldSignature for QuaternaryFieldCanonicalStructure {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (Natural::from(2u8), Natural::from(2u8))
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
            QuaternaryField::try_div(&QuaternaryField::Zero, &QuaternaryField::Zero),
            None
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Zero, &QuaternaryField::One),
            Some(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Zero, &QuaternaryField::Alpha),
            Some(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Zero, &QuaternaryField::Beta),
            Some(QuaternaryField::Zero)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::One, &QuaternaryField::Zero),
            None
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::One, &QuaternaryField::One),
            Some(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::One, &QuaternaryField::Alpha),
            Some(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::One, &QuaternaryField::Beta),
            Some(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Alpha, &QuaternaryField::Zero),
            None
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Alpha, &QuaternaryField::One),
            Some(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Alpha, &QuaternaryField::Alpha),
            Some(QuaternaryField::One)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Alpha, &QuaternaryField::Beta),
            Some(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Beta, &QuaternaryField::Zero),
            None
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Beta, &QuaternaryField::One),
            Some(QuaternaryField::Beta)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Beta, &QuaternaryField::Alpha),
            Some(QuaternaryField::Alpha)
        );
        assert_eq!(
            QuaternaryField::try_div(&QuaternaryField::Beta, &QuaternaryField::Beta),
            Some(QuaternaryField::One)
        );
    }
}
