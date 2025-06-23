use crate::algebraic_number_field::structure::{
    AlgebraicIntegerRingInAlgebraicNumberField, AlgebraicNumberFieldSignature,
};
use crate::polynomial::{Polynomial, PolynomialStructure, factorize_by_factorize_primitive_part};
use crate::structure::*;
use algebraeon_nzq::traits::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Cow;

impl AdditiveMonoidSignature for RationalCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Rational::ZERO
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
}

impl AdditiveGroupSignature for RationalCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a - b
    }
}

impl SemiRingSignature for RationalCanonicalStructure {
    fn one(&self) -> Self::Set {
        Rational::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl RingSignature for RationalCanonicalStructure {}

impl CharacteristicSignature for RationalCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl SemiRingUnitsSignature for RationalCanonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }
}

impl IntegralDomainSignature for RationalCanonicalStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if b == &Rational::ZERO {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}

impl OrderedRingSignature for RationalCanonicalStructure {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FieldSignature for RationalCanonicalStructure {}

impl CharZeroRingSignature for RationalCanonicalStructure {
    fn try_to_int(&self, x: &Rational) -> Option<Integer> {
        let (n, d) = x.numerator_and_denominator();
        debug_assert_ne!(&d, &Natural::ZERO);
        if d == Natural::ONE { Some(n) } else { None }
    }
}

impl CharZeroFieldSignature for RationalCanonicalStructure {
    fn try_to_rat(&self, x: &Rational) -> Option<Rational> {
        Some(x.clone())
    }
}

impl<'h, B: BorrowedStructure<RationalCanonicalStructure>>
    FinitelyGeneratedModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        RationalCanonicalStructure,
        PrincipalRationalSubfieldInclusion<RationalCanonicalStructure, B>,
    >
{
}

impl<'h, B: BorrowedStructure<RationalCanonicalStructure>>
    FreeModuleSignature<SingletonSetStructure, RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        RationalCanonicalStructure,
        PrincipalRationalSubfieldInclusion<RationalCanonicalStructure, B>,
    >
{
    fn basis_set(&self) -> impl std::borrow::Borrow<SingletonSetStructure> {
        SingletonSetStructure::default()
    }

    fn to_component<'a>(&self, _: &(), v: &'a Rational) -> Cow<'a, Rational> {
        Cow::Borrowed(v)
    }

    fn from_component(&self, _: &(), r: &Rational) -> Rational {
        r.clone()
    }
}

impl<'h, B: BorrowedStructure<RationalCanonicalStructure>>
    FinitelyFreeModuleSignature<SingletonSetStructure, RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        RationalCanonicalStructure,
        PrincipalRationalSubfieldInclusion<RationalCanonicalStructure, B>,
    >
{
}

impl ComplexSubsetSignature for RationalCanonicalStructure {
    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
        (self.as_f64(z), 0.0)
    }
}

impl RealSubsetSignature for RationalCanonicalStructure {}

impl RealToFloatSignature for RationalCanonicalStructure {
    fn as_f64(&self, x: &Rational) -> f64 {
        let fof = PrincipalSubringInclusion::new(self.clone());
        RealToFloatSignature::as_f64(&Integer::structure(), &fof.numerator(x))
            / RealToFloatSignature::as_f64(&Integer::structure(), &fof.denominator(x))
    }
}

impl<B: BorrowedStructure<RationalCanonicalStructure>>
    FieldOfFractionsInclusion<IntegerCanonicalStructure, RationalCanonicalStructure>
    for PrincipalSubringInclusion<RationalCanonicalStructure, B>
{
    fn numerator_and_denominator(&self, a: &Rational) -> (Integer, Integer) {
        (a.numerator(), a.denominator().into())
    }
}

impl RealRoundingSignature for RationalCanonicalStructure {
    fn floor(&self, x: &Self::Set) -> Integer {
        Floor::floor(x)
    }
    fn ceil(&self, x: &Self::Set) -> Integer {
        Ceil::ceil(x)
    }
    fn round(&self, x: &Self::Set) -> Integer {
        self.floor(&(x + Rational::ONE_HALF))
    }
}

impl RealFromFloatSignature for RationalCanonicalStructure {
    fn from_f64_approx(&self, x: f64) -> Self::Set {
        Rational::try_from_float_simplest(x).unwrap()
    }
}

impl AlgebraicIntegerRingInAlgebraicNumberField<RationalCanonicalStructure>
    for PrincipalSubringInclusion<RationalCanonicalStructure, RationalCanonicalStructure>
{
    fn discriminant(&self) -> Integer {
        Integer::ONE
    }
}

impl AlgebraicNumberFieldSignature for RationalCanonicalStructure {
    type Basis = SingletonSetStructure;
    type RingOfIntegers = IntegerCanonicalStructure;
    type RingOfIntegersInclusion = PrincipalSubringInclusion<Self, Self>;
    type RationalInclusion<B: BorrowedStructure<Self>> =
        PrincipalRationalSubfieldInclusion<Self, B>;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self> {
        self.rational_extension()
    }
    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        self.into_rational_extension()
    }

    fn into_ring_of_integers_extension(self) -> Self::RingOfIntegersInclusion {
        self.into_principal_subring_inclusion()
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool {
        self.try_to_int(a).is_some()
    }
}

impl<B: BorrowedStructure<RationalCanonicalStructure>> FactorableSignature
    for PolynomialStructure<RationalCanonicalStructure, B>
{
    fn factor(&self, p: &Self::Set) -> Option<FactoredRingElement<Polynomial<Rational>>> {
        factorize_by_factorize_primitive_part(
            &PrincipalSubringInclusion::new(self.coeff_ring().clone()),
            self,
            p,
        )
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
        assert!([Integer::from(-2), Integer::from(-1)].contains(&rat("-3/2").round()));
        assert_eq!(rat("-5/4").round(), Integer::from(-1));
        assert_eq!(rat("-1").round(), Integer::from(-1));
        assert_eq!(rat("-3/4").round(), Integer::from(-1));
        assert!([Integer::from(-1), Integer::from(0)].contains(&rat("-1/2").round()));
        assert_eq!(rat("-1/4").round(), Integer::from(0));
        assert_eq!(rat("0").round(), Integer::from(0));
        assert_eq!(rat("1/4").round(), Integer::from(0));
        assert!([Integer::from(0), Integer::from(1)].contains(&rat("1/2").round()));
        assert_eq!(rat("3/4").round(), Integer::from(1));
        assert_eq!(rat("1").round(), Integer::from(1));
        assert_eq!(rat("5/4").round(), Integer::from(1));
        assert!([Integer::from(1), Integer::from(2)].contains(&rat("3/2").round()));
        assert_eq!(rat("7/4").round(), Integer::from(2));
        assert_eq!(rat("2").round(), Integer::from(2));
    }
}
