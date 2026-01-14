use crate::algebraic_number_field::{AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature};
use crate::polynomial::{PolynomialStructure, factorize_by_factorize_primitive_part};
use crate::structure::*;
use algebraeon_nzq::traits::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use static_assertions::const_assert;
use std::borrow::Cow;

impl RinglikeSpecializationSignature for RationalCanonicalStructure {
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Set = Self::Set> + RingSignature> {
        Some(self.clone())
    }

    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Set = Self::Set> + CharZeroRingSignature> {
        Some(self.clone())
    }
}

impl ZeroSignature for RationalCanonicalStructure {
    fn zero(&self) -> Self::Set {
        Rational::ZERO
    }
}

impl AdditionSignature for RationalCanonicalStructure {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
}

impl CancellativeAdditionSignature for RationalCanonicalStructure {
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl TryNegateSignature for RationalCanonicalStructure {
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }
}

impl AdditiveMonoidSignature for RationalCanonicalStructure {}

impl AdditiveGroupSignature for RationalCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a - b
    }
}

impl MultiplicativeMonoidSignature for RationalCanonicalStructure {
    fn one(&self) -> Self::Set {
        Rational::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl SemiRingSignature for RationalCanonicalStructure {}

impl RingSignature for RationalCanonicalStructure {
    fn is_reduced(&self) -> Result<bool, String> {
        Ok(true)
    }
}

impl CharacteristicSignature for RationalCanonicalStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl MultiplicativeMonoidUnitsSignature for RationalCanonicalStructure {
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set> {
        self.try_div(&self.one(), a)
    }
}

impl MultiplicativeIntegralMonoidSignature for RationalCanonicalStructure {
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        if b == &Rational::ZERO {
            None
        } else {
            Some(a / b)
        }
    }
}

impl IntegralDomainSignature for RationalCanonicalStructure {}

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
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        RationalCanonicalStructure,
        PrincipalRationalMap<RationalCanonicalStructure, B>,
    >
{
    type Basis = SingletonSetStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        Self::Basis::default()
    }

    fn to_component<'a>(&self, _: &(), v: &'a Rational) -> Cow<'a, Rational> {
        Cow::Borrowed(v)
    }

    fn from_component(&self, _: &(), r: &Rational) -> Rational {
        r.clone()
    }
}

impl ComplexSubsetSignature for RationalCanonicalStructure {
    fn as_f32_real_and_imaginary_parts(&self, z: &Self::Set) -> (f32, f32) {
        (self.as_f32(z), 0.0)
    }

    fn as_f64_real_and_imaginary_parts(&self, z: &Self::Set) -> (f64, f64) {
        (self.as_f64(z), 0.0)
    }
}

impl RealSubsetSignature for RationalCanonicalStructure {
    fn as_f64(&self, x: &Rational) -> f64 {
        x.into()
    }

    fn as_f32(&self, x: &Self::Set) -> f32 {
        x.into()
    }
}

impl<B: BorrowedStructure<RationalCanonicalStructure>>
    FieldOfFractionsInclusion<IntegerCanonicalStructure, RationalCanonicalStructure>
    for PrincipalIntegerMap<RationalCanonicalStructure, B>
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

impl AlgebraicNumberFieldSignature for RationalCanonicalStructure {
    type Basis = SingletonSetStructure;
    type RationalInclusion<B: BorrowedStructure<Self>> = PrincipalRationalMap<Self, B>;

    fn inbound_finite_dimensional_rational_extension(&self) -> Self::RationalInclusion<&Self> {
        self.inbound_principal_rational_map()
    }
    fn into_inbound_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        self.into_inbound_principal_rational_map()
    }

    fn generator(&self) -> Rational {
        Rational::ONE
    }

    fn discriminant(&self) -> Integer {
        Integer::ONE
    }

    fn integral_basis(&self) -> Vec<Self::Set> {
        vec![Rational::ONE]
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool {
        self.try_to_int(a).is_some()
    }
}

const_assert!(
    impls::impls!(IntegerCanonicalStructure : AlgebraicIntegerRingSignature<RationalCanonicalStructure>)
);

impl<B: BorrowedStructure<RationalCanonicalStructure>> FactoringMonoidSignature
    for PolynomialStructure<RationalCanonicalStructure, B>
{
    fn factor_unchecked(
        &self,
        p: &Self::Set,
    ) -> Factored<Self::Set, <Self::FactoredExponent as SetSignature>::Set> {
        factorize_by_factorize_primitive_part(
            &PrincipalIntegerMap::new(self.coeff_ring().clone()),
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
