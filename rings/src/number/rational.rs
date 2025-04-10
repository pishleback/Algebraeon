use crate::polynomial::{PolynomialStructure, factorize_by_factorize_primitive_part};
use crate::structure::*;
use algebraeon_nzq::traits::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

impl SemiRingStructure for RationalCanonicalStructure {
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

impl RingStructure for RationalCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }
}

impl UnitsStructure for RationalCanonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
    }
}

impl IntegralDomainStructure for RationalCanonicalStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if b == &Rational::ZERO {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}

impl CharZeroRingStructure for RationalCanonicalStructure {
    fn try_to_int(&self, x: &Rational) -> Option<Integer> {
        let (n, d) = x.numerator_and_denominator();
        debug_assert_ne!(&d, &Natural::ZERO);
        if d == Natural::ONE { Some(n) } else { None }
    }
}

impl OrderedRingStructure for RationalCanonicalStructure {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FieldStructure for RationalCanonicalStructure {}

impl ComplexSubsetStructure for RationalCanonicalStructure {}

impl RealSubsetStructure for RationalCanonicalStructure {}

impl RealToFloatStructure for RationalCanonicalStructure {
    fn as_f64(&self, x: &Rational) -> f64 {
        let fof = PrincipalSubringInclusion::new(self.clone());
        RealToFloatStructure::as_f64(&Integer::structure(), &fof.numerator(x))
            / RealToFloatStructure::as_f64(&Integer::structure(), &fof.denominator(x))
    }
}

impl FieldOfFractionsInclusion<IntegerCanonicalStructure, RationalCanonicalStructure>
    for PrincipalSubringInclusion<RationalCanonicalStructure>
{
    fn numerator_and_denominator(&self, a: &Rational) -> (Integer, Integer) {
        (a.numerator(), a.denominator().into())
    }
}

impl RealRoundingStructure for RationalCanonicalStructure {
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

impl RealFromFloatStructure for RationalCanonicalStructure {
    fn from_f64_approx(&self, x: f64) -> Self::Set {
        Rational::try_from_float_simplest(x).unwrap()
    }
}

impl FactorableStructure for PolynomialStructure<RationalCanonicalStructure> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
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
