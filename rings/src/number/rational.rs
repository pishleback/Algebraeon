use crate::polynomial::{PolynomialStructure, factorize_by_factorize_primitive_part};
use crate::structure::*;
use algebraeon_nzq::traits::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

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

impl UnitsStructure for CannonicalStructure<Rational> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
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

impl CharZeroStructure for CannonicalStructure<Rational> {
    fn try_to_int(&self, x: &Rational) -> Option<Integer> {
        let (n, d) = x.numerator_and_denominator();
        debug_assert_ne!(&d, &Natural::ZERO);
        if d == Natural::ONE { Some(n) } else { None }
    }
}

impl OrderedRingStructure for CannonicalStructure<Rational> {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        Self::Set::cmp(a, b)
    }
}

impl FieldStructure for CannonicalStructure<Rational> {}

impl ComplexSubsetStructure for CannonicalStructure<Rational> {}

impl RealSubsetStructure for CannonicalStructure<Rational> {}

impl RealToFloatStructure for CannonicalStructure<Rational> {
    fn as_f64(&self, x: &Rational) -> f64 {
        let fof = PrincipalSubringInclusion::new(self.clone());
        RealToFloatStructure::as_f64(Integer::structure().as_ref(), &fof.numerator(x))
            / RealToFloatStructure::as_f64(Integer::structure().as_ref(), &fof.denominator(x))
    }
}

impl FieldOfFractionsInclusionStructure<CannonicalStructure<Integer>, CannonicalStructure<Rational>>
    for PrincipalSubringInclusion<CannonicalStructure<Rational>>
{
    fn numerator_and_denominator(&self, a: &Rational) -> (Integer, Integer) {
        (a.numerator(), a.denominator().into())
    }
}

impl RealRoundingStructure for CannonicalStructure<Rational> {
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

impl RealFromFloatStructure for CannonicalStructure<Rational> {
    fn from_f64_approx(&self, x: f64) -> Self::Set {
        Rational::try_from_float_simplest(x).unwrap()
    }
}

impl FactorableStructure for PolynomialStructure<CannonicalStructure<Rational>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
        factorize_by_factorize_primitive_part(
            &PrincipalSubringInclusion::new(self.coeff_ring().as_ref().clone()),
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
