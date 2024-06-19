use std::rc::Rc;

use malachite_base::num::arithmetic::traits::Floor;
use malachite_base::num::basic::traits::One;
use malachite_base::num::basic::traits::OneHalf;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::integer::Integer;
use malachite_q::Rational;

use crate::rings::structure::*;

use super::super::polynomial::polynomial::*;

use super::super::ring_structure::cannonical::*;
use super::super::ring_structure::factorization::*;
use super::super::ring_structure::structure::*;
use super::algebraic::isolated_roots::ComplexAlgebraic;
use super::algebraic::isolated_roots::RealAlgebraic;

impl StructuredType for Rational {
    type Structure = CannonicalStructure<Self>;

    fn structure() -> Rc<Self::Structure> {
        Self::Structure::new().into()
    }
}

impl EqualityStructure for CannonicalStructure<Rational> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl RingStructure for CannonicalStructure<Rational> {
    fn zero(&self) -> Self::Set {
        Rational::ZERO
    }

    fn one(&self) -> Self::Set {
        Rational::ONE
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        -a
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
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
        //malachite returns a natural for the numerator for some
        if elem >= &0 {
            Integer::from(elem.numerator_ref())
        } else {
            -Integer::from(elem.numerator_ref())
        }
    }

    fn denominator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set {
        Integer::from(elem.denominator_ref())
    }
}

impl RealRoundingStructure for CannonicalStructure<Rational> {
    fn floor(&self, x: &Self::Set) -> Integer {
        <Rational as malachite_base::num::arithmetic::traits::Floor>::floor(x.clone())
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
        Rational::try_from_float_simplest(x).unwrap()
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Rational>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
        self.factorize_by_factorize_primitive_part(p)
    }
}

impl AlgebraicClosureStructure for CannonicalStructure<Rational> {
    type ACFS = CannonicalStructure<ComplexAlgebraic>;

    fn algebraic_closure_field(&self) -> Rc<Self::ACFS> {
        ComplexAlgebraic::structure()
    }

    fn algebraic_closure_inclusion(&self, x: &Self::Set) -> <Self::ACFS as Structure>::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(x.clone()))
    }

    fn all_roots_list(
        &self,
        poly: &Polynomial<Self::Set>,
    ) -> Option<Vec<<Self::ACFS as Structure>::Set>> {
        if poly.is_zero() {
            None
        } else {
            Some(poly.primitive_part_fof().all_complex_roots())
        }
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::super::super::ring_structure::cannonical::*;
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
