use malachite_base::num::basic::traits::One;
use malachite_base::num::basic::traits::Zero;
use malachite_nz::integer::Integer;
use malachite_q::Rational;

use super::super::super::structure::*;

use super::super::ring_structure::cannonical::*;
use super::super::ring_structure::factorization::*;
use super::super::ring_structure::structure::*;

impl StructuredType for Rational {
    type Structure = CannonicalStructure<Self>;

    fn structure() -> Self::Structure {
        Self::Structure::new()
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

    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if b == &Rational::ZERO {
            Err(RingDivisionError::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}

impl IntegralDomainStructure for CannonicalStructure<Rational> {}

impl FieldStructure for CannonicalStructure<Rational> {}

impl FieldOfFractionsStructure for CannonicalStructure<Rational> {
    type RS = CannonicalStructure<Integer>;

    fn base_ring_structure(&self) -> Self::RS {
        Integer::structure()
    }

    fn from_base_ring(&self, elem: <Self::RS as Structure>::Set) -> Self::Set {
        Rational::from(elem)
    }

    fn numerator(&self, elem: &Self::Set) -> <Self::RS as Structure>::Set {
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
