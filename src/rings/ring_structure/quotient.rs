use std::rc::Rc;

use super::super::structure::*;

use super::structure::*;

#[derive(Debug, Clone)]
pub struct QuotientStructure<
    RS: EuclideanDivisionStructure + UniqueFactorizationStructure,
    const IS_FIELD: bool,
> {
    ring: Rc<RS>,
    modulus: RS::Set,
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool>
    QuotientStructure<RS, IS_FIELD>
{
    pub fn new(ring: Rc<RS>, modulus: RS::Set) -> Self {
        assert!(!ring.is_zero(&modulus));
        if IS_FIELD {
            if !ring.is_irreducible(&modulus) {
                panic!(
                    "The modulus must be irreducible to form a quotient field. Got {:?}.",
                    modulus
                );
            }
        }
        Self { ring, modulus }
    }

    pub fn ring(&self) -> Rc<RS> {
        self.ring.clone()
    }

    pub fn modulus(&self) -> &RS::Set {
        &self.modulus
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure> QuotientStructure<RS, false> {
    pub fn new_ring(ring: Rc<RS>, modulus: RS::Set) -> Self {
        Self::new(ring, modulus)
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure> QuotientStructure<RS, true> {
    pub fn new_field(ring: Rc<RS>, modulus: RS::Set) -> Self {
        Self::new(ring, modulus)
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool> PartialEq
    for QuotientStructure<RS, IS_FIELD>
{
    fn eq(&self, other: &Self) -> bool {
        if self.ring == other.ring {
            debug_assert_eq!(
                self.ring.equal(&self.modulus, &other.modulus),
                other.ring.equal(&self.modulus, &other.modulus)
            );
            self.ring.equal(&self.modulus, &other.modulus)
        } else {
            false
        }
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool> Eq
    for QuotientStructure<RS, IS_FIELD>
{
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool> Structure
    for QuotientStructure<RS, IS_FIELD>
{
    type Set = RS::Set;
}

impl<
        RS: EuclideanDivisionStructure + UniqueFactorizationStructure + DisplayableStructure,
        const IS_FIELD: bool,
    > DisplayableStructure for QuotientStructure<RS, IS_FIELD>
{
    fn elem_to_string(&self, elem: &Self::Set) -> String {
        self.ring.elem_to_string(elem)
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool>
    EqualityStructure for QuotientStructure<RS, IS_FIELD>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.ring.is_zero(
            &self
                .ring
                .rem(&self.ring.add(a, &self.ring.neg(b)), &self.modulus),
        )
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure, const IS_FIELD: bool>
    RingStructure for QuotientStructure<RS, IS_FIELD>
{
    fn zero(&self) -> Self::Set {
        self.ring.zero()
    }

    fn one(&self) -> Self::Set {
        self.ring.one()
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.ring.neg(a)
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring.rem(&self.ring.add(a, b), &self.modulus)
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ring.rem(&self.ring.mul(a, b), &self.modulus)
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure> IntegralDomainStructure
    for QuotientStructure<RS, true>
{
    fn div(
        &self,
        top: &Self::Set,
        bot: &Self::Set,
    ) -> Result<Self::Set, super::structure::RingDivisionError> {
        if self.is_zero(bot) {
            Err(RingDivisionError::DivideByZero)
        } else {
            let (g, a, b) = self.ring.euclidean_xgcd(bot.clone(), self.modulus.clone());
            debug_assert!(self.ring.equal(
                &g,
                &self
                    .ring
                    .add(&self.ring.mul(&a, &bot), &self.ring.mul(&b, &self.modulus))
            ));
            //g = a * bot mod N and g divides N
            //want ?*bot = top
            match self.ring.div(top, &g) {
                Ok(d) => Ok(self.ring.rem(&self.ring.mul(&a, &d), &self.modulus)),
                Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
                Err(RingDivisionError::DivideByZero) => unreachable!(),
            }
        }
    }
}

impl<RS: EuclideanDivisionStructure + UniqueFactorizationStructure> FieldStructure
    for QuotientStructure<RS, true>
{
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use super::super::super::number::rational::*;
    use super::super::super::ring_structure::cannonical::*;
    use super::super::super::ring_structure::structure::*;

    use super::*;

    #[test]
    fn test() {
        let mod5 = QuotientStructure::new_field(Integer::structure(), Integer::from(5));

        assert!(mod5.equal(
            &mod5.div(&Integer::from(3), &Integer::from(2)).unwrap(),
            &Integer::from(9)
        ));
    }
}
