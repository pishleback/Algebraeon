use super::*;
use algebraeon_sets::structure::*;
use std::rc::Rc;

pub trait RingHomomorphismStructure<Domain: RingStructure, Range: RingStructure>:
    FunctionStructure<Domain, Range>
{
}

pub trait FieldOfFractionsStructure<RS: IntegralDomainStructure>: FieldStructure {
    fn base_ring_structure(&self) -> Rc<RS>;
    fn from_base_ring(&self, elem: <RS as SetStructure>::Set) -> Self::Set;
    fn numerator(&self, elem: &Self::Set) -> <RS as SetStructure>::Set;
    fn denominator(&self, elem: &Self::Set) -> <RS as SetStructure>::Set;
    fn as_base_ring(&self, elem: Self::Set) -> Option<<RS as SetStructure>::Set> {
        let base_ring = self.base_ring_structure();
        if base_ring.equal(&self.denominator(&elem), &base_ring.one()) {
            Some(self.numerator(&elem))
        } else {
            None
        }
    }
}
pub trait MetaFieldOfFractions<RS: IntegralDomainStructure>: MetaRing
where
    Self::Structure: FieldOfFractionsStructure<RS>,
{
    fn base_ring_structure() -> Rc<RS> {
        Self::structure().base_ring_structure()
    }
    fn from_base_ring(a: <RS as SetStructure>::Set) -> Self {
        Self::structure().from_base_ring(a)
    }
    fn numerator(&self) -> <RS as SetStructure>::Set {
        Self::structure().numerator(self)
    }
    fn denominator(&self) -> <RS as SetStructure>::Set {
        Self::structure().denominator(self)
    }
    fn as_base_ring(self) -> Option<<RS as SetStructure>::Set> {
        Self::structure().as_base_ring(self)
    }
}
impl<R: MetaRing, RS: IntegralDomainStructure> MetaFieldOfFractions<RS> for R where
    Self::Structure: FieldOfFractionsStructure<RS, Set = R>
{
}
