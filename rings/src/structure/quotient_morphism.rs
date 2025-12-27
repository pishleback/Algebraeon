use crate::structure::{
    EuclideanDomainSignature, EuclideanRemainderQuotientStructure, RingHomomorphism,
};
use algebraeon_sets::structure::{Function, Morphism, SetSignature};

#[derive(Clone, Debug)]
pub struct EuclideanDomainQuotientRing<RS, const IS_FIELD: bool>
where
    RS: EuclideanDomainSignature,
{
    source: RS,
    target: EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD>,
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> EuclideanDomainQuotientRing<RS, IS_FIELD> {
    pub fn new(source: RS, target: EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD>) -> Self {
        Self { source, target }
    }
}

impl<RS, const IS_FIELD: bool> Morphism<RS, EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD>>
    for EuclideanDomainQuotientRing<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
    fn domain(&self) -> &RS {
        &self.source
    }

    fn range(&self) -> &EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD> {
        &self.target
    }
}

impl<RS, const IS_FIELD: bool> Function<RS, EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD>>
    for EuclideanDomainQuotientRing<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
    fn image(
        &self,
        x: &<RS as SetSignature>::Set,
    ) -> <EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD> as SetSignature>::Set {
        self.target.reduce(x)
    }
}

impl<RS, const IS_FIELD: bool>
    RingHomomorphism<RS, EuclideanRemainderQuotientStructure<RS, RS, IS_FIELD>>
    for EuclideanDomainQuotientRing<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
}
