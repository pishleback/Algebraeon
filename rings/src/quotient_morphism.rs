use algebraeon_sets::structure::{Function, Morphism, SetSignature};

use crate::{
    rings::quotient::QuotientStructure,
    structure::{EuclideanDomainSignature, RingHomomorphism},
};

#[derive(Clone, Debug)]
pub struct EuclideanDomainQuotienting<RS, const IS_FIELD: bool>
where
    RS: EuclideanDomainSignature,
{
    source: RS,
    target: QuotientStructure<RS, IS_FIELD>,
}

impl<RS: EuclideanDomainSignature, const IS_FIELD: bool> EuclideanDomainQuotienting<RS, IS_FIELD> {
    pub fn new(source: RS, target: QuotientStructure<RS, IS_FIELD>) -> Self {
        Self { source, target }
    }
}

impl<RS, const IS_FIELD: bool> Morphism<RS, QuotientStructure<RS, IS_FIELD>>
    for EuclideanDomainQuotienting<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
    fn domain(&self) -> &RS {
        &self.source
    }

    fn range(&self) -> &QuotientStructure<RS, IS_FIELD> {
        &self.target
    }
}

impl<RS, const IS_FIELD: bool> Function<RS, QuotientStructure<RS, IS_FIELD>>
    for EuclideanDomainQuotienting<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
    fn image(
        &self,
        x: &<RS as SetSignature>::Set,
    ) -> <QuotientStructure<RS, IS_FIELD> as SetSignature>::Set {
        self.target.reduce(x)
    }
}

impl<RS, const IS_FIELD: bool> RingHomomorphism<RS, QuotientStructure<RS, IS_FIELD>>
    for EuclideanDomainQuotienting<RS, IS_FIELD>
where
    RS: EuclideanDomainSignature,
{
}
