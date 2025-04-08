use super::*;

pub trait MorphismStructure<Domain: Structure, Range: Structure> {
    fn domain(&self) -> &Domain;
    fn range(&self) -> &Range;
}

pub trait FunctionStructure<Domain: SetStructure, Range: SetStructure>:
    MorphismStructure<Domain, Range>
{
    fn image(&self, x: &Domain::Set) -> Range::Set;
}

pub trait InjectiveFunctionStructure<Domain: SetStructure, Range: SetStructure>:
    FunctionStructure<Domain, Range>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<Domain::Set>;
}

pub struct Morphism<Domain: Structure, Range: Structure> {
    domain: Domain,
    range: Range,
}

impl<Domain: Structure, Range: Structure> Morphism<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: Structure, Range: Structure> MorphismStructure<Domain, Range>
    for Morphism<Domain, Range>
{
    fn domain(&self) -> &Domain {
        &self.domain
    }

    fn range(&self) -> &Range {
        &self.range
    }
}
