use super::*;

pub struct Morphism<Domain: Structure, Range: Structure> {
    domain: Domain,
    range: Range,
}

impl<Domain: Structure, Range: Structure> Morphism<Domain, Range> {
    pub fn domain(&self) -> &Domain {
        &self.domain
    }

    pub fn range(&self) -> &Range {
        &self.range
    }
}

pub trait FunctionStructure<Domain: SetStructure, Range: SetStructure> {
    fn image(&self, x: &Domain::Set) -> Range::Set;
}

pub trait InjectiveFunctionStructure<Domain: SetStructure, Range: SetStructure>:
    FunctionStructure<Domain, Range>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<Domain::Set>;
}
