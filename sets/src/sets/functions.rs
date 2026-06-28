use algebraeon_structures::*;
use itertools::Itertools;
use std::fmt::Debug;

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Functions<Domain: SetSignature, Range: SetSignature> {
    domain: Domain,
    range: Range,
}

impl<Domain: SetSignature, Range: SetSignature> Functions<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: SetSignature, Range: SetSignature> Signature for Functions<Domain, Range> {}

impl<Domain: FiniteSetSignature, Range: EqSignature> SetSignature for Functions<Domain, Range> {
    type Elem = Vec<Range::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if Natural::from(x.len()) != self.domain.size() {
            return Err("Incorrect vector length".to_string());
        }
        for y in x {
            self.range.validate_element(y)?;
        }
        Ok(())
    }
}

impl<Domain: FiniteSetSignature, Range: EqSignature + FiniteSetSignature> CountableSetSignature
    for Functions<Domain, Range>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.domain.size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.range.list_all_elements())
            .multi_cartesian_product()
    }
}

impl<Domain: FiniteSetSignature, Range: EqSignature + FiniteSetSignature> FiniteSetSignature
    for Functions<Domain, Range>
{
}
