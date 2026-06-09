use algebraeon_structures::Signature;
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
        if x.len() != self.domain.size() {
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
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        (0..self.domain.size())
            .map(|_| self.range.list_all_elements())
            .multi_cartesian_product()
    }
}

impl<Domain: FiniteSetSignature, Range: EqSignature + FiniteSetSignature> FiniteSetSignature
    for Functions<Domain, Range>
{
}

/// The set of all endofunctions on a finite set X: functions X → X
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSetEndofunctions<X: FiniteSetSignature + EqSignature> {
    set: X,
}

impl<X: FiniteSetSignature + EqSignature> FiniteSetEndofunctions<X> {
    pub fn new(set: X) -> Self {
        Self { set }
    }
}

impl<X: FiniteSetSignature + EqSignature> Signature for FiniteSetEndofunctions<X> {}

impl<X: FiniteSetSignature + EqSignature> SetSignature for FiniteSetEndofunctions<X> {
    type Elem = Vec<X::Elem>;

    fn validate_element(&self, f: &Self::Elem) -> Result<(), String> {
        if f.len() != self.set.size() {
            return Err("Function must have one value per element in the domain.".to_string());
        }
        for y in f {
            self.set.validate_element(y)?;
        }
        Ok(())
    }
}

impl<X: FiniteSetSignature + EqSignature> CountableSetSignature for FiniteSetEndofunctions<X> {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        (0..self.set.size())
            .map(|_| self.set.list_all_elements())
            .multi_cartesian_product()
    }
}

impl<X: FiniteSetSignature + EqSignature> FiniteSetSignature for FiniteSetEndofunctions<X> {}
