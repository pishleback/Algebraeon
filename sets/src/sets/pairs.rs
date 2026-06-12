use algebraeon_structures::*;
use std::fmt::Debug;

/// The set of Pairs
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PairsStructure<S> {
    set: S,
}

impl<S: SetSignature> PairsStructure<S> {
    pub fn new(set: S) -> Self {
        Self { set }
    }

    /// Construct a new pair from two elements of a set
    pub fn new_pair(&self, a: S::Elem, b: S::Elem) -> Result<(S::Elem, S::Elem), String> {
        self.set.validate_element(&a)?;
        self.set.validate_element(&b)?;
        Ok((a, b))
    }
}

impl<S: SetSignature> Signature for PairsStructure<S> {}

impl<S: SetSignature> SetSignature for PairsStructure<S> {
    type Elem = (S::Elem, S::Elem);

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.set.validate_element(&x.0)?;
        self.set.validate_element(&x.1)?;
        Ok(())
    }
}

impl<S: SetSignature + EqSignature> EqSignature for PairsStructure<S> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set.equal(&a.0, &b.0) && self.set.equal(&a.1, &b.1)
    }
}

/// The set of unordered Pairs of distinct elements
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnorderedPairs<Set> {
    set: Set,
}

impl<Set: SetSignature> UnorderedPairs<Set> {
    pub fn new(set: Set) -> Self {
        Self { set }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct UnorderedPair<T>(T, T);

impl<S: SetSignature + EqSignature> UnorderedPairs<S> {
    pub fn new_pair(&self, a: &S::Elem, b: &S::Elem) -> Result<UnorderedPair<S::Elem>, String> {
        if self.set.equal(a, b) {
            Err("UnorderedPair elements must be distinct".to_string())
        } else {
            Ok(UnorderedPair(a.clone(), b.clone()))
        }
    }
}

impl<S: SetSignature> Signature for UnorderedPairs<S> {}

impl<S: SetSignature + EqSignature> SetSignature for UnorderedPairs<S> {
    type Elem = UnorderedPair<S::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.set.validate_element(&x.0)?;
        self.set.validate_element(&x.1)?;
        if self.set.equal(&x.0, &x.1) {
            Err("UnorderedPair elements must be distinct".to_string())
        } else {
            Ok(())
        }
    }
}

impl<S: SetSignature + EqSignature> EqSignature for UnorderedPairs<S> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set.equal(&a.0, &b.0) && self.set.equal(&a.1, &b.1)
            || self.set.equal(&a.0, &b.1) && self.set.equal(&a.1, &b.0)
    }
}
