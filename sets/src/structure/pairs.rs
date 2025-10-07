use crate::structure::Signature;

use super::{EqSignature, SetSignature};
use std::fmt::Debug;

/// The set of Pairs
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pairs<T> {
    set: T,
}

impl<S: SetSignature> Pairs<S> {
    pub fn new(set: S) -> Self {
        Self { set }
    }

    /// Construct a new pair from two elements of a set
    /// # Errors
    /// If either `a` or `b` is not an element of `S`, then `(a,b)` is not one of `Pairs<S>`
    pub fn new_pair(&self, a: S::Set, b: S::Set) -> Result<(S::Set, S::Set), String> {
        self.set.is_element(&a)?;
        self.set.is_element(&b)?;
        Ok((a, b))
    }
}

impl<S: SetSignature> Signature for Pairs<S> {}

impl<S: SetSignature> SetSignature for Pairs<S> {
    type Set = (S::Set, S::Set);

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.set.is_element(&x.0)?;
        self.set.is_element(&x.1)?;
        Ok(())
    }
}

impl<S: SetSignature + EqSignature> EqSignature for Pairs<S> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
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
    /// # Errors
    /// `a` and `b` are expected to not be equal as elements of `S`
    pub fn new_pair(&self, a: &S::Set, b: &S::Set) -> Result<UnorderedPair<S::Set>, String> {
        if self.set.equal(a, b) {
            Err("UnorderedPair elements must be distinct".to_string())
        } else {
            Ok(UnorderedPair(a.clone(), b.clone()))
        }
    }
}

impl<S: SetSignature> Signature for UnorderedPairs<S> {}

impl<S: SetSignature + EqSignature> SetSignature for UnorderedPairs<S> {
    type Set = UnorderedPair<S::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.set.is_element(&x.0)?;
        self.set.is_element(&x.1)?;
        if self.set.equal(&x.0, &x.1) {
            Err("UnorderedPair elements must be distinct".to_string())
        } else {
            Ok(())
        }
    }
}

impl<S: SetSignature + EqSignature> EqSignature for UnorderedPairs<S> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.set.equal(&a.0, &b.0) && self.set.equal(&a.1, &b.1)
            || self.set.equal(&a.0, &b.1) && self.set.equal(&a.1, &b.0)
    }
}
