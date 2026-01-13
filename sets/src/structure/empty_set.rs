use crate::structure::orderings::PartialOrdSignature;

use super::{EqSignature, OrdSignature, SetSignature, Signature};
use std::fmt::Debug;
use std::marker::PhantomData;

pub struct EmptySetStructure<Set> {
    _set: PhantomData<Set>,
}

impl<Set> Clone for EmptySetStructure<Set> {
    fn clone(&self) -> Self {
        Self { _set: PhantomData }
    }
}

impl<Set> Debug for EmptySetStructure<Set> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EmptySetStructure").finish()
    }
}

impl<Set> PartialEq for EmptySetStructure<Set> {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl<Set> Eq for EmptySetStructure<Set> {}

impl<Set> Default for EmptySetStructure<Set> {
    fn default() -> Self {
        Self { _set: PhantomData }
    }
}

impl<Set: Send + Sync> Signature for EmptySetStructure<Set> {}

impl<Set: Debug + Clone + Send + Sync> SetSignature for EmptySetStructure<Set> {
    type Set = Set;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Err("Empty set has no elements".to_string())
    }
}

impl<Set: Debug + Clone + Send + Sync> EqSignature for EmptySetStructure<Set> {
    fn equal(&self, _: &Self::Set, _: &Self::Set) -> bool {
        panic!("Empty set had no elements to compare for equality")
    }
}

impl<Set: Debug + Clone + Send + Sync> PartialOrdSignature for EmptySetStructure<Set> {
    fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: Debug + Clone + Send + Sync> OrdSignature for EmptySetStructure<Set> {
    fn cmp(&self, _: &Self::Set, _: &Self::Set) -> std::cmp::Ordering {
        panic!("Empty set had no elements to compare for ordering")
    }
}
