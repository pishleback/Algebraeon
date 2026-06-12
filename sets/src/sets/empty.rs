use algebraeon_structures::*;
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Clone)]
pub struct EmptySetStructure<Elem> {
    _set: PhantomData<Elem>,
}

impl<Elem> Debug for EmptySetStructure<Elem> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("EmptySetStructure").finish()
    }
}

impl<Elem> PartialEq for EmptySetStructure<Elem> {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl<Elem> Eq for EmptySetStructure<Elem> {}

impl<Elem> Default for EmptySetStructure<Elem> {
    fn default() -> Self {
        Self { _set: PhantomData }
    }
}

impl<Elem: Clone + Send + Sync> Signature for EmptySetStructure<Elem> {}

impl<Elem: Debug + Clone + Send + Sync> SetSignature for EmptySetStructure<Elem> {
    type Elem = Elem;

    fn validate_element(&self, _: &Self::Elem) -> Result<(), String> {
        Err("Empty set has no elements".to_string())
    }
}

impl<Elem: Debug + Clone + Send + Sync> EqSignature for EmptySetStructure<Elem> {
    fn equal(&self, _: &Self::Elem, _: &Self::Elem) -> bool {
        panic!("Empty set had no elements to compare for equality")
    }
}

impl<Elem: Debug + Clone + Send + Sync> PartialOrdSignature for EmptySetStructure<Elem> {
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Elem: Debug + Clone + Send + Sync> OrdSignature for EmptySetStructure<Elem> {
    fn cmp(&self, _: &Self::Elem, _: &Self::Elem) -> std::cmp::Ordering {
        panic!("Empty set had no elements to compare for ordering")
    }
}
