use algebraeon_structures::*;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::rc::Rc;

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

impl<Elem: Clone> Signature for EmptySetStructure<Elem> {}

impl<Elem: Debug + Clone> SetSignature for EmptySetStructure<Elem> {
    type Elem = Elem;

    fn validate_element(self: &Rc<Self>, _: &Self::Elem) -> Result<(), String> {
        Err("Empty set has no elements".to_string())
    }
}

impl<Elem: Debug + Clone> EqSignature for EmptySetStructure<Elem> {
    fn equal(self: &Rc<Self>, _: &Self::Elem, _: &Self::Elem) -> bool {
        panic!("Empty set had no elements to compare for equality")
    }
}

impl<Elem: Debug + Clone> PartialOrdSignature for EmptySetStructure<Elem> {
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Elem: Debug + Clone> OrdSignature for EmptySetStructure<Elem> {
    fn cmp(self: &Rc<Self>, _: &Self::Elem, _: &Self::Elem) -> std::cmp::Ordering {
        panic!("Empty set had no elements to compare for ordering")
    }
}
