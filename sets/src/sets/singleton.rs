use algebraeon_structures::*;
use std::{fmt::Debug, rc::Rc};

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct SingletonSetStructure {}

impl SingletonSetStructure {
    pub fn new() -> Rc<Self> {
        Self {}.into()
    }
}

impl Signature for SingletonSetStructure {}

impl SetSignature for SingletonSetStructure {
    type Elem = ();

    fn validate_element(self: &Rc<Self>, _: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for SingletonSetStructure {
    fn equal(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl PartialOrdSignature for SingletonSetStructure {
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for SingletonSetStructure {
    fn cmp(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for SingletonSetStructure {
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        [()].into_iter()
    }
}

impl FiniteSetSignature for SingletonSetStructure {
    fn size(self: &Rc<Self>) -> Natural {
        Natural::ONE
    }
}
