use algebraeon_structures::*;
use std::{fmt::Debug, sync::Arc};

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct SingletonSetStructure {}

impl SingletonSetStructure {
    pub fn new() -> Arc<Self> {
        Self {}.into()
    }
}

impl Signature for SingletonSetStructure {}

impl SetSignature for SingletonSetStructure {
    type Elem = ();

    fn validate_element(self: &Arc<Self>, _: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for SingletonSetStructure {
    fn equal(self: &Arc<Self>, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl PartialOrdSignature for SingletonSetStructure {
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for SingletonSetStructure {
    fn cmp(self: &Arc<Self>, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for SingletonSetStructure {
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        [()].into_iter()
    }
}

impl FiniteSetSignature for SingletonSetStructure {
    fn size(self: &Arc<Self>) -> Natural {
        Natural::ONE
    }
}
