use super::{EqSignature, OrdSignature, SetSignature, Signature};
use crate::structure::{CountableSetSignature, FiniteSetSignature, orderings::PartialOrdSignature};
use std::fmt::Debug;

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct SingletonSetStructure {}

impl Signature for SingletonSetStructure {}

impl SetSignature for SingletonSetStructure {
    type Set = ();

    fn validate_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for SingletonSetStructure {
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        x == y
    }
}

impl PartialOrdSignature for SingletonSetStructure {
    fn partial_cmp(&self, a: &Self::Set, b: &Self::Set) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for SingletonSetStructure {
    fn cmp(&self, x: &Self::Set, y: &Self::Set) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for SingletonSetStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        [()].into_iter()
    }
}

impl FiniteSetSignature for SingletonSetStructure {
    fn size(&self) -> usize {
        1
    }
}
