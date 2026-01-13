use super::{EqSignature, OrdSignature, SetSignature, Signature};
use crate::structure::{CountableSetSignature, FiniteSetSignature, orderings::PartialOrdSignature};
use std::fmt::Debug;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnumeratedFiniteSetStructure {
    n: usize,
}

impl EnumeratedFiniteSetStructure {
    pub fn new(n: usize) -> Self {
        Self { n }
    }
}

impl Signature for EnumeratedFiniteSetStructure {}

impl SetSignature for EnumeratedFiniteSetStructure {
    type Set = usize;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if x >= &self.n {
            return Err("Too big to be an element".to_string());
        }
        Ok(())
    }
}

impl EqSignature for EnumeratedFiniteSetStructure {
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        x == y
    }
}

impl PartialOrdSignature for EnumeratedFiniteSetStructure {
    fn partial_cmp(&self, x: &Self::Set, y: &Self::Set) -> Option<std::cmp::Ordering> {
        Some(self.cmp(x, y))
    }
}

impl OrdSignature for EnumeratedFiniteSetStructure {
    fn cmp(&self, x: &Self::Set, y: &Self::Set) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for EnumeratedFiniteSetStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        0..self.n
    }
}

impl FiniteSetSignature for EnumeratedFiniteSetStructure {
    fn size(&self) -> usize {
        self.n
    }
}
