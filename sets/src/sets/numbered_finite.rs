use algebraeon_structures::*;
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
    type Elem = usize;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if x >= &self.n {
            return Err("Too big to be an element".to_string());
        }
        Ok(())
    }
}

impl EqSignature for EnumeratedFiniteSetStructure {
    fn equal(&self, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl PartialOrdSignature for EnumeratedFiniteSetStructure {
    fn partial_cmp(&self, x: &Self::Elem, y: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(x, y))
    }
}

impl OrdSignature for EnumeratedFiniteSetStructure {
    fn cmp(&self, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for EnumeratedFiniteSetStructure {
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        0..self.n
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl FiniteSetSignature for EnumeratedFiniteSetStructure {
    fn size(&self) -> Natural {
        Natural::from(self.n)
    }
}

impl EnumeratedOrdFiniteSetSignature for EnumeratedFiniteSetStructure {
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        (0..self.n).collect()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        Natural::from(*elem)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if let Ok(num) = TryInto::<usize>::try_into(num) {
            if num < self.n { Some(num) } else { None }
        } else {
            None
        }
    }
}
