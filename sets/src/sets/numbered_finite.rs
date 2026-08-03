use algebraeon_structures::*;
use std::{fmt::Debug, rc::Rc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnumeratedFiniteSetStructure {
    n: usize,
}

impl EnumeratedFiniteSetStructure {
    pub fn new(n: usize) -> Rc<Self> {
        Self { n }.into()
    }
}

impl Signature for EnumeratedFiniteSetStructure {}

impl SetSignature for EnumeratedFiniteSetStructure {
    type Elem = usize;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if *x >= self.n {
            return Err("Too big to be an element".to_string());
        }
        Ok(())
    }
}

impl EqSignature for EnumeratedFiniteSetStructure {
    fn equal(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl PartialOrdSignature for EnumeratedFiniteSetStructure {
    fn partial_cmp(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(x, y))
    }
}

impl OrdSignature for EnumeratedFiniteSetStructure {
    fn cmp(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl CountableSetSignature for EnumeratedFiniteSetStructure {
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        0..self.n
    }
}

impl FiniteSetSignature for EnumeratedFiniteSetStructure {
    fn size(self: &Rc<Self>) -> Natural {
        Natural::from(self.n)
    }
}

impl OrderedFiniteSetSignature for EnumeratedFiniteSetStructure {
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        (0..self.n).collect()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        Natural::from(*elem)
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        if let Ok(num) = TryInto::<usize>::try_into(num) {
            if num < self.n { Some(num) } else { None }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            EnumeratedFiniteSetStructure::new(123),
            123
        );
    }
}
