use crate::*;
use std::{cmp::Ordering, rc::Rc};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OptionalStructure<Set: SetSignature> {
    set: Rc<Set>,
}

impl<Set: SetSignature> OptionalStructure<Set> {
    pub fn new(set: Rc<Set>) -> Rc<Self> {
        Self { set }.into()
    }

    pub fn set(self: &Rc<Self>) -> &Rc<Set> {
        &self.set
    }
}

impl<Set: SetSignature> Signature for OptionalStructure<Set> {}

impl<Set: SetSignature> SetSignature for OptionalStructure<Set> {
    type Elem = Option<Set::Elem>;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if let Some(x) = x {
            self.set().validate_element(x)
        } else {
            Ok(())
        }
    }
}

impl<Set: EqSignature> EqSignature for OptionalStructure<Set> {
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        match (a, b) {
            (None, None) => true,
            (None, Some(_)) | (Some(_), None) => false,
            (Some(a), Some(b)) => self.set().equal(a, b),
        }
    }
}

impl<Set: OrdSignature> PartialOrdSignature for OptionalStructure<Set> {
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

// take None to be +inf in terms of ordering
impl<Set: OrdSignature> OrdSignature for OptionalStructure<Set> {
    fn cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Greater,
            (Some(_), None) => Ordering::Less,
            (Some(a), Some(b)) => self.set().cmp(a, b),
        }
    }
}

impl<Set: CountableSetSignature> CountableSetSignature for OptionalStructure<Set> {
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        [None]
            .into_iter()
            .chain(self.set().clone().generate_all_elements().map(Some))
    }
}

impl<Set: FiniteSetSignature> FiniteSetSignature for OptionalStructure<Set> {
    fn size(self: &Rc<Self>) -> Natural {
        self.set().size() + Natural::ONE
    }
}

impl<Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature for OptionalStructure<Set> {
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        self.set()
            .list_all_elements_ordered()
            .into_iter()
            .map(Some)
            .chain([None])
            .collect()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        match elem {
            Some(elem) => self.set().element_to_enumeration(elem),
            None => self.set().size(),
        }
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        if *num == self.set().size() {
            Some(None)
        } else {
            self.set().enumeration_to_element(num).map(Some)
        }
    }
}

impl<Elem: MetaType> MetaType for Option<Elem> {
    type Signature = OptionalStructure<Elem::Signature>;

    fn structure() -> Rc<Self::Signature> {
        OptionalStructure::new(Elem::structure())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn optional_enumeration() {
        assert_enumerated_ord_finite_set!(
            OptionalStructure::new(i32::structure().const_size_finite_subset([1, 2, 3, 4, 5])),
            6
        );
    }
}
