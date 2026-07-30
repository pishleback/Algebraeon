use crate::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OptionalStructure<Set: SetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> OptionalStructure<Set, SetB> {
    pub fn new(set: SetB) -> Self {
        Self {
            _set: PhantomData,
            set,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> Signature for OptionalStructure<Set, SetB> {}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for OptionalStructure<Set, SetB>
{
    type Elem = Option<Set::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if let Some(x) = x {
            self.set().validate_element(x)
        } else {
            Ok(())
        }
    }
}

impl<Set: EqSignature, SetB: BorrowedStructure<Set>> EqSignature for OptionalStructure<Set, SetB> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        match (a, b) {
            (None, None) => true,
            (None, Some(_)) | (Some(_), None) => false,
            (Some(a), Some(b)) => self.set().equal(a, b),
        }
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for OptionalStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

// take None to be +inf in terms of ordering
impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for OptionalStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        match (a, b) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Greater,
            (Some(_), None) => Ordering::Less,
            (Some(a), Some(b)) => self.set().cmp(a, b),
        }
    }
}

impl<Set: CountableSetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for OptionalStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        [None].into_iter().chain(
            self.set
                .borrow()
                .clone()
                .into_generate_all_elements()
                .map(Some),
        )
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        [None]
            .into_iter()
            .chain(self.set().generate_all_elements().map(Some))
    }
}

impl<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for OptionalStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        self.set().size() + Natural::ONE
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for OptionalStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.set()
            .list_all_elements_ordered()
            .into_iter()
            .map(Some)
            .chain([None])
            .collect()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        match elem {
            Some(elem) => self.set().element_to_enumeration(elem),
            None => self.set().size(),
        }
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if *num == self.set().size() {
            Some(None)
        } else {
            self.set().enumeration_to_element(num).map(Some)
        }
    }
}

impl<Elem: MetaType> MetaType for Option<Elem> {
    type Signature = OptionalStructure<Elem::Signature, Elem::Signature>;

    fn structure() -> Self::Signature {
        OptionalStructure::new(Elem::structure())
    }
}
