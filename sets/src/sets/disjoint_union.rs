use algebraeon_structures::*;
use std::marker::PhantomData;

/// A sized finite set from an unsized finite set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DisjointUnionSetStructure<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> {
    _set_0: PhantomData<Set0>,
    set_0: Set0B,
    _set_1: PhantomData<Set1>,
    set_1: Set1B,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DisjointUnionElem2<Elem0, Elem1> {
    Elem0(Elem0),
    Elem1(Elem1),
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> DisjointUnionSetStructure<Set0, Set0B, Set1, Set1B>
{
    pub fn new(set_0: Set0B, set_1: Set1B) -> Self {
        Self {
            _set_0: PhantomData,
            set_0,
            _set_1: PhantomData,
            set_1,
        }
    }

    pub fn set_0(&self) -> &Set0 {
        self.set_0.borrow()
    }

    pub fn set_1(&self) -> &Set1 {
        self.set_1.borrow()
    }
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> Signature for DisjointUnionSetStructure<Set0, Set0B, Set1, Set1B>
{
}

impl<
    Set0: SetSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: SetSignature,
    Set1B: BorrowedStructure<Set1>,
> SetSignature for DisjointUnionSetStructure<Set0, Set0B, Set1, Set1B>
{
    type Elem = DisjointUnionElem2<Set0::Elem, Set1::Elem>;

    fn validate_element(&self, elem: &Self::Elem) -> Result<(), String> {
        match elem {
            DisjointUnionElem2::Elem0(elem) => self.set_0().validate_element(elem),
            DisjointUnionElem2::Elem1(elem) => self.set_1().validate_element(elem),
        }
    }
}

impl<
    Set0: EqSignature,
    Set0B: BorrowedStructure<Set0>,
    Set1: EqSignature,
    Set1B: BorrowedStructure<Set1>,
> EqSignature for DisjointUnionSetStructure<Set0, Set0B, Set1, Set1B>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        match (a, b) {
            (DisjointUnionElem2::Elem0(a), DisjointUnionElem2::Elem0(b)) => {
                self.set_0().equal(a, b)
            }
            (DisjointUnionElem2::Elem1(a), DisjointUnionElem2::Elem1(b)) => {
                self.set_1().equal(a, b)
            }
            (DisjointUnionElem2::Elem0(_), DisjointUnionElem2::Elem1(_))
            | (DisjointUnionElem2::Elem1(_), DisjointUnionElem2::Elem0(_)) => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sets::SetToFiniteSubsetByOrdSizedSignature;

    #[test]
    fn test() {
        let set_0 = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5]);
        let set_1 = i32::structure().into_finite_subset_sized([6, 7, 8]);

        let set_01 = DisjointUnionSetStructure::new(set_0, set_1);

        assert!(set_01.is_element(&DisjointUnionElem2::Elem0(1)));
        assert!(set_01.is_element(&DisjointUnionElem2::Elem0(5)));
        assert!(!set_01.is_element(&DisjointUnionElem2::Elem0(6)));
        assert!(!set_01.is_element(&DisjointUnionElem2::Elem0(8)));
        assert!(!set_01.is_element(&DisjointUnionElem2::Elem1(1)));
        assert!(!set_01.is_element(&DisjointUnionElem2::Elem1(5)));
        assert!(set_01.is_element(&DisjointUnionElem2::Elem1(6)));
        assert!(set_01.is_element(&DisjointUnionElem2::Elem1(8)));
    }
}
