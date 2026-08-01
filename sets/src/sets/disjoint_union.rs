use algebraeon_structures::*;
use std::sync::Arc;

/// A sized finite set from an unsized finite set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DisjointUnionSetStructure<Set0: SetSignature, Set1: SetSignature> {
    set_0: Arc<Set0>,
    set_1: Arc<Set1>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DisjointUnionElem2<Elem0, Elem1> {
    Elem0(Elem0),
    Elem1(Elem1),
}

impl<Set0: SetSignature, Set1: SetSignature> DisjointUnionSetStructure<Set0, Set1> {
    pub fn new(set_0: Arc<Set0>, set_1: Arc<Set1>) -> Arc<Self> {
        Self { set_0, set_1 }.into()
    }

    pub fn set_0(&self) -> &Arc<Set0> {
        &self.set_0
    }

    pub fn set_1(&self) -> &Arc<Set1> {
        &self.set_1
    }
}

impl<Set0: SetSignature, Set1: SetSignature> Signature for DisjointUnionSetStructure<Set0, Set1> {}

impl<Set0: SetSignature, Set1: SetSignature> SetSignature
    for DisjointUnionSetStructure<Set0, Set1>
{
    type Elem = DisjointUnionElem2<Set0::Elem, Set1::Elem>;

    fn validate_element(self: &Arc<Self>, elem: &Self::Elem) -> Result<(), String> {
        match elem {
            DisjointUnionElem2::Elem0(elem) => self.set_0().validate_element(elem),
            DisjointUnionElem2::Elem1(elem) => self.set_1().validate_element(elem),
        }
    }
}

impl<Set0: EqSignature, Set1: EqSignature> EqSignature for DisjointUnionSetStructure<Set0, Set1> {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
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

    #[test]
    fn test() {
        let set_0 = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_1 = i32::structure().into_const_size_finite_subset([6, 7, 8]);

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
