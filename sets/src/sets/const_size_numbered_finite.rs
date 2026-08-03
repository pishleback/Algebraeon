use algebraeon_structures::*;
use std::{fmt::Debug, rc::Rc};

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ConstSizeEnumeratedFiniteSet<const N: usize>(usize);

impl<const N: usize> TryFrom<usize> for ConstSizeEnumeratedFiniteSet<N> {
    type Error = ();

    fn try_from(value: usize) -> Result<Self, Self::Error> {
        if value < N { Ok(Self(value)) } else { Err(()) }
    }
}

impl<const N: usize> From<ConstSizeEnumeratedFiniteSet<N>> for usize {
    fn from(value: ConstSizeEnumeratedFiniteSet<N>) -> Self {
        value.0
    }
}

impl<const N: usize> MetaType for ConstSizeEnumeratedFiniteSet<N> {
    type Signature = ConstSizeEnumeratedFiniteSetStructure<N>;

    fn structure() -> Rc<Self::Signature> {
        ConstSizeEnumeratedFiniteSetStructure::new()
    }
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeEnumeratedFiniteSetStructure<const N: usize> {}

impl<const N: usize> ConstSizeEnumeratedFiniteSetStructure<N> {
    pub fn new() -> Rc<Self> {
        Self {}.into()
    }
}

impl<const N: usize> Signature for ConstSizeEnumeratedFiniteSetStructure<N> {}

impl<const N: usize> SetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    type Elem = ConstSizeEnumeratedFiniteSet<N>;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if usize::from(*x) >= N {
            return Err("Too big to be an element".to_string());
        }
        Ok(())
    }
}

impl<const N: usize> EqSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn equal(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl<const N: usize> PartialOrdSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn partial_cmp(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(x, y))
    }
}

impl<const N: usize> OrdSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn cmp(self: &Rc<Self>, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl<const N: usize> CountableSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        (0..N).map(|x| x.try_into().unwrap())
    }
}

impl<const N: usize> FiniteSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn size(self: &Rc<Self>) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize> OrderedFiniteSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        (0..N).map(|x| x.try_into().unwrap()).collect()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        Natural::from(usize::from(*elem))
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        if let Ok(num) = TryInto::<usize>::try_into(num) {
            if num < N {
                Some(num.try_into().unwrap())
            } else {
                None
            }
        } else {
            None
        }
    }
}

impl<const N: usize> ConstSizeFiniteSetSignature<N> for ConstSizeEnumeratedFiniteSetStructure<N> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            ConstSizeEnumeratedFiniteSetStructure::<123>::new(),
            123
        );
    }
}
