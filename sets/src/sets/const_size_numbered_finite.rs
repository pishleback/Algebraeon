use algebraeon_structures::*;
use std::fmt::Debug;

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeEnumeratedFiniteSetStructure<const N: usize> {}

impl<const N: usize> ConstSizeEnumeratedFiniteSetStructure<N> {
    pub fn new() -> Self {
        Self {}
    }
}

impl<const N: usize> Signature for ConstSizeEnumeratedFiniteSetStructure<N> {}

impl<const N: usize> SetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    type Elem = usize;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if *x >= N {
            return Err("Too big to be an element".to_string());
        }
        Ok(())
    }
}

impl<const N: usize> EqSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn equal(&self, x: &Self::Elem, y: &Self::Elem) -> bool {
        x == y
    }
}

impl<const N: usize> PartialOrdSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn partial_cmp(&self, x: &Self::Elem, y: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(x, y))
    }
}

impl<const N: usize> OrdSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn cmp(&self, x: &Self::Elem, y: &Self::Elem) -> std::cmp::Ordering {
        x.cmp(y)
    }
}

impl<const N: usize> CountableSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        0..N
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<const N: usize> FiniteSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn size(&self) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize> OrderedFiniteSetSignature for ConstSizeEnumeratedFiniteSetStructure<N> {
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        (0..N).collect()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        Natural::from(*elem)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if let Ok(num) = TryInto::<usize>::try_into(num) {
            if num < N { Some(num) } else { None }
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
