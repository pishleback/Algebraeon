use crate::*;
use std::{cmp::Ordering, sync::Arc};

/// A sized finite set from an unsized finite set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeFiniteSetStructure<const N: usize, Set: FiniteSetSignature> {
    // such that self.inner.subset.size() == N
    set: Arc<Set>,
}

pub trait FiniteSetToFiniteSetSizedSignature: FiniteSetSignature {
    fn try_to_const_sized<const N: usize>(
        self: &Arc<Self>,
    ) -> Option<Arc<ConstSizeFiniteSetStructure<N, Self>>> {
        ConstSizeFiniteSetStructure::try_new(self.clone())
    }
}
impl<Set: FiniteSetSignature> FiniteSetToFiniteSetSizedSignature for Set {}

impl<const N: usize, Set: FiniteSetSignature> ConstSizeFiniteSetStructure<N, Set> {
    pub fn forget_const_sized(self: &Arc<Self>) -> &Arc<Set> {
        &self.set
    }

    pub fn try_new(set: Arc<Set>) -> Option<Arc<Self>> {
        if set.size() == Natural::from(N) {
            Some(Self { set }.into())
        } else {
            None
        }
    }
}

impl<const N: usize, Set: FiniteSetSignature> Signature for ConstSizeFiniteSetStructure<N, Set> {}

impl<const N: usize, Set: FiniteSetSignature> SetSignature for ConstSizeFiniteSetStructure<N, Set> {
    type Elem = Set::Elem;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.forget_const_sized().validate_element(x)?;
        Ok(())
    }
}

impl<const N: usize, Set: FiniteSetSignature + EqSignature> EqSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.forget_const_sized().equal(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature + PartialOrdSignature> PartialOrdSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.forget_const_sized().partial_cmp(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature + OrdSignature> OrdSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.forget_const_sized().cmp(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature> CountableSetSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.forget_const_sized().clone().generate_all_elements()
    }
}

impl<const N: usize, Set: FiniteSetSignature> FiniteSetSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn size(self: &Arc<Self>) -> Natural {
        #[cfg(debug_assertions)]
        {
            let n = self.forget_const_sized().size();
            assert_eq!(n, Natural::from(N));
        }
        Natural::from(N)
    }
}

impl<const N: usize, Set: FiniteSetSignature> ConstSizeFiniteSetSignature<N>
    for ConstSizeFiniteSetStructure<N, Set>
{
}

impl<const N: usize, Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for ConstSizeFiniteSetStructure<N, Set>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        self.forget_const_sized().list_all_elements_ordered()
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        self.forget_const_sized().element_to_enumeration(elem)
    }

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
        self.forget_const_sized().enumeration_to_element(num)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumeration() {
        assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![1, 2, 3, 4, 5])
                .try_to_const_sized::<5>()
                .unwrap(),
            5
        );
    }
}
