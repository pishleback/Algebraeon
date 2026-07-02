use algebraeon_structures::*;
use std::cmp::Ordering;

/// A sized finite set from an unsized finite set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSetSizedStructure<const N: usize, Set: FiniteSetSignature> {
    // such that self.inner.subset.size() == N
    set: Set,
}

pub trait FiniteSetToFiniteSetSizedSignature: FiniteSetSignature {
    fn try_into_sized<const N: usize>(self) -> Option<FiniteSetSizedStructure<N, Self>> {
        FiniteSetSizedStructure::try_new(self)
    }
}
impl<Set: FiniteSetSignature> FiniteSetToFiniteSetSizedSignature for Set {}

impl<const N: usize, Set: FiniteSetSignature> FiniteSetSizedStructure<N, Set> {
    pub fn forget_sized(&self) -> &Set {
        &self.set
    }

    pub fn into_forget_sized(self) -> Set {
        self.set
    }

    pub fn try_new(set: Set) -> Option<Self> {
        if set.size() == Natural::from(N) {
            Some(Self { set })
        } else {
            None
        }
    }
}

impl<const N: usize, Set: FiniteSetSignature> Signature for FiniteSetSizedStructure<N, Set> {}

impl<const N: usize, Set: FiniteSetSignature> SetSignature for FiniteSetSizedStructure<N, Set> {
    type Elem = Set::Elem;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.forget_sized().validate_element(x)?;
        Ok(())
    }
}

impl<const N: usize, Set: FiniteSetSignature + EqSignature> EqSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.forget_sized().equal(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature + PartialOrdSignature> PartialOrdSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.forget_sized().partial_cmp(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature + OrdSignature> OrdSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.forget_sized().cmp(a, b)
    }
}

impl<const N: usize, Set: FiniteSetSignature> CountableSetSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.into_forget_sized().into_generate_all_elements()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.forget_sized().generate_all_elements()
    }
}

impl<const N: usize, Set: FiniteSetSignature> FiniteSetSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn size(&self) -> Natural {
        #[cfg(debug_assertions)]
        {
            let n = self.forget_sized().size();
            assert_eq!(n, Natural::from(N));
        }
        Natural::from(N)
    }
}

impl<const N: usize, Set: FiniteSetSignature> FiniteSetSizedSignature<N>
    for FiniteSetSizedStructure<N, Set>
{
}

impl<const N: usize, Set: EnumeratedOrdFiniteSetSignature> EnumeratedOrdFiniteSetSignature
    for FiniteSetSizedStructure<N, Set>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.forget_sized().list_all_elements_ordered()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.forget_sized().element_to_enumeration(elem)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.forget_sized().enumeration_to_element(num)
    }
}
