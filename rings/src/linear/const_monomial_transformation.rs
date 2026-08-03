//! For a ring R, a monomial transformation of R^n is a composition of a permutation of coordinates with a scalar multiplication for each coordinate

use algebraeon_sets::sets::{ConstSizePermutation, SetToConstSizePermutationsStructure};
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ConstSizeMonomialTransformation<const N: usize, Elem> {
    // represents the opperation of multiplying by the scalars component-wise followed by the permutation of coordinates
    // note that the order matters here
    scalars: [Elem; N],
    permutation: ConstSizePermutation<N, Elem>,
}

impl<const N: usize, Elem> ConstSizeMonomialTransformation<N, Elem> {
    pub fn validate(&self) -> Result<(), String> {
        self.permutation.validate()?;
        Ok(())
    }
}

impl<const N: usize, Elem: MetaType> MetaType for ConstSizeMonomialTransformation<N, Elem>
where
    Elem::Signature: ConstSizeFiniteSetSignature<N>,
{
    type Signature = ConstSizeMonomialTransformationsStructure<N, Elem::Signature>;

    fn structure() -> Arc<Self::Signature> {
        ConstSizeMonomialTransformationsStructure::new(Elem::structure())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeMonomialTransformationsStructure<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N>,
> {
    set: Arc<Set>,
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>>
    ConstSizeMonomialTransformationsStructure<N, Set>
{
    pub fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }

    pub fn set(&self) -> &Arc<Set> {
        &self.set
    }
}

pub trait SetToConstSizeMonomialTransformationsStructure<const N: usize>:
    ConstSizeFiniteSetSignature<N>
{
    fn const_size_monomial_transformations(
        self: &Arc<Self>,
    ) -> Arc<ConstSizeMonomialTransformationsStructure<N, Self>> {
        ConstSizeMonomialTransformationsStructure::new(self.clone())
    }
}
impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>>
    SetToConstSizeMonomialTransformationsStructure<N> for Set
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> Signature
    for ConstSizeMonomialTransformationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> SetSignature
    for ConstSizeMonomialTransformationsStructure<N, Set>
{
    type Elem = ConstSizeMonomialTransformation<N, Set::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        x.validate()
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + EqSignature> EqSignature
    for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        (0..N).all(|i| self.set().equal(&a.scalars[i], &b.scalars[i]))
            && self
                .set()
                .const_size_permutations()
                .equal(&a.permutation, &b.permutation)
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    CompositionSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let s = ConstSizeMonomialTransformation {
            scalars: todo!(),
            permutation: todo!(),
        };
        debug_assert!(self.is_element(&s));
        s
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    AssociativeCompositionSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    LeftCancellativeCompositionSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    RightCancellativeCompositionSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    IdentitySignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        ConstSizeMonomialTransformation {
            scalars: todo!(),
            permutation: self.set().const_size_permutations().identity(),
        }
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    MonoidSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryLeftInverseSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryRightInverseSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryInverseSignature for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature> GroupSignature
    for ConstSizeMonomialTransformationsStructure<N, Set>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        ConstSizeMonomialTransformation {
            scalars: todo!(),
            permutation: self.set().const_size_permutations().identity(),
        }
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_sets::sets::{
        ConstSizeEnumeratedFiniteSet, ConstSizeEnumeratedFiniteSetStructure,
        SetToConstSizePermutationsStructure,
    };
    use algebraeon_structures::{
        MetaCompositionSignature, MetaPermutationsSignature, MetaType, PermutationsSignature,
    };

    use crate::{
        finite_fields::quaternary_field::QuaternaryField,
        linear::const_finitely_free_module::RingToConstFinitelyFreeModuleSignature,
    };

    #[test]
    fn test() {
        let f4 = QuaternaryField::structure();
        let basis = ConstSizeEnumeratedFiniteSetStructure::<6>::new();
        let space = f4.free_module(&basis);

        let b = |i: usize| -> ConstSizeEnumeratedFiniteSet<6> { i.try_into().unwrap() };

        let perm1 = basis
            .const_size_permutations()
            .new_cycle(vec![b(1), b(2), b(3)])
            .unwrap();

        let perm2 = basis
            .const_size_permutations()
            .new_cycle(vec![b(2), b(3), b(4)])
            .unwrap();

        println!("{:?}", perm1.compose(&perm2).disjoint_cycles());
        todo!();
    }
}
