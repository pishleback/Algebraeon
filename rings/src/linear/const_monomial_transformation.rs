//! For a ring R, a monomial transformation of R^n is a composition of a permutation of coordinates with a scalar multiplication for each coordinate

use crate::{
    linear::{
        const_finitely_free_module::ConstFinitelyFreeModuleStructure,
        monomial_transformations::MonomialTransformationsSignature,
    },
    structure::{RingSignature, TryReciprocalSignature},
};
use algebraeon_sets::sets::{
    ConstSizePermutation, ConstSizePermutationsStructure, SetToConstSizePermutationsStructure,
};
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ConstSizeMonomialTransformation<const N: usize, BasisElem, RingElem> {
    // The order of operations matters here
    // Represents the operations of:
    //  - Applying the permutation to the coordinates
    //  - Then applying the scalar multiplications component-wise
    permutation: ConstSizePermutation<N, BasisElem>,
    scalars: [RingElem; N],
}

impl<const N: usize, BasisElem: MetaType, RingElem: MetaType> MetaType
    for ConstSizeMonomialTransformation<N, BasisElem, RingElem>
where
    BasisElem::Signature: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    RingElem::Signature: RingSignature + TryReciprocalSignature,
{
    type Signature =
        ConstSizeMonomialTransformationsStructure<N, BasisElem::Signature, RingElem::Signature>;

    fn structure() -> Arc<Self::Signature> {
        ConstSizeMonomialTransformationsStructure::new(
            BasisElem::structure(),
            RingElem::structure(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeMonomialTransformationsStructure<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> {
    basis: Arc<Basis>,
    ring: Arc<Ring>,
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    pub fn new(basis: Arc<Basis>, ring: Arc<Ring>) -> Arc<Self> {
        Self { basis, ring }.into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> ConstFinitelyFreeModuleStructure<N, Basis, Ring>
{
    pub fn monomial_transformations(
        &self,
    ) -> Arc<ConstSizeMonomialTransformationsStructure<N, Basis, Ring>> {
        ConstSizeMonomialTransformationsStructure::new(self.set().clone(), self.ring().clone())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> Signature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> SetSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    type Elem = ConstSizeMonomialTransformation<N, Basis::Elem, Ring::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.basis
            .const_size_permutations()
            .validate_element(&x.permutation)?;
        if !x.scalars.iter().all(|lambda| self.ring.is_unit(lambda)) {
            return Err("Scalars are not all units".to_string());
        }
        Ok(())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature + EqSignature,
    Ring: RingSignature + TryReciprocalSignature + EqSignature,
> EqSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        (0..N).all(|i| self.ring().equal(&a.scalars[i], &b.scalars[i]))
            && self
                .basis()
                .const_size_permutations()
                .equal(&a.permutation, &b.permutation)
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> CompositionSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let s = ConstSizeMonomialTransformation {
            scalars: todo!(),
            permutation: self
                .basis()
                .const_size_permutations()
                .compose(&a.permutation, &b.permutation),
        };
        debug_assert!(self.is_element(&s));
        s
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> AssociativeCompositionSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> LeftCancellativeCompositionSignature
    for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> RightCancellativeCompositionSignature
    for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> IdentitySignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        ConstSizeMonomialTransformation {
            scalars: std::array::from_fn(|_| self.ring().one()),
            permutation: self.basis().const_size_permutations().identity(),
        }
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> MonoidSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryLeftInverseSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryRightInverseSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryInverseSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> GroupSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        ConstSizeMonomialTransformation {
            scalars: todo!(),
            permutation: self
                .basis()
                .const_size_permutations()
                .inverse(&a.permutation),
        }
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> MonomialTransformationsSignature<Basis, Ring>
    for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    type Permutations = ConstSizePermutationsStructure<N, Basis>;

    fn basis(self: &Arc<Self>) -> &Arc<Basis> {
        &self.basis
    }

    fn ring(self: &Arc<Self>) -> &Arc<Ring> {
        &self.ring
    }

    fn new_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem {
        todo!()
    }

    fn new_scalars(self: &Arc<Self>, scalars: Vec<Ring::Elem>) -> Self::Elem {
        todo!()
    }

    fn permutation_part(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> <Self::Permutations as SetSignature>::Elem {
        todo!()
    }

    fn permutation_then_scalars(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (<Self::Permutations as SetSignature>::Elem, Vec<Ring::Elem>) {
        todo!()
    }

    fn scalars_then_permutation(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (Vec<Ring::Elem>, <Self::Permutations as SetSignature>::Elem) {
        todo!()
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
    fn monomial_transformations_group_operations_and_decomposition() {
        let f4 = QuaternaryField::structure();
        let basis = ConstSizeEnumeratedFiniteSetStructure::<4>::new();
        let space = f4.free_module(&basis);
        let mon_trans = space.monomial_transformations();

        let b = |i: usize| -> ConstSizeEnumeratedFiniteSet<4> { i.try_into().unwrap() };

        let perm1 = basis
            .const_size_permutations()
            .new_cycle(vec![b(0), b(1), b(2)])
            .unwrap();

        let perm2 = basis
            .const_size_permutations()
            .new_cycle(vec![b(1), b(2), b(3)])
            .unwrap();

        let scalars1 = [
            QuaternaryField::One,
            QuaternaryField::One,
            QuaternaryField::Alpha,
            QuaternaryField::Alpha,
        ];

        let scalars2 = [
            QuaternaryField::Zero,
            QuaternaryField::Beta,
            QuaternaryField::Alpha,
            QuaternaryField::One,
        ];

        println!("{:?}", perm1.compose(&perm2).disjoint_cycles());
        todo!();
    }
}
