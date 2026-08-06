//! For a ring R, a monomial transformation of R^n is a composition of a permutation of coordinates with a scalar multiplication for each coordinate

use crate::{
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        monomial_transformations::MonomialTransformationsSignature,
    },
    structure::{RingSignature, TryReciprocalSignature},
};
use algebraeon_sets::sets::{
    ConstSizePermutation, ConstSizePermutationsStructure, Function,
    SetToConstSizeFunctionsToSignature, SetToConstSizePermutationsStructure,
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
    scalars: Function<N, BasisElem, RingElem>,
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

    fn identity_scalars(self: &Arc<Self>) -> Function<N, Basis::Elem, Ring::Elem> {
        self.basis()
            .const_size_functions_to(self.ring())
            .function(|_| self.ring().one())
            .unwrap()
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
        let mod_fns = self.module().functions_restructure();
        let s = ConstSizeMonomialTransformation {
            scalars: mod_fns
                .function(|i| {
                    self.ring().mul(
                        mod_fns.image(&a.scalars, i),
                        mod_fns.image(
                            &b.scalars,
                            &self
                                .basis()
                                .const_size_permutations()
                                .preimage(&a.permutation, i),
                        ),
                    )
                })
                .unwrap(),
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
            scalars: self.identity_scalars(),
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
        let mod_fns = self.module().functions_restructure();
        ConstSizeMonomialTransformation {
            scalars: self
                .module()
                .functions_restructure()
                .output_const_size_permutation_action()
                .apply_inverse(
                    &a.permutation,
                    &mod_fns
                        .function(|i| {
                            self.ring()
                                .try_reciprocal(mod_fns.image(&a.scalars, i))
                                .unwrap()
                        })
                        .unwrap(),
                ),
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
    type FinitelyFreeModule = ConstFinitelyFreeModuleStructure<N, Basis, Ring>;

    fn basis(self: &Arc<Self>) -> &Arc<Basis> {
        &self.basis
    }

    fn basis_permutations(self: &Arc<Self>) -> Arc<Self::Permutations> {
        self.basis().const_size_permutations()
    }

    fn ring(self: &Arc<Self>) -> &Arc<Ring> {
        &self.ring
    }

    fn module(self: &Arc<Self>) -> Arc<Self::FinitelyFreeModule> {
        self.ring().free_module(self.basis())
    }

    fn new_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem {
        debug_assert!(self.basis_permutations().is_element(permutation));
        ConstSizeMonomialTransformation {
            permutation: permutation.clone(),
            scalars: self.identity_scalars(),
        }
    }

    fn new_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem {
        debug_assert!(self.module().is_element(scalars));
        ConstSizeMonomialTransformation {
            permutation: self.basis_permutations().identity(),
            scalars: scalars.clone(),
        }
    }

    fn permutation_part(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> <Self::Permutations as SetSignature>::Elem {
        debug_assert!(self.is_element(monomial_transformation));
        monomial_transformation.permutation.clone()
    }

    fn permutation_then_scalars(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::FinitelyFreeModule as SetSignature>::Elem,
        <Self::Permutations as SetSignature>::Elem,
    ) {
        debug_assert!(self.is_element(monomial_transformation));
        (
            monomial_transformation.scalars.clone(),
            monomial_transformation.permutation.clone(),
        )
    }

    fn scalars_then_permutation(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::Permutations as SetSignature>::Elem,
        <Self::FinitelyFreeModule as SetSignature>::Elem,
    ) {
        debug_assert!(self.is_element(monomial_transformation));
        (
            monomial_transformation.permutation.clone(),
            self.module()
                .functions_restructure()
                .output_const_size_permutation_action()
                .apply_inverse(
                    &monomial_transformation.permutation,
                    &monomial_transformation.scalars,
                ),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        finite_fields::quaternary_field::QuaternaryField,
        linear::{
            const_finitely_free_module::RingToConstFinitelyFreeModuleSignature,
            monomial_transformations::MonomialTransformationsSignature,
        },
    };
    use algebraeon_sets::sets::{
        ConstSizeEnumeratedFiniteSet, ConstSizeEnumeratedFiniteSetStructure,
        SetToConstSizePermutationsStructure,
    };
    use algebraeon_structures::{
        CompositionSignature, EqSignature, GroupSignature, IdentitySignature,
        LeftGroupActionSignature, MetaType, PermutationsSignature,
    };

    #[test]
    fn monomial_transformations_group_operations_and_decomposition() {
        let f4 = QuaternaryField::structure();
        let basis = ConstSizeEnumeratedFiniteSetStructure::<3>::new();
        let space = f4.free_module(&basis);
        let mon_trans = space.monomial_transformations();

        let b = |i: usize| -> ConstSizeEnumeratedFiniteSet<3> { i.try_into().unwrap() };

        let perm = mon_trans.new_permutation(
            &basis
                .const_size_permutations()
                .new_cycle(vec![b(0), b(1), b(2)])
                .unwrap(),
        );

        let scalars0 = mon_trans.new_scalars(
            &[
                QuaternaryField::One,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
            ]
            .into(),
        );

        let scalars1 = mon_trans.new_scalars(
            &[
                QuaternaryField::Beta,
                QuaternaryField::One,
                QuaternaryField::Alpha,
            ]
            .into(),
        );

        let scalars2 = mon_trans.new_scalars(
            &[
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::One,
            ]
            .into(),
        );

        for scalars in [&scalars0, &scalars1, &scalars2] {
            assert!(
                mon_trans.equal(
                    &mon_trans.new_scalars(
                        &mon_trans
                            .permutation_then_scalars(&mon_trans.compose(scalars, &perm))
                            .0
                    ),
                    scalars
                )
            );
        }

        for scalars in [&scalars0, &scalars1, &scalars2] {
            assert!(
                mon_trans.equal(
                    &mon_trans.new_scalars(
                        &mon_trans
                            .scalars_then_permutation(&mon_trans.compose(&perm, scalars))
                            .1
                    ),
                    scalars
                )
            );
        }

        assert!(
            mon_trans.equal(
                &mon_trans.new_scalars(
                    &mon_trans
                        .permutation_then_scalars(&mon_trans.compose(&perm, &scalars0))
                        .0
                ),
                &scalars1
            )
        );

        for scalars in [&scalars0, &scalars1, &scalars2] {
            let t = mon_trans.compose(scalars, &perm);
            let t_inv = mon_trans.inverse(&t);
            assert!(mon_trans.equal(&mon_trans.identity(), &mon_trans.compose(&t, &t_inv)));
        }
    }

    #[test]
    fn monomial_transformation_application() {
        type F4 = QuaternaryField;

        let f4 = F4::structure();
        let basis = ConstSizeEnumeratedFiniteSetStructure::<3>::new();
        let space = f4.free_module(&basis);
        let mon_trans = space.monomial_transformations();

        let b = |i: usize| -> ConstSizeEnumeratedFiniteSet<3> { i.try_into().unwrap() };

        let trans = mon_trans.new_permutation_then_scalars(
            &basis
                .const_size_permutations()
                .new_cycle(vec![b(0), b(1), b(2)])
                .unwrap(),
            &[F4::Alpha, F4::Alpha, F4::Beta].into(),
        );

        let mon_trans_action = space.monomial_transformation_action();

        assert!(space.equal(
            &mon_trans_action.apply(&trans, &[F4::One, F4::Zero, F4::Zero].into()),
            &[F4::Zero, F4::Alpha, F4::Zero].into(),
        ));

        assert!(space.equal(
            &mon_trans_action.apply(&trans, &[F4::Zero, F4::One, F4::Zero].into()),
            &[F4::Zero, F4::Zero, F4::Beta].into(),
        ));

        assert!(space.equal(
            &mon_trans_action.apply(&trans, &[F4::Zero, F4::Zero, F4::One].into()),
            &[F4::Alpha, F4::Zero, F4::Zero].into(),
        ));
    }
}
