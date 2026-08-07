//! For a ring R, a monomial transformation of R^n is a composition of a permutation of coordinates with a scalar multiplication for each coordinate

use crate::{
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        const_size_functions_to_ring_units_group::{
            ConstSizeFunctionsToRingUnitsGroup,
            ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure,
        },
        monomial_transformations::MonomialTransformationsSignature,
    },
    structure::{GaloisFieldWithGroupSignature, RingSignature, TryReciprocalSignature},
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
    repr:
        SemidirectProductElem<Function<N, BasisElem, RingElem>, ConstSizePermutation<N, BasisElem>>,
}

impl<const N: usize, BasisElem, RingElem>
    From<
        SemidirectProductElem<Function<N, BasisElem, RingElem>, ConstSizePermutation<N, BasisElem>>,
    > for ConstSizeMonomialTransformation<N, BasisElem, RingElem>
{
    fn from(
        repr: SemidirectProductElem<
            Function<N, BasisElem, RingElem>,
            ConstSizePermutation<N, BasisElem>,
        >,
    ) -> Self {
        Self { repr }
    }
}

impl<const N: usize, BasisElem, RingElem> ConstSizeMonomialTransformation<N, BasisElem, RingElem> {
    pub fn scalars(&self) -> &Function<N, BasisElem, RingElem> {
        &self.repr.n
    }

    pub fn permutation(&self) -> &ConstSizePermutation<N, BasisElem> {
        &self.repr.h
    }
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

    fn semidirect_product_structure(
        self: &Arc<Self>,
    ) -> Arc<
        SemidirectProductStructure<
            ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>,
            ConstSizePermutationsStructure<N, Basis>,
            ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure<N, Basis, Ring>,
        >,
    > {
        SemidirectProductStructure::new(
            ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure::new(
                self.basis.clone(),
                self.ring.clone(),
            ),
        )
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
        self.semidirect_product_structure()
            .validate_element(&x.repr)
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature + EqSignature,
    Ring: RingSignature + TryReciprocalSignature + EqSignature,
> EqSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.semidirect_product_structure().equal(&a.repr, &b.repr)
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
        let s = self
            .semidirect_product_structure()
            .compose(&a.repr, &b.repr)
            .into();
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
        self.semidirect_product_structure()
            .try_left_difference(&a.repr, &b.repr)
            .map(|repr| repr.into())
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
        self.semidirect_product_structure()
            .try_right_difference(&a.repr, &b.repr)
            .map(|repr| repr.into())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> IdentitySignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        self.semidirect_product_structure().identity().into()
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
        self.semidirect_product_structure()
            .try_left_inverse(&a.repr)
            .map(|repr| repr.into())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryRightInverseSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.semidirect_product_structure()
            .try_right_inverse(&a.repr)
            .map(|repr| repr.into())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryInverseSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.semidirect_product_structure()
            .try_inverse(&a.repr)
            .map(|repr| repr.into())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> GroupSignature for ConstSizeMonomialTransformationsStructure<N, Basis, Ring>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        let s = self.semidirect_product_structure().inverse(&a.repr).into();
        debug_assert!(self.is_element(&s));
        s
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
        self.semidirect_product_structure()
            .new_h(permutation)
            .into()
    }

    fn new_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem {
        debug_assert!(self.module().is_element(scalars));
        self.semidirect_product_structure().new_n(scalars).into()
    }

    fn new_permutation_then_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem {
        self.semidirect_product_structure()
            .new_n_compose_h(scalars, permutation)
            .into()
    }

    fn new_scalars_then_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem {
        self.semidirect_product_structure()
            .new_h_compose_n(permutation, scalars)
            .into()
    }

    fn permutation_part(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> <Self::Permutations as SetSignature>::Elem {
        debug_assert!(self.is_element(monomial_transformation));
        self.semidirect_product_structure()
            .h_quotient_project(&monomial_transformation.repr)
    }

    fn permutation_then_scalars(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::FinitelyFreeModule as SetSignature>::Elem,
        <Self::Permutations as SetSignature>::Elem,
    ) {
        debug_assert!(self.is_element(monomial_transformation));
        self.semidirect_product_structure()
            .n_compose_h(&monomial_transformation.repr)
    }

    fn scalars_then_permutation(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::Permutations as SetSignature>::Elem,
        <Self::FinitelyFreeModule as SetSignature>::Elem,
    ) {
        debug_assert!(self.is_element(monomial_transformation));
        self.semidirect_product_structure()
            .h_compose_n(&monomial_transformation.repr)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeGaloisActionOnMonomialTransformationsStructure<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> {
    basis: Arc<Basis>,
    field: Arc<Field>,
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> ConstSizeGaloisActionOnMonomialTransformationsStructure<N, Basis, Field>
{
    pub fn new(basis: Arc<Basis>, field: Arc<Field>) -> Arc<Self> {
        Self { basis, field }.into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> Signature for ConstSizeGaloisActionOnMonomialTransformationsStructure<N, Basis, Field>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
>
    LeftGroupActionSignature<
        Field::GaloisGroup,
        ConstSizeMonomialTransformationsStructure<N, Basis, Field>,
    > for ConstSizeGaloisActionOnMonomialTransformationsStructure<N, Basis, Field>
{
    fn group(self: &Arc<Self>) -> Arc<Field::GaloisGroup> {
        self.field.clone().galois_group_action().group()
    }

    fn set(self: &Arc<Self>) -> Arc<ConstSizeMonomialTransformationsStructure<N, Basis, Field>> {
        self.field
            .free_module(&self.basis)
            .monomial_transformations()
    }

    fn apply(
        self: &Arc<Self>,
        g: &<Field::GaloisGroup as SetSignature>::Elem,
        x: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) -> ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem> {
        let fns = self.basis.const_size_functions_to(&self.field);
        self.set().new_scalars_then_permutation(
            x.permutation(),
            &fns.function(|i| {
                self.field
                    .clone()
                    .galois_group_action()
                    .apply(g, fns.image(x.scalars(), i))
            })
            .unwrap(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_fields::quaternary_field::QuaternaryField;
    use algebraeon_sets::sets::{
        ConstSizeEnumeratedFiniteSet, ConstSizeEnumeratedFiniteSetStructure,
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
            &[F4::Alpha, F4::Alpha, F4::Beta].into(),
            &basis
                .const_size_permutations()
                .new_cycle(vec![b(0), b(1), b(2)])
                .unwrap(),
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
