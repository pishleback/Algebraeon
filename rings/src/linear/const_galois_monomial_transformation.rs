//! For a Galois field F a galois monomial transformation of F^n is a monomial transformation and an element of the Galois group of F

use crate::{
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        const_monomial_transformation::{
            ConstSizeGaloisActionOnMonomialTransformationsStructure,
            ConstSizeMonomialTransformation, ConstSizeMonomialTransformationsStructure,
        },
        monomial_transformations::{
            GaloisMonomialTransformationsSignature, MonomialTransformationsSupersetSignature,
        },
    },
    structure::{GaloisFieldWithGroupSignature, TryReciprocalSignature},
};
use algebraeon_sets::sets::{
    ConstSizePermutationsStructure, Function, SetToConstSizePermutationsStructure,
};
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ConstSizeGaloisMonomialTransformation<const N: usize, BasisElem, FieldElem, GaloisElem> {
    // The order of operations matters here
    // Represents the operations of:
    //  - Applying the Galois group element to all components
    //  - Then applying the monomial transformation
    repr:
        SemidirectProductElem<ConstSizeMonomialTransformation<N, BasisElem, FieldElem>, GaloisElem>,
}

impl<const N: usize, BasisElem, FieldElem, GaloisElem>
    From<
        SemidirectProductElem<ConstSizeMonomialTransformation<N, BasisElem, FieldElem>, GaloisElem>,
    > for ConstSizeGaloisMonomialTransformation<N, BasisElem, FieldElem, GaloisElem>
{
    fn from(
        repr: SemidirectProductElem<
            ConstSizeMonomialTransformation<N, BasisElem, FieldElem>,
            GaloisElem,
        >,
    ) -> Self {
        Self { repr }
    }
}

impl<const N: usize, BasisElem: MetaType, FieldElem: MetaType> MetaType
    for ConstSizeGaloisMonomialTransformation<N, BasisElem, FieldElem,
     <<FieldElem::Signature as GaloisFieldWithGroupSignature>::GaloisGroup as SetSignature>::Elem>
where
    BasisElem::Signature: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    FieldElem::Signature: GaloisFieldWithGroupSignature + TryReciprocalSignature,
{
    type Signature = ConstSizeGaloisMonomialTransformationsStructure<
        N,
        BasisElem::Signature,
        FieldElem::Signature,
    >;

    fn structure() -> Arc<Self::Signature> {
        ConstSizeGaloisMonomialTransformationsStructure::new(
            BasisElem::structure(),
            FieldElem::structure(),
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeGaloisMonomialTransformationsStructure<
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
> ConstFinitelyFreeModuleStructure<N, Basis, Field>
{
    pub fn galois_monomial_transformations(
        &self,
    ) -> Arc<ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>> {
        ConstSizeGaloisMonomialTransformationsStructure::new(
            self.set().clone(),
            self.ring().clone(),
        )
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    pub fn new(basis: Arc<Basis>, field: Arc<Field>) -> Arc<Self> {
        Self { basis, field }.into()
    }

    fn semidirect_product_structure(
        self: &Arc<Self>,
    ) -> Arc<
        SemidirectProductStructure<
            ConstSizeMonomialTransformationsStructure<N, Basis, Field>,
            Field::GaloisGroup,
            ConstSizeGaloisActionOnMonomialTransformationsStructure<N, Basis, Field>,
        >,
    > {
        SemidirectProductStructure::new(
            ConstSizeGaloisActionOnMonomialTransformationsStructure::new(
                self.basis.clone(),
                self.field.clone(),
            ),
        )
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> Signature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> SetSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    type Elem = ConstSizeGaloisMonomialTransformation<
        N,
        Basis::Elem,
        Field::Elem,
        <Field::GaloisGroup as SetSignature>::Elem,
    >;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.semidirect_product_structure()
            .validate_element(&x.repr)
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> EqSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
where
    Field::GaloisGroup: EqSignature,
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.semidirect_product_structure().equal(&a.repr, &b.repr)
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> CompositionSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> AssociativeCompositionSignature
    for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> LeftCancellativeCompositionSignature
    for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> RightCancellativeCompositionSignature
    for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> IdentitySignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        self.semidirect_product_structure().identity().into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> MonoidSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> TryLeftInverseSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> TryRightInverseSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> TryInverseSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> GroupSignature for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
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
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> MonomialTransformationsSupersetSignature<Basis, Field>
    for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    type Permutations = ConstSizePermutationsStructure<N, Basis>;
    type FinitelyFreeModule = ConstFinitelyFreeModuleStructure<N, Basis, Field>;

    fn basis(self: &Arc<Self>) -> Arc<Basis> {
        self.basis.clone()
    }

    fn basis_permutations(self: &Arc<Self>) -> Arc<Self::Permutations> {
        self.basis().const_size_permutations()
    }

    fn ring(self: &Arc<Self>) -> Arc<Field> {
        self.field.clone()
    }

    fn module(self: &Arc<Self>) -> Arc<Self::FinitelyFreeModule> {
        self.field.free_module(&self.basis)
    }

    fn new_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem {
        self.new_monomial_transformation(
            &self
                .semidirect_product_structure()
                .group_n()
                .new_permutation(permutation),
        )
    }

    fn new_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem {
        self.new_monomial_transformation(
            &self
                .semidirect_product_structure()
                .group_n()
                .new_scalars(scalars),
        )
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> GaloisMonomialTransformationsSignature<Basis, Field>
    for ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    type MonomialTransformations = ConstSizeMonomialTransformationsStructure<N, Basis, Field>;

    fn new_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_h(automorphism)
            .into()
    }

    fn new_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_n(monomial_transformation)
            .into()
    }

    fn new_galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_n_compose_h(monomial_transformation, automorphism)
            .into()
    }

    fn new_monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_h_compose_n(automorphism, monomial_transformation)
            .into()
    }

    fn galois_automorphism_part(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> <Field::GaloisGroup as SetSignature>::Elem {
        self.semidirect_product_structure()
            .h_quotient_project(&elem.repr)
    }

    fn galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> (
        ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
        <Field::GaloisGroup as SetSignature>::Elem,
    ) {
        self.semidirect_product_structure().n_compose_h(&elem.repr)
    }

    fn monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> (
        <Field::GaloisGroup as SetSignature>::Elem,
        ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) {
        self.semidirect_product_structure().h_compose_n(&elem.repr)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LeftGaloisMonomialTransformationActionOnConstFinitelyFreeModuleStructure<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> {
    module: Arc<ConstFinitelyFreeModuleStructure<N, Basis, Field>>,
    group: Arc<ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>>,
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> LeftGaloisMonomialTransformationActionOnConstFinitelyFreeModuleStructure<N, Basis, Field>
{
    fn new(
        module: Arc<ConstFinitelyFreeModuleStructure<N, Basis, Field>>,
        group: Arc<ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>>,
    ) -> Arc<Self> {
        Self { module, group }.into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> Signature
    for LeftGaloisMonomialTransformationActionOnConstFinitelyFreeModuleStructure<N, Basis, Field>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
> ConstFinitelyFreeModuleStructure<N, Basis, Field>
{
    pub fn galois_monomial_transformation_action(
        self: &Arc<Self>,
    ) -> Arc<
        impl LeftGroupActionSignature<
            ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>,
            Self,
        >,
    > {
        LeftGaloisMonomialTransformationActionOnConstFinitelyFreeModuleStructure::new(
            self.clone(),
            self.galois_monomial_transformations(),
        )
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
>
    LeftGroupActionSignature<
        ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>,
        ConstFinitelyFreeModuleStructure<N, Basis, Field>,
    >
    for LeftGaloisMonomialTransformationActionOnConstFinitelyFreeModuleStructure<N, Basis, Field>
{
    fn group(
        self: &Arc<Self>,
    ) -> Arc<ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>> {
        self.group.clone()
    }

    fn set(self: &Arc<Self>) -> Arc<ConstFinitelyFreeModuleStructure<N, Basis, Field>> {
        self.module.clone()
    }

    fn apply(
        self: &Arc<Self>,
        g: &ConstSizeGaloisMonomialTransformation<
            N,
            Basis::Elem,
            Field::Elem,
            <Field::GaloisGroup as SetSignature>::Elem,
        >,
        vec: &Function<N, Basis::Elem, Field::Elem>,
    ) -> Function<N, Basis::Elem, Field::Elem> {
        let mod_fns = self.module.functions_restructure();
        let field = self.module.ring();
        let (monomial, galois_aut) = self
            .group
            .galois_automorphism_then_monomial_transformation(g);
        self.module.monomial_transformation_action().apply(
            &monomial,
            &mod_fns
                .function(|i| {
                    field
                        .clone()
                        .galois_group_action()
                        .apply(&galois_aut, mod_fns.image(vec, i))
                })
                .unwrap(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::finite_fields::quaternary_field::QuaternaryField;
    use algebraeon_groups::examples::c2::C2;
    use algebraeon_sets::sets::ConstSizeEnumeratedFiniteSetStructure;

    #[test]
    fn monomial_transformations_group_operations_and_decomposition() {
        type F4 = QuaternaryField;
        let f4 = QuaternaryField::structure();
        let basis = ConstSizeEnumeratedFiniteSetStructure::<3>::new();
        let space = f4.free_module(&basis);
        let mon_trans = space.monomial_transformations();
        let gal_mon_trans = space.galois_monomial_transformations();

        let mon1 = mon_trans.new_scalars(&[F4::Alpha, F4::Alpha, F4::Alpha].into());
        let mon2 = mon_trans.new_scalars(&[F4::Beta, F4::Beta, F4::Beta].into());

        assert!(gal_mon_trans.equal(
            &gal_mon_trans.new_monomial_transformation_then_galois_automorphism(&C2::Flip, &mon1),
            &gal_mon_trans.new_galois_automorphism_then_monomial_transformation(&mon2, &C2::Flip),
        ));
    }
}
