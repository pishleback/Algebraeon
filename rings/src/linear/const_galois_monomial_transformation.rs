//! For a Galois field F a galois monomial transformation of F^n is a monomial transformation and an element of the Galois group of F

use crate::{
    linear::{
        const_finitely_free_module::ConstFinitelyFreeModuleStructure,
        const_monomial_transformation::{
            ConstSizeGaloisActionOnMonomialTransformationsStructure,
            ConstSizeMonomialTransformation, ConstSizeMonomialTransformationsStructure,
        },
    },
    structure::{GaloisFieldWithGroupSignature, TryReciprocalSignature},
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
> ConstSizeGaloisMonomialTransformationsStructure<N, Basis, Field>
{
    pub fn new_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_h(automorphism)
            .into()
    }

    pub fn new_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_n(monomial_transformation)
            .into()
    }

    pub fn new_galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_n_compose_h(monomial_transformation, automorphism)
            .into()
    }

    pub fn new_monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
        monomial_transformation: &ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) -> <Self as SetSignature>::Elem {
        self.semidirect_product_structure()
            .new_h_compose_n(automorphism, monomial_transformation)
            .into()
    }

    pub fn galois_automorphism_part(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> <Field::GaloisGroup as SetSignature>::Elem {
        self.semidirect_product_structure()
            .h_quotient_project(&elem.repr)
    }

    pub fn galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> (
        ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
        <Field::GaloisGroup as SetSignature>::Elem,
    ) {
        self.semidirect_product_structure().n_compose_h(&elem.repr)
    }

    pub fn monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        elem: &<Self as SetSignature>::Elem,
    ) -> (
        <Field::GaloisGroup as SetSignature>::Elem,
        ConstSizeMonomialTransformation<N, Basis::Elem, Field::Elem>,
    ) {
        self.semidirect_product_structure().h_compose_n(&elem.repr)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finite_fields::quaternary_field::QuaternaryField,
        linear::{
            const_finitely_free_module::RingToConstFinitelyFreeModuleSignature,
            monomial_transformations::MonomialTransformationsSignature,
        },
    };
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
