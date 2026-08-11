use crate::structure::{
    FinitelyFreeModuleSignature, GaloisFieldWithGroupSignature, RingSignature,
    TryReciprocalSignature,
};
use algebraeon_macros::{signature_meta_trait, skip_meta};
use algebraeon_structures::*;
use std::sync::Arc;

#[signature_meta_trait]
pub trait MonomialTransformationsSupersetSignature<
    Basis: OrderedFiniteSetSignature,
    Ring: RingSignature,
>: GroupSignature
{
    type Permutations: PermutationsSignature<Basis>;
    type FinitelyFreeModule: FinitelyFreeModuleSignature<Basis, Ring>;

    #[skip_meta]
    fn basis(self: &Arc<Self>) -> Arc<Basis>;

    #[skip_meta]
    fn basis_permutations(self: &Arc<Self>) -> Arc<Self::Permutations>;

    #[skip_meta]
    fn ring(self: &Arc<Self>) -> Arc<Ring>;

    #[skip_meta]
    fn module(self: &Arc<Self>) -> Arc<Self::FinitelyFreeModule>;

    /// Construct a monomial transformation which acts purely as a permutation
    fn new_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem;

    /// Construt a monomial transformation which acts purely by scalar multiplications
    fn new_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem;

    /// Construct a monomial transformation which acts as a permutation followed by scalar multiplications
    fn new_permutation_then_scalars(
        self: &Arc<Self>,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
        permutation: &<Self::Permutations as SetSignature>::Elem,
    ) -> Self::Elem {
        debug_assert!(self.module().is_element(scalars));
        debug_assert!(self.basis_permutations().is_element(permutation));
        self.compose(
            &self.new_scalars(scalars),
            &self.new_permutation(permutation),
        )
    }

    /// Construct a monomial transformation which acts as scalar multiplications followed by a permutation
    fn new_scalars_then_permutation(
        self: &Arc<Self>,
        permutation: &<Self::Permutations as SetSignature>::Elem,
        scalars: &<Self::FinitelyFreeModule as SetSignature>::Elem,
    ) -> Self::Elem {
        debug_assert!(self.module().is_element(scalars));
        debug_assert!(self.basis_permutations().is_element(permutation));
        self.compose(
            &self.new_permutation(permutation),
            &self.new_scalars(scalars),
        )
    }
}

#[signature_meta_trait]
pub trait MonomialTransformationsSignature<Basis: OrderedFiniteSetSignature, Ring: RingSignature>:
    MonomialTransformationsSupersetSignature<Basis, Ring>
{
    /// Get the permutation part of a monomial transformation
    fn permutation_part(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> <<Self as MonomialTransformationsSupersetSignature<Basis, Ring>>::Permutations as SetSignature>::Elem;

    /// Decompose a monomial transformation as a permutation followed by scalar multiplications
    fn permutation_then_scalars(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <<Self as MonomialTransformationsSupersetSignature<Basis, Ring>>::FinitelyFreeModule as SetSignature>::Elem,
        <<Self as MonomialTransformationsSupersetSignature<Basis, Ring>>::Permutations as SetSignature>::Elem,
    );

    /// Decompose a monomial transformation as scalar multiplications followed by a permutation
    fn scalars_then_permutation(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <<Self as MonomialTransformationsSupersetSignature<Basis, Ring>>::Permutations as SetSignature>::Elem,
        <<Self as MonomialTransformationsSupersetSignature<Basis, Ring>>::FinitelyFreeModule as SetSignature>::Elem,
    );
}

#[signature_meta_trait]
pub trait GaloisMonomialTransformationsSignature<
    Basis: OrderedFiniteSetSignature,
    Field: GaloisFieldWithGroupSignature + TryReciprocalSignature,
>: MonomialTransformationsSupersetSignature<Basis, Field>
{
    type MonomialTransformations: MonomialTransformationsSignature<Basis, Field>;

    fn new_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem;

    fn new_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &<Self::MonomialTransformations as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem;

    fn new_galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        monomial_transformation: &<Self::MonomialTransformations as SetSignature>::Elem,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.compose(
            &self.new_monomial_transformation(monomial_transformation),
            &self.new_galois_automorphism(automorphism),
        )
    }

    fn new_monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        automorphism: &<Field::GaloisGroup as SetSignature>::Elem,
        monomial_transformation: &<Self::MonomialTransformations as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        self.compose(
            &self.new_galois_automorphism(automorphism),
            &self.new_monomial_transformation(monomial_transformation),
        )
    }

    fn galois_automorphism_part(
        self: &Arc<Self>,
        elem: &Self::Elem,
    ) -> <Field::GaloisGroup as SetSignature>::Elem;

    fn galois_automorphism_then_monomial_transformation(
        self: &Arc<Self>,
        elem: &Self::Elem,
    ) -> (
        <Self::MonomialTransformations as SetSignature>::Elem,
        <Field::GaloisGroup as SetSignature>::Elem,
    );

    fn monomial_transformation_then_galois_automorphism(
        self: &Arc<Self>,
        elem: &Self::Elem,
    ) -> (
        <Field::GaloisGroup as SetSignature>::Elem,
        <Self::MonomialTransformations as SetSignature>::Elem,
    );
}
