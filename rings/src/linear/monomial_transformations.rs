use crate::structure::{FinitelyFreeModuleSignature, RingSignature};
use algebraeon_macros::{signature_meta_trait, skip_meta};
use algebraeon_structures::*;
use std::sync::Arc;

#[signature_meta_trait]
pub trait MonomialTransformationsSignature<Basis: OrderedFiniteSetSignature, Ring: RingSignature>:
    GroupSignature
{
    type Permutations: PermutationsSignature<Basis>;
    type FinitelyFreeModule: FinitelyFreeModuleSignature<Basis, Ring>;

    #[skip_meta]
    fn basis(self: &Arc<Self>) -> &Arc<Basis>;

    #[skip_meta]
    fn basis_permutations(self: &Arc<Self>) -> Arc<Self::Permutations>;

    #[skip_meta]
    fn ring(self: &Arc<Self>) -> &Arc<Ring>;

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

    /// Get the permutation part of a monomial transformation
    fn permutation_part(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> <Self::Permutations as SetSignature>::Elem;

    /// Decompose a monomial transformation as a permutation followed by scalar multiplications
    fn permutation_then_scalars(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::FinitelyFreeModule as SetSignature>::Elem,
        <Self::Permutations as SetSignature>::Elem,
    );

    /// Decompose a monomial transformation as scalar multiplications followed by a permutation
    fn scalars_then_permutation(
        self: &Arc<Self>,
        monomial_transformation: &Self::Elem,
    ) -> (
        <Self::Permutations as SetSignature>::Elem,
        <Self::FinitelyFreeModule as SetSignature>::Elem,
    );
}
