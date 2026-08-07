//! For a ring `R` a product of `n` coppies of the unit group `R^*` as a subset of `R^n`.

use crate::structure::{RingSignature, TryReciprocalSignature};
use algebraeon_sets::sets::{
    ConstSizePermutation, ConstSizePermutationsStructure, Function,
    SetToConstSizeFunctionsToSignature, SetToConstSizePermutationsStructure,
};
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeFunctionsToRingUnitsGroup<
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
> ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    pub fn new(basis: Arc<Basis>, ring: Arc<Ring>) -> Arc<Self> {
        Self { basis, ring }.into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> Signature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> SetSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    type Elem = Function<N, Basis::Elem, Ring::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !x.iter().all(|lambda| self.ring.is_unit(lambda)) {
            return Err("Scalars are not all units".to_string());
        }
        Ok(())
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature + EqSignature,
    Ring: RingSignature + TryReciprocalSignature + EqSignature,
> EqSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        (0..N).all(|i| self.ring.equal(&a[i], &b[i]))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> CompositionSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let fns = self.basis.const_size_functions_to(&self.ring);
        let s = fns
            .function(|i| self.ring.mul(fns.image(a, i), fns.image(b, i)))
            .unwrap();
        debug_assert!(self.is_element(&s));
        s
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> AssociativeCompositionSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> LeftCancellativeCompositionSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> RightCancellativeCompositionSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
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
> IdentitySignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        self.basis
            .const_size_functions_to(&self.ring)
            .function(|_| self.ring.one())
            .unwrap()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> MonoidSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryLeftInverseSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryRightInverseSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> TryInverseSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> GroupSignature for ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        let fns = self.basis.const_size_functions_to(&self.ring);
        let s = fns
            .function(|i| self.ring.units_group().inverse(fns.image(a, i)))
            .unwrap();
        debug_assert!(self.is_element(&s));
        s
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure<
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
> ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure<N, Basis, Ring>
{
    pub fn new(basis: Arc<Basis>, ring: Arc<Ring>) -> Arc<Self> {
        Self { basis, ring }.into()
    }
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
> Signature for ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure<N, Basis, Ring>
{
}

impl<
    const N: usize,
    Basis: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Ring: RingSignature + TryReciprocalSignature,
>
    LeftGroupActionSignature<
        ConstSizePermutationsStructure<N, Basis>,
        ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>,
    > for ConstSizeLeftPermutationActionOnFunctionsToRingUnitsGroupStructure<N, Basis, Ring>
{
    fn group(self: &Arc<Self>) -> Arc<ConstSizePermutationsStructure<N, Basis>> {
        self.basis.const_size_permutations()
    }

    fn set(self: &Arc<Self>) -> Arc<ConstSizeFunctionsToRingUnitsGroup<N, Basis, Ring>> {
        ConstSizeFunctionsToRingUnitsGroup::new(self.basis.clone(), self.ring.clone())
    }

    fn apply(
        self: &Arc<Self>,
        g: &ConstSizePermutation<N, Basis::Elem>,
        x: &Function<N, Basis::Elem, Ring::Elem>,
    ) -> Function<N, Basis::Elem, Ring::Elem> {
        self.basis
            .const_size_functions_to(&self.ring)
            .output_const_size_permutation_action()
            .apply(g, x)
    }
}
