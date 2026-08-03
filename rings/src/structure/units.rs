use crate::structure::TryReciprocalSignature;
use algebraeon_structures::*;
use std::sync::Arc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiplicativeMonoidUnitsStructure<M: TryReciprocalSignature> {
    monoid: Arc<M>,
}

impl<M: TryReciprocalSignature> MultiplicativeMonoidUnitsStructure<M> {
    pub fn new(monoid: Arc<M>) -> Arc<Self> {
        Self { monoid }.into()
    }

    pub fn monoid(&self) -> &Arc<M> {
        &self.monoid
    }
}

impl<M: TryReciprocalSignature> Signature for MultiplicativeMonoidUnitsStructure<M> {}

impl<M: TryReciprocalSignature> SetSignature for MultiplicativeMonoidUnitsStructure<M> {
    type Elem = M::Elem;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if self.monoid().validate_element(x).is_ok() {
            if self.monoid().is_unit(x) {
                Ok(())
            } else {
                Err("not a unit".to_string())
            }
        } else {
            Err("not an element of the monoid".to_string())
        }
    }
}

impl<M: TryReciprocalSignature> IdentitySignature for MultiplicativeMonoidUnitsStructure<M> {
    fn identity(self: &Arc<Self>) -> Self::Elem {
        self.monoid().one()
    }
}

impl<M: TryReciprocalSignature> CompositionSignature for MultiplicativeMonoidUnitsStructure<M> {
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.monoid().mul(a, b)
    }
}

impl<M: TryReciprocalSignature> AssociativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M>
{
}

impl<M: TryReciprocalSignature> CommutativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M>
{
}

impl<M: TryReciprocalSignature> TryInverseSignature for MultiplicativeMonoidUnitsStructure<M> {
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<M: TryReciprocalSignature> LeftCancellativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<M: TryReciprocalSignature> RightCancellativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<M: TryReciprocalSignature> MonoidSignature for MultiplicativeMonoidUnitsStructure<M> {}

impl<M: TryReciprocalSignature> GroupSignature for MultiplicativeMonoidUnitsStructure<M> {
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        self.monoid().try_reciprocal(a).unwrap()
    }
}
