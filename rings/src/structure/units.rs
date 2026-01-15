use crate::structure::TryReciprocalSignature;
use algebraeon_groups::structure::{
    AssociativeCompositionSignature, CommutativeCompositionSignature, CompositionSignature,
    GroupSignature, IdentitySignature, LeftCancellativeCompositionSignature, MonoidSignature,
    RightCancellativeCompositionSignature, TryInverseSignature,
};
use algebraeon_sets::structure::{BorrowedStructure, SetSignature, Signature};
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiplicativeMonoidUnitsStructure<M: TryReciprocalSignature, MB: BorrowedStructure<M>> {
    _monoid: PhantomData<M>,
    monoid: MB,
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>>
    MultiplicativeMonoidUnitsStructure<M, MB>
{
    pub fn new(monoid: MB) -> Self {
        Self {
            _monoid: PhantomData,
            monoid,
        }
    }

    pub fn monoid(&self) -> &M {
        self.monoid.borrow()
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> Signature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> SetSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    type Set = M::Set;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if self.monoid().is_element(x).is_ok() {
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

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> IdentitySignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn identity(&self) -> Self::Set {
        self.monoid().one()
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> CompositionSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.monoid().mul(a, b)
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> AssociativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> CommutativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> TryInverseSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn try_inverse(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.inverse(a))
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> LeftCancellativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn try_left_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> RightCancellativeCompositionSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn try_right_difference(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> MonoidSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: TryReciprocalSignature, MB: BorrowedStructure<M>> GroupSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn inverse(&self, a: &Self::Set) -> Self::Set {
        self.monoid().try_reciprocal(a).unwrap()
    }
}
