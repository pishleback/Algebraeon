use crate::structure::{
    MultiplicativeGroupSignature, MultiplicativeMonoidSignature,
    MultiplicativeMonoidUnitsSignature, RinglikeSpecializationSignature,
};
use algebraeon_sets::structure::{BorrowedStructure, SetSignature, Signature};
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiplicativeMonoidUnitsStructure<
    M: MultiplicativeMonoidUnitsSignature,
    MB: BorrowedStructure<M>,
> {
    _monoid: PhantomData<M>,
    monoid: MB,
}

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>>
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

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>> Signature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>> SetSignature
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

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>>
    RinglikeSpecializationSignature for MultiplicativeMonoidUnitsStructure<M, MB>
{
}

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>> MultiplicativeMonoidSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn one(&self) -> Self::Set {
        self.monoid().one()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.monoid().mul(a, b)
    }
}

impl<M: MultiplicativeMonoidUnitsSignature, MB: BorrowedStructure<M>> MultiplicativeGroupSignature
    for MultiplicativeMonoidUnitsStructure<M, MB>
{
    fn inv(&self, a: &Self::Set) -> Self::Set {
        self.monoid().try_inv(a).unwrap()
    }
}
