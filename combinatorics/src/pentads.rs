use std::marker::PhantomData;

use algebraeon_structures::*;

#[derive(Debug, Clone, PartialEq, Eq)]
struct PentadsStructure<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for PentadsStructure<Set, SetB>
{
}
