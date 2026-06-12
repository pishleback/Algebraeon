use algebraeon_structures::*;
use std::marker::PhantomData;

pub struct Duad<Elem> {
    // a and b must be distinct
    pub a: Elem,
    pub b: Elem,
}

pub struct Syntheme<Elem> {
    // the three duads must be disjoint
    // their order does not matter
    pub duads: [Duad<Elem>; 3],
}

pub struct Pentad<Elem> {
    pub synthemes: [Syntheme<Elem>; 5],
}

/// The set of pentads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
struct PentadsStructure<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> PentadsStructure<Set, SetB> {
    pub fn new(set: SetB) -> Option<Self> {
        if set.borrow().size() == 6 {
            Some(Self {
                _set: PhantomData,
                set,
            })
        } else {
            None
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: FiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for PentadsStructure<Set, SetB>
{
}
