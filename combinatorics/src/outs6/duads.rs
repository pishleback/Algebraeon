//! Given a 6-element set make the following definitions
//!  - A duad is a 2-element subset
//!  - A syntheme is an unordered set of 3 disjoint duads
//!  - A pentad is an unordered set of 5 disjoint synthemes
//!
//! Then there are 6 pentads, and any permutation of the 6 points induces a permutation of the 6 pentads via an outer automorphism of S6

use crate::outs6::*;
use algebraeon_structures::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct Duad<Set: EnumeratedOrdFiniteSetSignature> {
    // must have point_1 < point_2
    point_1: Set::Elem,
    point_2: Set::Elem,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuadsStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for DuadsStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for DuadsStructure<Set, SetB>
{
    type Elem = Duad<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        todo!()
    }
}
