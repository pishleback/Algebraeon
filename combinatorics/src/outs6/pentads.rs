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
pub struct Pentad<Set: EnumeratedOrdFiniteSetSignature> {
    // must have syntheme_1 < syntheme_2 < syntheme_3 < syntheme_4 < syntheme_5 and all disjoint
    syntheme_1: Syntheme<Set>,
    syntheme_2: Syntheme<Set>,
    syntheme_3: Syntheme<Set>,
    syntheme_4: Syntheme<Set>,
    syntheme_5: Syntheme<Set>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PentadsStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for PentadsStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for PentadsStructure<Set, SetB>
{
    type Elem = Pentad<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        todo!()
    }
}
