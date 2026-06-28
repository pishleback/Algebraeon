use crate::outs6::*;
use algebraeon_structures::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct Syntheme<Set: EnumeratedOrdFiniteSetSignature> {
    // must have duad_1 < duad_2 < duad_3 and all disjoint
    duad_1: Duad<Set>,
    duad_2: Duad<Set>,
    duad_3: Duad<Set>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SynthemesStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for SynthemesStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for SynthemesStructure<Set, SetB>
{
    type Elem = Syntheme<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        todo!()
    }
}
