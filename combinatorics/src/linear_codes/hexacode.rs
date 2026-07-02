use algebraeon_rings::{
    finite_fields::quaternary_field::{QuaternaryField, QuaternaryFieldCanonicalStructure},
    linear::{
        finitely_free_module::{FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature},
        finitely_free_submodule::FinitelyFreeSubmodule,
    },
};
use algebraeon_structures::*;
use std::marker::PhantomData;

/// The hexacode on a 6-element set
///     0 1  2 3  4 5
/// It is the 3-dimensional vector subspace of the free vectorspace over F4 with basis
///     a b  b a  b a
///     b a  a b  b a
///     b a  b a  a b
/// where
///     F4 = {0, 1, a, b}
#[derive(Debug, Clone)]
pub struct HexacodeStructure<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> {
    full_space: FinitelyFreeModuleStructure<
        Set,
        SetB,
        QuaternaryFieldCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
    >,
    hexacode_subspace: FinitelyFreeSubmodule<QuaternaryField>,
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> PartialEq for HexacodeStructure<Set, SetB>
{
    fn eq(&self, other: &Self) -> bool {
        self.full_space == other.full_space
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> Eq for HexacodeStructure<Set, SetB>
{
}

pub trait SetToHexacodeSignature:
    EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>
{
    fn hexacode(&self) -> HexacodeStructure<Self, &Self> {
        HexacodeStructure::new(self)
    }

    fn into_hexacode(self) -> HexacodeStructure<Self, Self> {
        HexacodeStructure::new(self)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>> SetToHexacodeSignature
    for Set
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> HexacodeStructure<Set, SetB>
{
    pub fn new(set: SetB) -> Self {
        let full_space = QuaternaryField::structure().into_free_module(set);
        let hexacode_subspace = full_space.generated_submodule(vec![
            &vec![
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
            ],
            &vec![
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
            ],
            &vec![
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
            ],
        ]);
        Self {
            full_space,
            hexacode_subspace,
        }
    }

    pub fn set(&self) -> &Set {
        self.full_space.set()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> Signature for HexacodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> SetSignature for HexacodeStructure<Set, SetB>
{
    type Elem = Vec<QuaternaryField>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self
            .full_space
            .submodules()
            .contains_element(&self.hexacode_subspace, x)
        {
            return Err("not a hexacodeword".to_string());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::{
        FiniteSetToFinitelySupportedPermutationsStructure, SetToFiniteSubsetByOrdSizedSignature,
    };

    #[test]
    fn test() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        let set_perms = set.permutations();

        let hexacode = set.hexacode();

        println!("{:?}", hexacode);
    }
}
