use algebraeon_rings::{
    finite_fields::quaternary_field::{QuaternaryField, QuaternaryFieldCanonicalStructure},
    linear::{
        finitely_free_module::{FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature},
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
    },
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, RinglikeSpecializationSignature, SemiModuleSignature,
        TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;

/// The hexacode on a 6-element set
///     0 1  2 3  4 5
/// It is the 3-dimensional vector subspace of the free vectorspace over F4 with basis
///     a b  b a  b a
///     b a  a b  b a
///     b a  b a  a b
/// where
///     F4 = {0, 1, a, b}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HexacodeStructure<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> {
    hexacode: FinitelyFreeSubmoduleStructure<
        Set,
        SetB,
        QuaternaryFieldCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
        FinitelyFreeModuleStructure<
            Set,
            SetB,
            QuaternaryFieldCanonicalStructure,
            QuaternaryFieldCanonicalStructure,
        >,
    >,
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
            hexacode: FinitelyFreeSubmoduleStructure::new(full_space, hexacode_subspace),
        }
    }

    pub fn set(&self) -> &Set {
        self.hexacode.set()
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
        self.hexacode.validate_element(x)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> RinglikeSpecializationSignature for HexacodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> ZeroSignature for HexacodeStructure<Set, SetB>
{
    fn zero(&self) -> Self::Elem {
        self.hexacode.zero()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> AdditionSignature for HexacodeStructure<Set, SetB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.hexacode.add(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> CancellativeAdditionSignature for HexacodeStructure<Set, SetB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.hexacode.try_sub(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> TryNegateSignature for HexacodeStructure<Set, SetB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.hexacode.try_neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> AdditiveMonoidSignature for HexacodeStructure<Set, SetB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> AdditiveGroupSignature for HexacodeStructure<Set, SetB>
{
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.hexacode.neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> SemiModuleSignature<QuaternaryFieldCanonicalStructure> for HexacodeStructure<Set, SetB>
{
    fn ring(&self) -> &QuaternaryFieldCanonicalStructure {
        self.hexacode.ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &QuaternaryField) -> Self::Elem {
        self.hexacode.scalar_mul(a, x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::SetToFiniteSubsetByOrdSizedSignature;

    #[test]
    fn test() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        // let set_perms = set.permutations();

        let hexacode = set.hexacode();

        println!("{:?}", hexacode);
    }
}
