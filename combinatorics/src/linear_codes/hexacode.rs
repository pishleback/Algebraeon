use std::cmp::Ordering;

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

// not necessarily a codeword
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct HexacodeVector {
    // a length 6 vector
    vec: Vec<QuaternaryField>,
}

impl TryFrom<Vec<QuaternaryField>> for HexacodeVector {
    type Error = ();

    fn try_from(vec: Vec<QuaternaryField>) -> Result<Self, Self::Error> {
        if vec.len() == 6 {
            Ok(Self { vec })
        } else {
            Err(())
        }
    }
}

impl From<HexacodeVector> for Vec<QuaternaryField> {
    fn from(val: HexacodeVector) -> Self {
        val.vec
    }
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
    type Elem = HexacodeVector;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.hexacode.validate_element(&x.vec)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> EqSignature for HexacodeStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.hexacode.equal(&a.vec, &b.vec)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> PartialOrdSignature for HexacodeStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.hexacode.partial_cmp(&a.vec, &b.vec)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> OrdSignature for HexacodeStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.hexacode.cmp(&a.vec, &b.vec)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> CountableSetSignature for HexacodeStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.hexacode
            .into_generate_all_elements()
            .map(|vec| HexacodeVector { vec })
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.hexacode
            .generate_all_elements()
            .map(|vec| HexacodeVector { vec })
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> FiniteSetSignature for HexacodeStructure<Set, SetB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.hexacode
            .list_all_elements()
            .into_iter()
            .map(|vec| HexacodeVector { vec })
            .collect()
    }

    fn size(&self) -> Natural {
        self.hexacode.size()
    }

    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        self.hexacode
            .generate_random_elements(seed)
            .map(|vec| HexacodeVector { vec })
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> EnumeratedOrdFiniteSetSignature for HexacodeStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.hexacode
            .list_all_elements_ordered()
            .into_iter()
            .map(|vec| HexacodeVector { vec })
            .collect()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.hexacode.element_to_enumeration(&elem.vec)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.hexacode
            .enumeration_to_element(num)
            .map(|vec| HexacodeVector { vec })
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
        HexacodeVector {
            vec: self.hexacode.zero(),
        }
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> AdditionSignature for HexacodeStructure<Set, SetB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        HexacodeVector {
            vec: self.hexacode.add(&a.vec, &b.vec),
        }
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> CancellativeAdditionSignature for HexacodeStructure<Set, SetB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(HexacodeVector {
            vec: self.hexacode.try_sub(&a.vec, &b.vec)?,
        })
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature + FiniteSetSizedSignature<6>,
    SetB: BorrowedStructure<Set>,
> TryNegateSignature for HexacodeStructure<Set, SetB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(HexacodeVector {
            vec: self.hexacode.try_neg(&a.vec)?,
        })
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
        HexacodeVector {
            vec: self.hexacode.neg(&a.vec),
        }
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
        HexacodeVector {
            vec: self.hexacode.scalar_mul(&a.vec, x),
        }
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

        for x in hexacode.list_all_elements() {
            println!("{:?}", x);
        }
    }

    #[test]
    fn test_enumeration() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        let hexacode = set.hexacode();
        let codewords = hexacode.list_all_elements_ordered();
        assert_eq!(codewords.len(), 64);
        assert_eq!(hexacode.size(), Natural::from(64usize));

        // codewords are all valid
        for v in &codewords {
            println!("{:?}", v);
            assert!(hexacode.validate_element(v).is_ok());
        }

        // enumeration is correct
        for (i, v) in codewords.iter().enumerate() {
            assert_eq!(Natural::from(i), hexacode.element_to_enumeration(v));
            assert!(hexacode.equal(
                &hexacode.enumeration_to_element(&Natural::from(i)).unwrap(),
                v
            ));
        }
        assert!(
            hexacode
                .enumeration_to_element(&Natural::from(64usize))
                .is_none()
        );
    }
}
