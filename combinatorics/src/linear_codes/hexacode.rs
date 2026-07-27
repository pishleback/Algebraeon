//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::{
    OrderedSynthemePoint, OrderedSynthemePointCanonicalStructure,
};
use algebraeon_macros::CanonicalStructure;
use algebraeon_rings::{
    finite_fields::quaternary_field::{QuaternaryField, QuaternaryFieldCanonicalStructure},
    linear::{
        finitely_free_module::{FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature},
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::FinitelyFreeSubmodule,
    },
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, RinglikeSpecializationSignature, SemiModuleSignature,
        TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;
use cantor::{ArrayMap, Finite};
use std::{cmp::Ordering, sync::OnceLock};

struct HexacodeCache {
    subspace: FinitelyFreeSubmoduleStructure<
        OrderedSynthemePointCanonicalStructure,
        OrderedSynthemePointCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
        FinitelyFreeModuleStructure<
            OrderedSynthemePointCanonicalStructure,
            OrderedSynthemePointCanonicalStructure,
            QuaternaryFieldCanonicalStructure,
            QuaternaryFieldCanonicalStructure,
        >,
    >,
}

static HEXACODE_CACHE: OnceLock<HexacodeCache> = OnceLock::new();

fn cache() -> &'static HexacodeCache {
    HEXACODE_CACHE.get_or_init(|| {
        let points = OrderedSynthemePoint::structure();
        let space = QuaternaryField::structure().into_free_module(points);
        let subspace = space.generated_submodule(vec![
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
        let subspace_structure = FinitelyFreeSubmoduleStructure::new(space, subspace);
        HexacodeCache {
            subspace: subspace_structure,
        }
    })
}

/// The 6 dimensional vector space structure over F4 with basis given by the points of an ordered syntheme
pub fn space_structure() -> &'static FinitelyFreeModuleStructure<
    OrderedSynthemePointCanonicalStructure,
    OrderedSynthemePointCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
> {
    cache().subspace.module()
}

/// The 3 dimensional vector subspace given by the hexacode
pub fn subspace() -> &'static FinitelyFreeSubmodule<QuaternaryField> {
    cache().subspace.submodule()
}

/// The 3 dimensional vector subspace structure given by the hexacode
pub fn subspace_structure() -> &'static FinitelyFreeSubmoduleStructure<
    OrderedSynthemePointCanonicalStructure,
    OrderedSynthemePointCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    FinitelyFreeModuleStructure<
        OrderedSynthemePointCanonicalStructure,
        OrderedSynthemePointCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
    >,
> {
    &cache().subspace
}

/// Any vector in F4^6
#[derive(CanonicalStructure, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Finite)]
#[canonical_structure(eq, partial_ord, ord, finite, ord_finite)]
pub struct HexacodeVector {
    a: QuaternaryField,
    b: QuaternaryField,
    c: QuaternaryField,
    d: QuaternaryField,
    e: QuaternaryField,
    f: QuaternaryField,
}

impl From<[QuaternaryField; 6]> for HexacodeVector {
    fn from(value: [QuaternaryField; 6]) -> Self {
        Self {
            a: value[0],
            b: value[1],
            c: value[2],
            d: value[3],
            e: value[4],
            f: value[5],
        }
    }
}

impl From<HexacodeVector> for [QuaternaryField; 6] {
    fn from(value: HexacodeVector) -> Self {
        [value.a, value.b, value.c, value.d, value.e, value.f]
    }
}

impl TryFrom<Vec<QuaternaryField>> for HexacodeVector {
    type Error = ();

    fn try_from(value: Vec<QuaternaryField>) -> Result<Self, Self::Error> {
        if value.len() == 6 {
            let mut value = value.into_iter();
            Ok(HexacodeVector {
                a: value.next().unwrap(),
                b: value.next().unwrap(),
                c: value.next().unwrap(),
                d: value.next().unwrap(),
                e: value.next().unwrap(),
                f: value.next().unwrap(),
            })
        } else {
            Err(())
        }
    }
}

impl From<HexacodeVector> for Vec<QuaternaryField> {
    fn from(value: HexacodeVector) -> Self {
        vec![value.a, value.b, value.c, value.d, value.e, value.f]
    }
}

impl<'a> From<&'a HexacodeVector> for Vec<&'a QuaternaryField> {
    fn from(value: &'a HexacodeVector) -> Self {
        vec![&value.a, &value.b, &value.c, &value.d, &value.e, &value.f]
    }
}

impl RinglikeSpecializationSignature for HexacodeVectorCanonicalStructure {}

impl ZeroSignature for HexacodeVectorCanonicalStructure {
    fn zero(&self) -> Self::Elem {
        space_structure().zero().try_into().unwrap()
    }
}

impl AdditionSignature for HexacodeVectorCanonicalStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        space_structure()
            .add(&a.clone().into(), &b.clone().into())
            .try_into()
            .unwrap()
    }
}

impl CancellativeAdditionSignature for HexacodeVectorCanonicalStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(
            space_structure()
                .try_sub(&a.clone().into(), &b.clone().into())?
                .try_into()
                .unwrap(),
        )
    }
}

impl TryNegateSignature for HexacodeVectorCanonicalStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(
            space_structure()
                .try_neg(&a.clone().into())?
                .try_into()
                .unwrap(),
        )
    }
}

impl AdditiveMonoidSignature for HexacodeVectorCanonicalStructure {}

impl AdditiveGroupSignature for HexacodeVectorCanonicalStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        space_structure().neg(&a.clone().into()).try_into().unwrap()
    }
}

impl SemiModuleSignature<QuaternaryFieldCanonicalStructure> for HexacodeVectorCanonicalStructure {
    fn ring(&self) -> &QuaternaryFieldCanonicalStructure {
        space_structure().ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &QuaternaryField) -> Self::Elem {
        space_structure()
            .scalar_mul(&a.clone().into(), x)
            .try_into()
            .unwrap()
    }
}

/// A hexacode codeword in F4^6
#[derive(CanonicalStructure, Debug, Clone, PartialEq, Eq)]
#[canonical_structure(eq, partial_ord, ord, ord_finite)]
pub struct Hexacodeword {
    vec: HexacodeVector,
}
cantor::impl_concrete_finite!(Hexacodeword);

impl PartialOrd for Hexacodeword {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Hexacodeword {
    fn cmp(&self, other: &Self) -> Ordering {
        subspace_structure().cmp(&self.clone().into(), &other.clone().into())
    }
}

unsafe impl Finite for Hexacodeword {
    const COUNT: usize = 64;

    fn index_of(value: Self) -> usize {
        subspace_structure()
            .element_to_enumeration(&value.into())
            .try_into()
            .unwrap()
    }

    fn nth(index: usize) -> Option<Self> {
        subspace_structure()
            .enumeration_to_element(&Natural::from(index))
            .map(|vec| vec.try_into().unwrap())
    }
}

impl TryFrom<Vec<QuaternaryField>> for Hexacodeword {
    type Error = ();

    fn try_from(value: Vec<QuaternaryField>) -> Result<Self, Self::Error> {
        Ok(Self {
            vec: HexacodeVector::try_from(value)?,
        })
    }
}

impl From<Hexacodeword> for Vec<QuaternaryField> {
    fn from(value: Hexacodeword) -> Self {
        value.vec.into()
    }
}

impl<'a> From<&'a Hexacodeword> for Vec<&'a QuaternaryField> {
    fn from(value: &'a Hexacodeword) -> Self {
        (&value.vec).into()
    }
}

impl CountableSetSignature for HexacodewordCanonicalStructure {
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        subspace_structure()
            .generate_all_elements()
            .map(|vec| Hexacodeword {
                vec: vec.try_into().unwrap(),
            })
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        subspace_structure()
            .generate_all_elements()
            .map(|vec| Hexacodeword {
                vec: vec.try_into().unwrap(),
            })
    }
}

impl FiniteSetSignature for HexacodewordCanonicalStructure {
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.generate_all_elements().collect()
    }

    fn size(&self) -> Natural {
        Natural::from(64usize)
    }
}

impl RinglikeSpecializationSignature for HexacodewordCanonicalStructure {}

impl ZeroSignature for HexacodewordCanonicalStructure {
    fn zero(&self) -> Self::Elem {
        space_structure().zero().try_into().unwrap()
    }
}

impl AdditionSignature for HexacodewordCanonicalStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        space_structure()
            .add(&a.clone().into(), &b.clone().into())
            .try_into()
            .unwrap()
    }
}

impl CancellativeAdditionSignature for HexacodewordCanonicalStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(
            space_structure()
                .try_sub(&a.clone().into(), &b.clone().into())?
                .try_into()
                .unwrap(),
        )
    }
}

impl TryNegateSignature for HexacodewordCanonicalStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(
            space_structure()
                .try_neg(&a.clone().into())?
                .try_into()
                .unwrap(),
        )
    }
}

impl AdditiveMonoidSignature for HexacodewordCanonicalStructure {}

impl AdditiveGroupSignature for HexacodewordCanonicalStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        space_structure().neg(&a.clone().into()).try_into().unwrap()
    }
}

impl SemiModuleSignature<QuaternaryFieldCanonicalStructure> for HexacodewordCanonicalStructure {
    fn ring(&self) -> &QuaternaryFieldCanonicalStructure {
        space_structure().ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &QuaternaryField) -> Self::Elem {
        space_structure()
            .scalar_mul(&a.clone().into(), x)
            .try_into()
            .unwrap()
    }
}

/// Given 3 coordinates find the unique completion to a hexacodeword
///
/// Returns None if the number of given coordinates is not 3
pub fn complete_hexacodeword_from_3(
    given: ArrayMap<OrderedSynthemePoint, Option<QuaternaryField>>,
) -> Option<ArrayMap<OrderedSynthemePoint, QuaternaryField>> {
    todo!()
}

/// Given 5 coordinates find the unique completion to a hexacodeword allowing at most 1 of the given coordinates to change
///
/// Returns None if the number of given coordinates is not 5
pub fn correct_hexacodeword_from_5(
    given: ArrayMap<OrderedSynthemePoint, Option<QuaternaryField>>,
) -> Option<ArrayMap<OrderedSynthemePoint, QuaternaryField>> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use cantor::ArrayMap;

    #[test]
    fn test() {
        println!("{:?}", space_structure());
        println!("{:?}", subspace());

        let z = ArrayMap::new(|vec: Hexacodeword| 4);

        for x in Hexacodeword::list_all_elements() {
            println!("{:?} -> {:?}", x, z[x.clone()]);
        }
    }
}
