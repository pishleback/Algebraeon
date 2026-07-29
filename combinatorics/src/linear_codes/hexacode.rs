//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::{
    OrderedSynthemePoint, OrderedSynthemePointCanonicalStructure,
};
use algebraeon_macros::CanonicalStructure;
use algebraeon_rings::{
    finite_fields::quaternary_field::{QuaternaryField, QuaternaryFieldCanonicalStructure},
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        finitely_free_module::FinitelyFreeModuleStructure,
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::FinitelyFreeSubmodule,
    },
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, FinitelyFreeModuleSignature,
        RinglikeSpecializationSignature, SemiModuleSignature, TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;
use cantor::{ArrayMap, Finite};
use std::{cmp::Ordering, sync::OnceLock};

type AmbientSpace = ConstFinitelyFreeModuleStructure<
    6,
    OrderedSynthemePointCanonicalStructure,
    OrderedSynthemePointCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
>;

struct HexacodeCache {
    subspace: FinitelyFreeSubmoduleStructure<
        OrderedSynthemePointCanonicalStructure,
        QuaternaryFieldCanonicalStructure,
        AmbientSpace,
        AmbientSpace,
    >,
}

static HEXACODE_CACHE: OnceLock<HexacodeCache> = OnceLock::new();

fn cache() -> &'static HexacodeCache {
    HEXACODE_CACHE.get_or_init(|| {
        let points = OrderedSynthemePoint::structure();
        let space = QuaternaryField::structure().into_free_module(points);
        let subspace = space.generated_submodule(vec![
            &[
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
            ],
            &[
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
                QuaternaryField::Alpha,
                QuaternaryField::Beta,
                QuaternaryField::Beta,
                QuaternaryField::Alpha,
            ],
            &[
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
pub fn space_structure() -> &'static ConstFinitelyFreeModuleStructure<
    6,
    OrderedSynthemePointCanonicalStructure,
    OrderedSynthemePointCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
> {
    cache().subspace.module()
}

/// The 3 dimensional vector subspace given by the hexacode
pub fn hexacode_subspace() -> &'static FinitelyFreeSubmodule<QuaternaryField> {
    cache().subspace.submodule()
}

/// The 3 dimensional vector subspace structure given by the hexacode
pub fn hexacode_subspace_structure() -> &'static FinitelyFreeSubmoduleStructure<
    OrderedSynthemePointCanonicalStructure,
    QuaternaryFieldCanonicalStructure,
    AmbientSpace,
    AmbientSpace,
> {
    &cache().subspace
}

pub type HexacodeVector = [QuaternaryField; 6];

pub fn is_hexacodeword(vector: &HexacodeVector) -> bool {
    todo!()
}

/// Given 3 coordinates find the unique completion to a hexacodeword
///
/// Returns None if the number of given coordinates is not 3
pub fn complete_hexacodeword_from_3(
    given: [Option<QuaternaryField>; 6],
) -> Option<ArrayMap<OrderedSynthemePoint, QuaternaryField>> {
    todo!()
}

/// Given 5 coordinates find the unique completion to a hexacodeword allowing at most 1 of the given coordinates to change
///
/// Returns None if the number of given coordinates is not 5
pub fn correct_hexacodeword_from_5(
    given: [Option<QuaternaryField>; 6],
) -> Option<ArrayMap<OrderedSynthemePoint, QuaternaryField>> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        println!("{:?}", space_structure());
        println!("{:?}", hexacode_subspace());

        for x in hexacode_subspace_structure().list_all_elements() {
            println!("{:?}", x);
        }
    }
}
