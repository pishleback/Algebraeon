//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::{
    OrderedSynthemePoint, OrderedSynthemePointCanonicalStructure,
};
use algebraeon_rings::{
    finite_fields::quaternary_field::{QuaternaryField, QuaternaryFieldCanonicalStructure},
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::FinitelyFreeSubmodule,
    },
    structure::FinitelyFreeModuleSignature,
};
use algebraeon_structures::*;
use std::sync::OnceLock;

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
    hexacode_subspace_structure().is_element(vector)
}

/// Given 3 coordinates find the unique completion to a hexacodeword
///
/// Returns None if the number of given coordinates is not 3
pub fn complete_hexacodeword_from_3(given: [Option<QuaternaryField>; 6]) -> Option<HexacodeVector> {
    if given
        .iter()
        .map(|x| if x.is_some() { 1 } else { 0 })
        .sum::<usize>()
        != 3
    {
        None
    } else {
        for codeword in hexacode_subspace_structure().list_all_elements() {
            if given.iter().enumerate().all(|(i, x)| {
                if let Some(x) = x {
                    x == &codeword[i]
                } else {
                    true
                }
            }) {
                return Some(codeword);
            }
        }
        unreachable!()
    }
}

/// Given 5 coordinates find the unique completion to a hexacodeword allowing at most 1 of the given coordinates to change
///
/// Returns None if the number of given coordinates is not 5
pub fn correct_hexacodeword_from_5(
    given: [Option<QuaternaryField>; 6],
) -> Option<[QuaternaryField; 6]> {
    if given
        .iter()
        .map(|x| if x.is_some() { 1 } else { 0 })
        .sum::<usize>()
        != 5
    {
        None
    } else {
        for codeword in hexacode_subspace_structure().list_all_elements() {
            if given
                .iter()
                .enumerate()
                .map(|(i, x)| {
                    if let Some(x) = x
                        && x == &codeword[i]
                    {
                        1
                    } else {
                        0
                    }
                })
                .sum::<usize>()
                >= 4
            {
                return Some(codeword);
            }
        }
        unreachable!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_hexacodeword() {
        for x in hexacode_subspace_structure().list_all_elements() {
            assert!(is_hexacodeword(&x));
        }

        let mut c = 0;
        for x in space_structure().list_all_elements() {
            if is_hexacodeword(&x) {
                c += 1;
            }
        }
        assert_eq!(c, 64);
    }

    #[test]
    fn test_complete_3() {
        type F4 = QuaternaryField;

        assert_eq!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
                None,
                None,
            ])
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero]
        );

        assert_eq!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::One),
                Some(F4::Zero),
                None,
                None,
                None,
            ])
            .unwrap(),
            [F4::Zero, F4::One, F4::Zero, F4::One, F4::Alpha, F4::Beta]
        );

        assert_eq!(
            complete_hexacodeword_from_3([
                None,
                None,
                Some(F4::Alpha),
                None,
                Some(F4::Beta),
                Some(F4::Beta),
            ])
            .unwrap(),
            [F4::One, F4::One, F4::Alpha, F4::Alpha, F4::Beta, F4::Beta]
        );

        assert!(complete_hexacodeword_from_3([None, None, None, None, None, None]).is_none());
        assert!(
            complete_hexacodeword_from_3([Some(F4::Zero), None, None, None, None, None]).is_none()
        );
        assert!(
            complete_hexacodeword_from_3([Some(F4::Zero), Some(F4::Zero), None, None, None, None])
                .is_none()
        );
        assert!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
                None,
                None
            ])
            .is_some()
        );
        assert!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
                None
            ])
            .is_none()
        );
        assert!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None
            ])
            .is_none()
        );
        assert!(
            complete_hexacodeword_from_3([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero)
            ])
            .is_none()
        );
    }

    #[test]
    fn test_correct_5() {
        type F4 = QuaternaryField;

        assert_eq!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
            ])
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero]
        );

        assert_eq!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::One),
                Some(F4::Zero),
                None,
            ])
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero]
        );

        assert_eq!(
            correct_hexacodeword_from_5([
                Some(F4::One),
                Some(F4::Alpha),
                None,
                Some(F4::Zero),
                Some(F4::Beta),
                Some(F4::Beta),
            ])
            .unwrap(),
            [F4::One, F4::Alpha, F4::Beta, F4::Zero, F4::Beta, F4::Zero]
        );

        assert!(correct_hexacodeword_from_5([None, None, None, None, None, None]).is_none());
        assert!(
            correct_hexacodeword_from_5([Some(F4::Zero), None, None, None, None, None]).is_none()
        );
        assert!(
            correct_hexacodeword_from_5([Some(F4::Zero), Some(F4::Zero), None, None, None, None])
                .is_none()
        );
        assert!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
                None,
                None
            ])
            .is_none()
        );
        assert!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None,
                None
            ])
            .is_none()
        );
        assert!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                None
            ])
            .is_some()
        );
        assert!(
            correct_hexacodeword_from_5([
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero),
                Some(F4::Zero)
            ])
            .is_none()
        );
    }
}
