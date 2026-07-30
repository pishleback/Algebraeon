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
use algebraeon_sets::sets::Function;
use algebraeon_structures::*;
use std::sync::OnceLock;

type F4 = QuaternaryField;
type F4Structure = QuaternaryFieldCanonicalStructure;

pub type LabelledPoints<Elem> = Function<6, OrderedSynthemePoint, Elem>;
pub type HexacodeVector = LabelledPoints<F4>;

type AmbientSpace = ConstFinitelyFreeModuleStructure<
    6,
    OrderedSynthemePointCanonicalStructure,
    OrderedSynthemePointCanonicalStructure,
    F4Structure,
    F4Structure,
>;

struct HexacodeCache {
    subspace: FinitelyFreeSubmoduleStructure<
        OrderedSynthemePointCanonicalStructure,
        F4Structure,
        AmbientSpace,
        AmbientSpace,
    >,
}

static HEXACODE_CACHE: OnceLock<HexacodeCache> = OnceLock::new();

fn cache() -> &'static HexacodeCache {
    HEXACODE_CACHE.get_or_init(|| {
        let points = OrderedSynthemePoint::structure();
        let space = F4::structure().into_free_module(points);
        let subspace = space.generated_submodule(vec![
            &[
                F4::Alpha,
                F4::Beta,
                F4::Beta,
                F4::Alpha,
                F4::Beta,
                F4::Alpha,
            ]
            .into(),
            &[
                F4::Beta,
                F4::Alpha,
                F4::Alpha,
                F4::Beta,
                F4::Beta,
                F4::Alpha,
            ]
            .into(),
            &[
                F4::Beta,
                F4::Alpha,
                F4::Beta,
                F4::Alpha,
                F4::Alpha,
                F4::Beta,
            ]
            .into(),
        ]);
        let subspace_structure = FinitelyFreeSubmoduleStructure::new(space, subspace);
        HexacodeCache {
            subspace: subspace_structure,
        }
    })
}

/// The 6 dimensional vector space structure over F4 with basis given by the points of an ordered syntheme
pub fn space_structure() -> &'static AmbientSpace {
    cache().subspace.module()
}

/// The 3 dimensional vector subspace given by the hexacode
pub fn hexacode_subspace() -> &'static FinitelyFreeSubmodule<F4> {
    cache().subspace.submodule()
}

/// The 3 dimensional vector subspace structure given by the hexacode
pub fn hexacode_subspace_structure() -> &'static FinitelyFreeSubmoduleStructure<
    OrderedSynthemePointCanonicalStructure,
    <F4 as MetaType>::Signature,
    AmbientSpace,
    AmbientSpace,
> {
    &cache().subspace
}

pub fn is_hexacodeword(vector: &HexacodeVector) -> bool {
    hexacode_subspace_structure().is_element(vector)
}

/// Given 3 coordinates find the unique completion to a hexacodeword
///
/// Returns None if the number of given coordinates is not 3
pub fn complete_hexacodeword_from_3(given: LabelledPoints<Option<F4>>) -> Option<HexacodeVector> {
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

/// Given 5 coordinates find the unique correction to a hexacodeword allowing at most 1 of the given coordinates to change
///
/// Returns None if the number of given coordinates is not 5
pub fn correct_hexacodeword_from_5(given: LabelledPoints<Option<F4>>) -> Option<HexacodeVector> {
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
    use crate::linear_codes::ordered_syntheme::{OrderedSynthemePair, OrderedSynthemeSide};

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
    fn construct_vector() {
        let v = HexacodeVector::new(|p| match (p.pair, p.side) {
            (OrderedSynthemePair::Left, OrderedSynthemeSide::Left) => F4::Zero,
            (OrderedSynthemePair::Left, OrderedSynthemeSide::Right) => F4::One,
            (OrderedSynthemePair::Middle, OrderedSynthemeSide::Left) => F4::Alpha,
            (OrderedSynthemePair::Middle, OrderedSynthemeSide::Right) => F4::Beta,
            (OrderedSynthemePair::Right, OrderedSynthemeSide::Left) => F4::One,
            (OrderedSynthemePair::Right, OrderedSynthemeSide::Right) => F4::Zero,
        });
        assert_eq!(
            v,
            [F4::Zero, F4::One, F4::Alpha, F4::Beta, F4::One, F4::Zero].into()
        );
    }

    #[test]
    fn test_complete_3() {
        assert_eq!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                    None,
                    None,
                ]
                .into()
            )
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero].into()
        );

        assert_eq!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::One),
                    Some(F4::Zero),
                    None,
                    None,
                    None,
                ]
                .into()
            )
            .unwrap(),
            [F4::Zero, F4::One, F4::Zero, F4::One, F4::Alpha, F4::Beta].into()
        );

        assert_eq!(
            complete_hexacodeword_from_3(
                [
                    None,
                    None,
                    Some(F4::Alpha),
                    None,
                    Some(F4::Beta),
                    Some(F4::Beta),
                ]
                .into()
            )
            .unwrap(),
            [F4::One, F4::One, F4::Alpha, F4::Alpha, F4::Beta, F4::Beta].into()
        );

        assert!(
            complete_hexacodeword_from_3([None, None, None, None, None, None].into()).is_none()
        );
        assert!(
            complete_hexacodeword_from_3([Some(F4::Zero), None, None, None, None, None].into())
                .is_none()
        );
        assert!(
            complete_hexacodeword_from_3(
                [Some(F4::Zero), Some(F4::Zero), None, None, None, None].into()
            )
            .is_none()
        );
        assert!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                    None,
                    None
                ]
                .into()
            )
            .is_some()
        );
        assert!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                    None
                ]
                .into()
            )
            .is_none()
        );
        assert!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None
                ]
                .into()
            )
            .is_none()
        );
        assert!(
            complete_hexacodeword_from_3(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero)
                ]
                .into()
            )
            .is_none()
        );
    }

    #[test]
    fn test_correct_5() {
        assert_eq!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                ]
                .into()
            )
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero].into()
        );

        assert_eq!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::One),
                    Some(F4::Zero),
                    None,
                ]
                .into()
            )
            .unwrap(),
            [F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero, F4::Zero].into()
        );

        assert_eq!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::One),
                    Some(F4::Alpha),
                    None,
                    Some(F4::Zero),
                    Some(F4::Beta),
                    Some(F4::Beta),
                ]
                .into()
            )
            .unwrap(),
            [F4::One, F4::Alpha, F4::Beta, F4::Zero, F4::Beta, F4::Zero].into()
        );

        assert!(correct_hexacodeword_from_5([None, None, None, None, None, None].into()).is_none());
        assert!(
            correct_hexacodeword_from_5([Some(F4::Zero), None, None, None, None, None].into())
                .is_none()
        );
        assert!(
            correct_hexacodeword_from_5(
                [Some(F4::Zero), Some(F4::Zero), None, None, None, None].into()
            )
            .is_none()
        );
        assert!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                    None,
                    None
                ]
                .into()
            )
            .is_none()
        );
        assert!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None,
                    None
                ]
                .into()
            )
            .is_none()
        );
        assert!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    None
                ]
                .into()
            )
            .is_some()
        );
        assert!(
            correct_hexacodeword_from_5(
                [
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero),
                    Some(F4::Zero)
                ]
                .into()
            )
            .is_none()
        );
    }
}
