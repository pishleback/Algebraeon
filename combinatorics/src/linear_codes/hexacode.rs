//! The hexacode on the standard ordered syntheme

use crate::linear_codes::ordered_syntheme::{
    OrderedSynthemePair, OrderedSynthemePoint, OrderedSynthemePointCanonicalStructure,
    OrderedSynthemeSide,
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
    structure::{FinitelyFreeModuleSignature, MetaAdditionSignature},
};
use algebraeon_sets::sets::Function;
use algebraeon_structures::*;
use std::{ops::Add, sync::OnceLock};

type F4 = QuaternaryField;
type F4Structure = QuaternaryFieldCanonicalStructure;

pub type LabelledPoints<Elem> = Function<6, OrderedSynthemePoint, Elem>;

#[derive(Clone, PartialEq, Eq)]
pub struct HexacodeVector(LabelledPoints<F4>);

impl std::fmt::Debug for HexacodeVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write;
        let mut s = String::new();
        for (i, pair) in OrderedSynthemePair::list_all_elements_ordered()
            .into_iter()
            .enumerate()
        {
            if i != 0 {
                write!(&mut s, " ")?;
            }
            for side in OrderedSynthemeSide::list_all_elements_ordered() {
                let p = OrderedSynthemePoint { side, pair };
                write!(&mut s, "{}", self.at(&p))?;
            }
        }
        f.debug_tuple("Vector")
            .field(&format_args!("{}", s))
            .finish()
    }
}

impl From<LabelledPoints<F4>> for HexacodeVector {
    fn from(value: LabelledPoints<F4>) -> Self {
        Self(value)
    }
}

impl From<HexacodeVector> for LabelledPoints<F4> {
    fn from(value: HexacodeVector) -> Self {
        value.0
    }
}

impl<'a> From<&'a HexacodeVector> for &'a LabelledPoints<F4> {
    fn from(value: &'a HexacodeVector) -> Self {
        &value.0
    }
}

impl From<[F4; 6]> for HexacodeVector {
    fn from(value: [F4; 6]) -> Self {
        LabelledPoints::from(value).into()
    }
}

impl HexacodeVector {
    pub fn new(f: impl FnMut(OrderedSynthemePoint) -> F4) -> HexacodeVector {
        HexacodeVector(LabelledPoints::new(f))
    }

    pub fn at(&self, point: &OrderedSynthemePoint) -> &F4 {
        &self.0[point.element_to_enumeration().try_into().unwrap()]
    }

    pub fn at_mut(&mut self, point: &OrderedSynthemePoint) -> &mut F4 {
        &mut self.0[point.element_to_enumeration().try_into().unwrap()]
    }
}

impl Add<&HexacodeVector> for &HexacodeVector {
    type Output = HexacodeVector;

    fn add(self, other: &HexacodeVector) -> Self::Output {
        LabelledPoints::<F4>::new(|p| F4::add(self.0.image(&p), other.0.image(&p))).into()
    }
}

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

pub fn all_hexacodewords() -> Vec<HexacodeVector> {
    hexacode_subspace_structure()
        .list_all_elements()
        .into_iter()
        .map(HexacodeVector::from)
        .collect()
}

impl HexacodeVector {
    pub fn is_hexacodeword(&self) -> bool {
        hexacode_subspace_structure().is_element(&self.0)
    }
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
                return Some(HexacodeVector(codeword));
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
                return Some(HexacodeVector(codeword));
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
        for x in all_hexacodewords() {
            assert!(x.is_hexacodeword());
        }

        let mut c = 0;
        for x in space_structure()
            .list_all_elements()
            .into_iter()
            .map(HexacodeVector)
        {
            if x.is_hexacodeword() {
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
