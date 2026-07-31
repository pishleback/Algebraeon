//! The hexacode on the standard ordered syntheme

use crate::linear_codes::pointed_ordered_3cycle::{
    PointedOrdered3Cycle, PointedOrdered3CycleCanonicalStructure,
};
use algebraeon_rings::{
    linear::{
        const_finitely_free_module::{
            ConstFinitelyFreeModuleStructure, RingToConstFinitelyFreeModuleSignature,
        },
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::FinitelyFreeSubmodule,
    },
    num_theory::modulo::const_naive::Modulo,
    structure::{FinitelyFreeModuleSignature, MetaAdditionSignature},
};
use algebraeon_sets::sets::{Function, SetToConstSizeFunctionsToSignature};
use algebraeon_structures::*;
use std::{ops::Add, sync::OnceLock};

type F3 = Modulo<3>;
pub const ZERO: F3 = F3::new(0);
pub const PLUS: F3 = F3::new(1);
pub const MINUS: F3 = F3::new(2);

pub type LabelledPoints<Elem> = Function<4, PointedOrdered3Cycle, Elem>;

#[derive(Clone, PartialEq, Eq)]
pub struct TetracodeVector(LabelledPoints<F3>);

impl std::fmt::Debug for TetracodeVector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let f3_fmt = |x: &F3| -> &'static str {
            if *x == ZERO {
                "0"
            } else if *x == PLUS {
                "+"
            } else {
                debug_assert_eq!(*x, MINUS);
                "-"
            }
        };
        f.debug_tuple("Vector")
            .field(&format_args!(
                "{} {}{}{}",
                f3_fmt(self.at(&PointedOrdered3Cycle::P)),
                f3_fmt(self.at(&PointedOrdered3Cycle::C1)),
                f3_fmt(self.at(&PointedOrdered3Cycle::C2)),
                f3_fmt(self.at(&PointedOrdered3Cycle::C3))
            ))
            .finish()
    }
}

impl From<LabelledPoints<F3>> for TetracodeVector {
    fn from(value: LabelledPoints<F3>) -> Self {
        Self(value)
    }
}

impl From<TetracodeVector> for LabelledPoints<F3> {
    fn from(value: TetracodeVector) -> Self {
        value.0
    }
}

impl<'a> From<&'a TetracodeVector> for &'a LabelledPoints<F3> {
    fn from(value: &'a TetracodeVector) -> Self {
        &value.0
    }
}

impl From<[F3; 4]> for TetracodeVector {
    fn from(value: [F3; 4]) -> Self {
        LabelledPoints::from(value).into()
    }
}

impl TetracodeVector {
    pub fn new(f: impl FnMut(PointedOrdered3Cycle) -> F3) -> TetracodeVector {
        TetracodeVector(LabelledPoints::new(f))
    }

    pub fn at(&self, point: &PointedOrdered3Cycle) -> &F3 {
        &self.0[point.element_to_enumeration().try_into().unwrap()]
    }

    pub fn at_mut(&mut self, point: &PointedOrdered3Cycle) -> &mut F3 {
        &mut self.0[point.element_to_enumeration().try_into().unwrap()]
    }
}

impl Add<&TetracodeVector> for &TetracodeVector {
    type Output = TetracodeVector;

    fn add(self, other: &TetracodeVector) -> Self::Output {
        LabelledPoints::<F3>::new(|p| F3::add(self.0.image(&p), other.0.image(&p))).into()
    }
}

type AmbientSpace = ConstFinitelyFreeModuleStructure<
    4,
    PointedOrdered3CycleCanonicalStructure,
    PointedOrdered3CycleCanonicalStructure,
    <F3 as MetaType>::Signature,
    <F3 as MetaType>::Signature,
>;

struct TetracodeCache {
    subspace: FinitelyFreeSubmoduleStructure<
        PointedOrdered3CycleCanonicalStructure,
        <F3 as MetaType>::Signature,
        AmbientSpace,
        AmbientSpace,
    >,
}

static TETRACODE_CACHE: OnceLock<TetracodeCache> = OnceLock::new();

fn cache() -> &'static TetracodeCache {
    TETRACODE_CACHE.get_or_init(|| {
        let space = F3::structure().into_free_module(PointedOrdered3Cycle::structure());
        let subspace = space.generated_submodule(vec![
            &PointedOrdered3Cycle::structure()
                .const_size_functions_to(&F3::structure())
                .function(|x| match x {
                    PointedOrdered3Cycle::P => ZERO,
                    PointedOrdered3Cycle::C1 => PLUS,
                    PointedOrdered3Cycle::C2 => PLUS,
                    PointedOrdered3Cycle::C3 => PLUS,
                })
                .unwrap(),
            &PointedOrdered3Cycle::structure()
                .const_size_functions_to(&F3::structure())
                .function(|x| match x {
                    PointedOrdered3Cycle::P => PLUS,
                    PointedOrdered3Cycle::C1 => ZERO,
                    PointedOrdered3Cycle::C2 => PLUS,
                    PointedOrdered3Cycle::C3 => MINUS,
                })
                .unwrap(),
        ]);
        let subspace_structure = FinitelyFreeSubmoduleStructure::new(space, subspace);
        TetracodeCache {
            subspace: subspace_structure,
        }
    })
}

/// The 4 dimensional vector space structure over F3 with basis given by the points of a pointed ordered 3-cycle
pub fn space_structure() -> &'static AmbientSpace {
    cache().subspace.module()
}

/// The 2 dimensional vector subspace given by the tetracode
pub fn tetracode_subspace() -> &'static FinitelyFreeSubmodule<F3> {
    cache().subspace.submodule()
}

/// The 2 dimensional vector subspace structure given by the tetracode
pub fn tetracode_subspace_structure() -> &'static FinitelyFreeSubmoduleStructure<
    PointedOrdered3CycleCanonicalStructure,
    <F3 as MetaType>::Signature,
    AmbientSpace,
    AmbientSpace,
> {
    &cache().subspace
}

pub fn all_tetracodewords() -> Vec<TetracodeVector> {
    tetracode_subspace_structure()
        .list_all_elements()
        .into_iter()
        .map(TetracodeVector::from)
        .collect()
}

impl TetracodeVector {
    pub fn is_tetracodeword(&self) -> bool {
        tetracode_subspace_structure().is_element(&self.0)
    }
}

/// Given 2 coordinates find the unique completion to a tetracodeword
///
/// Returns None if the number of given coordinates is not 2
pub fn complete_tetracodeword_from_2(given: LabelledPoints<Option<F3>>) -> Option<TetracodeVector> {
    if given
        .iter()
        .map(|x| if x.is_some() { 1 } else { 0 })
        .sum::<usize>()
        != 2
    {
        None
    } else {
        for codeword in tetracode_subspace_structure().list_all_elements() {
            if given.iter().enumerate().all(|(i, x)| {
                if let Some(x) = x {
                    x == &codeword[i]
                } else {
                    true
                }
            }) {
                return Some(TetracodeVector(codeword));
            }
        }
        unreachable!()
    }
}

/// Given 4 coordinates find the unique correction to a tetracodeword allowing at most 1 of the given coordinates to change
pub fn correct_tetracodeword_from_4(given: LabelledPoints<F3>) -> TetracodeVector {
    for codeword in tetracode_subspace_structure().list_all_elements() {
        if given
            .iter()
            .enumerate()
            .map(|(i, x)| if x == &codeword[i] { 1 } else { 0 })
            .sum::<usize>()
            >= 3
        {
            return TetracodeVector(codeword);
        }
    }
    unreachable!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_tetracodeword() {
        for x in all_tetracodewords() {
            assert!(x.is_tetracodeword());
        }

        let mut c = 0;
        for x in space_structure()
            .list_all_elements()
            .into_iter()
            .map(TetracodeVector)
        {
            if x.is_tetracodeword() {
                c += 1;
            }
        }
        assert_eq!(c, 9);
    }

    #[test]
    fn test_complete_2() {
        assert!(complete_tetracodeword_from_2([None, None, None, None].into()).is_none());

        assert!(complete_tetracodeword_from_2([Some(ZERO), None, None, None].into()).is_none());

        assert!(
            complete_tetracodeword_from_2([Some(ZERO), Some(ZERO), None, None].into()).is_some()
        );

        assert!(
            complete_tetracodeword_from_2([Some(ZERO), Some(ZERO), Some(ZERO), None].into())
                .is_none()
        );

        assert!(
            complete_tetracodeword_from_2([Some(ZERO), Some(ZERO), Some(ZERO), Some(ZERO)].into())
                .is_none()
        );

        assert_eq!(
            complete_tetracodeword_from_2([Some(ZERO), Some(ZERO), None, None].into()).unwrap(),
            [ZERO, ZERO, ZERO, ZERO].into()
        );

        assert_eq!(
            complete_tetracodeword_from_2([None, Some(ZERO), None, Some(ZERO)].into()).unwrap(),
            [ZERO, ZERO, ZERO, ZERO].into()
        );

        assert_eq!(
            complete_tetracodeword_from_2([None, Some(PLUS), None, Some(MINUS)].into()).unwrap(),
            [MINUS, PLUS, ZERO, MINUS].into()
        );

        assert_eq!(
            complete_tetracodeword_from_2([Some(PLUS), None, None, Some(MINUS)].into()).unwrap(),
            [PLUS, ZERO, PLUS, MINUS].into()
        );
    }

    #[test]
    fn test_correct_4() {
        assert_eq!(
            correct_tetracodeword_from_4([ZERO, ZERO, ZERO, ZERO].into()),
            [ZERO, ZERO, ZERO, ZERO].into()
        );

        assert_eq!(
            correct_tetracodeword_from_4([ZERO, PLUS, ZERO, ZERO].into()),
            [ZERO, ZERO, ZERO, ZERO].into()
        );

        assert_eq!(
            correct_tetracodeword_from_4([MINUS, ZERO, ZERO, ZERO].into()),
            [ZERO, ZERO, ZERO, ZERO].into()
        );

        assert_eq!(
            correct_tetracodeword_from_4([MINUS, MINUS, PLUS, PLUS].into()),
            [MINUS, MINUS, PLUS, ZERO].into()
        );

        assert_eq!(
            correct_tetracodeword_from_4([MINUS, PLUS, PLUS, PLUS].into()),
            [ZERO, PLUS, PLUS, PLUS].into()
        );
    }
}
