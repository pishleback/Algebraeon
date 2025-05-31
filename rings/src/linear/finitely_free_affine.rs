use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_module::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use crate::structure::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;
use std::fmt::Debug;

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Set: Clone + Debug> {
    Empty,
    NonEmpty(FinitelyFreeSubmoduleCoset<Set>),
}

impl<Set: Clone + Debug> From<FinitelyFreeSubmoduleCoset<Set>>
    for FinitelyFreeSubmoduleAffineSubset<Set>
{
    fn from(submodule: FinitelyFreeSubmoduleCoset<Set>) -> Self {
        Self::NonEmpty(submodule)
    }
}

impl<Set: Clone + Debug> FinitelyFreeSubmoduleAffineSubset<Set> {
    pub fn rank(&self) -> Option<usize> {
        match self {
            FinitelyFreeSubmoduleAffineSubset::Empty => None,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => Some(coset.rank()),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty => true,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(..) => false,
        }
    }

    pub fn unwrap_to_coset(self) -> FinitelyFreeSubmoduleCoset<Set> {
        match self {
            FinitelyFreeSubmoduleAffineSubset::Empty => {
                panic!("unwrap called on empty affine subset")
            }
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset,
        }
    }

    pub fn affine_rank(&self) -> usize {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty => 0,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.submodule().rank() + 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleAffineSubsetStructure<
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
> {
    module: FinitelyFreeModuleStructure<Ring, RingB>,
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB>
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring, RingB>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> Signature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> SetSignature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB>
{
    type Set = FinitelyFreeSubmoduleAffineSubset<Ring::Set>;

    fn is_element(&self, _x: &Self::Set) -> bool {
        //TODO: better checks
        true
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB>
{
    pub fn ring(&self) -> &Ring {
        self.module().ring()
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Ring, RingB> {
        &self.module
    }

    pub fn from_affine_span(
        &self,
        mut span: Vec<&Vec<Ring::Set>>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        for v in &span {
            debug_assert_eq!(v.len(), self.module().rank());
        }
        if let Some(offset) = span.pop() {
            let linear_span = span
                .into_iter()
                .map(|v| self.module().add(v, &self.module().neg(offset)))
                .collect::<Vec<_>>();
            FinitelyFreeSubmoduleAffineSubset::from(
                self.module().cosets().from_offset_and_submodule(
                    offset,
                    self.module()
                        .submodules()
                        .span(linear_span.iter().collect()),
                ),
            )
        } else {
            FinitelyFreeSubmoduleAffineSubset::Empty
        }
    }

    pub fn affine_basis(
        &self,
        subset: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
    ) -> Vec<Vec<Ring::Set>> {
        match &subset {
            FinitelyFreeSubmoduleAffineSubset::Empty => vec![],
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                let mut affine_basis = vec![coset.offset().clone()];
                for v in coset.submodule().basis() {
                    affine_basis.push(self.module().add(&v, coset.offset()));
                }
                affine_basis
            }
        }
    }

    pub fn sum(
        &self,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        debug_assert!(self.is_element(x));
        debug_assert!(self.is_element(y));
        match (x, y) {
            (FinitelyFreeSubmoduleAffineSubset::Empty, _)
            | (_, FinitelyFreeSubmoduleAffineSubset::Empty) => {
                FinitelyFreeSubmoduleAffineSubset::Empty
            }
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                self.module().cosets().sum(&x_coset, &y_coset),
            ),
        }
    }

    pub fn intersect(
        &self,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        debug_assert!(self.is_element(x));
        debug_assert!(self.is_element(y));
        match (x, y) {
            (FinitelyFreeSubmoduleAffineSubset::Empty, _)
            | (_, FinitelyFreeSubmoduleAffineSubset::Empty) => {
                FinitelyFreeSubmoduleAffineSubset::Empty
            }
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().intersect(&x_coset, &y_coset),
        }
    }

    pub fn intersect_list(
        &self,
        xs: Vec<impl Borrow<FinitelyFreeSubmoduleAffineSubset<Ring::Set>>>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Set> {
        for x in &xs {
            debug_assert!(self.is_element(x.borrow()));
        }

        let mut i = FinitelyFreeSubmoduleAffineSubset::NonEmpty(
            self.module()
                .cosets()
                .from_submodule(self.module().submodules().full_submodule()),
        );
        for x in xs {
            i = self.intersect(&i, x.borrow());
        }
        i
    }

    pub fn contains_element(
        &self,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
        p: &Vec<Ring::Set>,
    ) -> bool {
        debug_assert!(self.is_element(x));
        debug_assert!(self.module().is_element(p));
        match x {
            FinitelyFreeSubmoduleAffineSubset::Empty => false,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                self.module().cosets().contains_element(coset, p)
            }
        }
    }

    pub fn equal_slow(
        &self,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
    ) -> bool {
        debug_assert!(self.is_element(x));
        debug_assert!(self.is_element(y));
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().equal_slow(x_coset, y_coset),
        }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>> EqSignature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB>
{
    fn equal(
        &self,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Set>,
    ) -> bool {
        debug_assert!(self.is_element(x));
        debug_assert!(self.is_element(y));
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().equal(x_coset, y_coset),
        }
    }
}
