use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use crate::structure::{FinitelyFreeModuleSignature, ModuleSignature};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Ring: ReducedHermiteAlgorithmSignature> {
    Empty(FinitelyFreeModuleStructure<Ring>),
    NonEmpty(FinitelyFreeSubmoduleCoset<Ring>),
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleAffineSubset<Ring> {
    pub fn ring(&self) -> &Ring {
        self.module().ring()
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty(module) => module,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.module(),
        }
    }

    pub fn module_rank(&self) -> usize {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty(module) => module.rank(),
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.module_rank(),
        }
    }

    pub fn subset_rank(&self) -> Option<usize> {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty { .. } => None,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => Some(coset.coset_rank()),
        }
    }

    pub fn new_empty(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self::Empty(module)
    }

    pub fn from_coset(coset: FinitelyFreeSubmoduleCoset<Ring>) -> Self {
        Self::NonEmpty(coset)
    }

    pub fn is_empty(&self) -> bool {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty(..) => true,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(..) => false,
        }
    }

    pub fn unwrap_to_coset(self) -> FinitelyFreeSubmoduleCoset<Ring> {
        match self {
            FinitelyFreeSubmoduleAffineSubset::Empty(..) => {
                panic!("unwrap called on empty affine subset")
            }
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset,
        }
    }

    pub fn affine_rank(&self) -> usize {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty(..) => 0,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                coset.submodule().submodule_rank() + 1
            }
        }
    }

    pub fn affine_basis(&self) -> Vec<Vec<Ring::Set>> {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty(..) => vec![],
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                let mut affine_basis = vec![coset.offset().clone()];
                for v in coset.submodule().basis() {
                    affine_basis.push(self.module().add(&v, coset.offset()));
                }
                affine_basis
            }
        }
    }

    pub fn add(x: &Self, y: &Self) -> Self {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        let module_rank = x.module_rank();
        assert_eq!(module_rank, y.module_rank());
        match (x, y) {
            (Self::Empty(_), _) | (_, Self::Empty(_)) => Self::Empty(module),
            (Self::NonEmpty(x_coset), Self::NonEmpty(y_coset)) => {
                Self::NonEmpty(FinitelyFreeSubmoduleCoset::add(&x_coset, &y_coset))
            }
        }
    }

    pub fn intersect(x: &Self, y: &Self) -> Self {
        let module = common_structure::<FinitelyFreeModuleStructure<Ring>>(x.module(), y.module());
        let module_rank = x.module_rank();
        assert_eq!(module_rank, y.module_rank());
        match (x, y) {
            (Self::Empty(_), _) | (_, Self::Empty(_)) => Self::Empty(module),
            (Self::NonEmpty(x_coset), Self::NonEmpty(y_coset)) => {
                FinitelyFreeSubmoduleCoset::intersect(&x_coset, &y_coset)
            }
        }
    }

    pub fn equal_slow(x: &Self, y: &Self) -> bool {
        debug_assert_eq!(x.module(), y.module());
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => FinitelyFreeSubmoduleCoset::equal_slow(x_coset, y_coset),
        }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleAffineSubset<Ring> {
    pub fn equal(x: &Self, y: &Self) -> bool {
        debug_assert_eq!(x.module(), y.module());
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => FinitelyFreeSubmoduleCoset::equal(x_coset, y_coset),
        }
    }
}
