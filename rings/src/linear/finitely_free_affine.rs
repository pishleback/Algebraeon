use crate::structure::FinitelyFreeModuleSignature;

use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
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
                Self::NonEmpty(FinitelyFreeSubmoduleCoset::intersect(&x_coset, &y_coset))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
