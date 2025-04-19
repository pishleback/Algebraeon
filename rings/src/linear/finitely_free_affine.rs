use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Ring: ReducedHermiteAlgorithmSignature> {
    Empty { ring: Ring, module_rank: usize },
    NonEmpty(FinitelyFreeSubmoduleCoset<Ring>),
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleAffineSubset<Ring> {
    pub fn ring(&self) -> &Ring {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty { ring, .. } => ring,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.ring(),
        }
    }

    pub fn module(&self) -> FinitelyFreeModuleStructure<Ring> {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty { module_rank, .. } => {
                FinitelyFreeModuleStructure::new(self.ring().clone(), *module_rank)
            }
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.module(),
        }
    }

    pub fn module_rank(&self) -> usize {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty { module_rank, .. } => *module_rank,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.module_rank(),
        }
    }

    pub fn subset_rank(&self) -> Option<usize> {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty { .. } => None,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => Some(coset.coset_rank()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
