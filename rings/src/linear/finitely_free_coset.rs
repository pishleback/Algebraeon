use crate::structure::ModuleSignature;

use super::{
    finitely_free_modules::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
    matrix::{Matrix, ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Ring: ReducedHermiteAlgorithmSignature> {
    // reduced coset offset
    offset: Vec<Ring::Set>,
    submodule: FinitelyFreeSubmodule<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCoset<Ring> {
    pub fn ring(&self) -> &Ring {
        self.submodule.ring()
    }

    pub fn module(&self) -> FinitelyFreeModuleStructure<Ring> {
        self.submodule.module()
    }

    pub fn module_rank(&self) -> usize {
        self.submodule.module_rank()
    }

    pub fn coset_rank(&self) -> usize {
        self.submodule.submodule_rank()
    }

    pub fn into_offset_and_row_basis_matrix(self) -> (Vec<Ring::Set>, Matrix<Ring::Set>) {
        (self.offset, self.submodule.into_row_basis_matrix())
    }

    pub fn into_offset_and_col_basis_matrix(self) -> (Vec<Ring::Set>, Matrix<Ring::Set>) {
        (self.offset, self.submodule.into_col_basis_matrix())
    }

    pub fn from_offset_and_module(
        offset: &Vec<Ring::Set>,
        submodule: FinitelyFreeSubmodule<Ring>,
    ) -> Self {
        let (_, offset) = submodule.reduce_element(&offset);
        Self { offset, submodule }
    }

    pub fn equal_slow(x: &Self, y: &Self) -> bool {
        let module = common_structure(x.module(), y.module());
        if !FinitelyFreeSubmodule::equal_slow(&x.submodule, &y.submodule) {
            return false;
        }
        x.submodule
            .contains_element(&module.add(&x.offset, &module.neg(&y.offset)))
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCoset<Ring> {
    pub fn equal(x: &Self, y: &Self) -> bool {
        let module = common_structure(x.module(), y.module());
        module.equal(&x.offset, &y.offset)
            && FinitelyFreeSubmodule::equal(&x.submodule, &y.submodule)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_equal_slow() {
        let matrix =
            Matrix::<Integer>::from_rows(vec![vec![1, 2, 3], vec![4, 5, 6], vec![7, 9, 9]]);
        let submodule = matrix.row_span();

        submodule.clone().into_row_basis_matrix().pprint();
    }
}
