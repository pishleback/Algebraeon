use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
use crate::{linear::matrix::*, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmodule<Ring: RingSignature> {
    // a matrix in hermite normal form with all non-zero rows whose rows form a basis for the submodule
    // may also be required to be reduced or unique depending on the structure used
    row_basis: Matrix<Ring::Set>,
    // the columns of the pivots of row_basis
    pivots: Vec<usize>,
}

impl<Ring: RingSignature> FinitelyFreeSubmodule<Ring> {
    pub fn into_row_basis(self) -> Matrix<Ring::Set> {
        self.row_basis
    }

    pub fn into_col_basis(self) -> Matrix<Ring::Set> {
        self.row_basis.transpose()
    }
}

/// The set of all submodules of some module represented by a basis in hermite normal form
/// REDUCED: Submodules are represented by matricies in _reduced_ hermite normal form
/// UNIQUE: The hermite normal form representing a submodule is _unique_
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleStructure<
    Ring: RingSignature,
    const REDUCED: bool,
    const UNIQUE: bool,
> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool>
    FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleAlgorithms<Ring>,
{
    pub fn basis(&self, sm: &FinitelyFreeSubmodule<Ring>) -> Vec<Vec<Ring::Set>> {
        debug_assert!(self.is_element(sm));
        (0..sm.row_basis.rows())
            .map(|r| {
                (0..self.module().rank())
                    .map(|c| sm.row_basis.at(r, c).unwrap().clone())
                    .collect()
            })
            .collect()
    }
}

impl<Ring: HermiteAlgorithmSignature> FinitelyFreeSubmoduleStructure<Ring, false, false> {
    pub fn new_unreduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleStructure<Ring, true, false> {
    pub fn new_reduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleStructure<Ring, true, true>
{
    pub fn new_uniquely_reduced(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

pub trait FinitelyFreeSubmoduleAlgorithms<Ring: RingSignature> {
    fn matrix_row_span(&self, row_span: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring>;

    fn matrix_col_span(&self, col_span: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring> {
        self.matrix_row_span(col_span.transpose())
    }
}

impl<Ring: HermiteAlgorithmSignature, const UNIQUE: bool> FinitelyFreeSubmoduleAlgorithms<Ring>
    for FinitelyFreeSubmoduleStructure<Ring, false, UNIQUE>
{
    fn matrix_row_span(&self, row_span: Matrix<<Ring>::Set>) -> FinitelyFreeSubmodule<Ring> {
        let mat_ring = MatrixStructure::new(self.ring().clone());
        let (h, _u, _det, pivots) = mat_ring.row_hermite_algorithm(row_span);
        let row_basis = h.submatrix(
            (0..pivots.len()).collect(),
            (0..self.module().rank()).collect(),
        );
        FinitelyFreeSubmodule { row_basis, pivots }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleAlgorithms<Ring> for FinitelyFreeSubmoduleStructure<Ring, true, UNIQUE>
{
    fn matrix_row_span(&self, row_span: Matrix<Ring::Set>) -> FinitelyFreeSubmodule<Ring> {
        let mat_ring = MatrixStructure::new(self.ring().clone());
        debug_assert_eq!(row_span.cols(), self.module().rank());
        let (h, _u, _det, pivots) = mat_ring.row_reduced_hermite_algorithm(row_span);
        let row_basis = h.submatrix(
            (0..pivots.len()).collect(),
            (0..self.module().rank()).collect(),
        );
        FinitelyFreeSubmodule { row_basis, pivots }
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> Signature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
{
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleAlgorithms<Ring>,
{
    type Set = FinitelyFreeSubmodule<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        if sm.row_basis.cols() != self.module().rank() {
            return false;
        }
        // todo check sm.row_basis is in hermite normal form with all non-zero rows
        if REDUCED {
            // todo check sm.row_basis is in _reduced_ hermite normal form
        }
        true
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleAlgorithms<Ring>,
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        self.contains(&x, &y) && self.contains(&y, &x)
    }
}

impl<Ring: RingSignature, const REDUCED: bool, const UNIQUE: bool>
    SubModuleSignature<Ring, FinitelyFreeModuleStructure<Ring>>
    for FinitelyFreeSubmoduleStructure<Ring, REDUCED, UNIQUE>
where
    Self: FinitelyFreeSubmoduleAlgorithms<Ring>,
{
    fn ring(&self) -> &Ring {
        self.module.ring()
    }

    fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }

    fn improper_submodule(&self) -> Self::Set {
        todo!()
    }

    fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        self.matrix_row_span(Matrix::join_rows(
            self.module().rank(),
            vec![&x.row_basis, &y.row_basis],
        ))
    }

    fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        todo!()
    }

    fn generated(&self, generators: Vec<&Vec<Ring::Set>>) -> Self::Set {
        for generator in &generators {
            debug_assert!(self.module.is_element(generator));
        }
        let row_span = Matrix::construct(generators.len(), self.module().rank(), |r, c| {
            generators[r][c].clone()
        });
        self.matrix_row_span(row_span)
    }

    fn contains_element(&self, x: &Self::Set, p: &Vec<Ring::Set>) -> bool {
        debug_assert!(self.is_element(x));
        debug_assert!(self.module().is_element(&p));
        todo!()
    }

    fn contains(&self, x: &Self::Set, y: &Self::Set) -> bool {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        for g in self.basis(y) {
            if !self.contains_element(x, &g) {
                return false;
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finitely_free_submodule_reduction() {
        let module = FinitelyFreeModuleStructure::new(Integer::structure(), 4);
        {
            let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
            a.pprint();
            let a_reduced = FinitelyFreeSubmoduleStructure::new_uniquely_reduced(module.clone())
                .matrix_row_span(a)
                .into_row_basis();
            let a_expected = Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]);
            a_reduced.pprint();
            a_expected.pprint();
            assert_eq!(a_reduced, a_expected);
        }
        println!();
        {
            let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
            a.pprint();
            let a_reduced = FinitelyFreeSubmoduleStructure::new_unreduced(module.clone())
                .matrix_row_span(a)
                .into_row_basis();
            let a_expected = Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]]);
            a_reduced.pprint();
            a_expected.pprint();
            assert_eq!(a_reduced, a_expected);
        }
    }
}

// use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
// use crate::structure::*;
// use algebraeon_nzq::IntegerCanonicalStructure;
// use algebraeon_sets::structure::*;

// #[derive(Debug, Clone)]
// pub struct SubmoduleOfFinitelyFreeModule<Ring: RingSignature> {
//     // rows are a basis for the submodule
//     basis: Matrix<Ring::Set>,
// }

// pub trait SubmoduleOfFinitelyFreeModuleSignature<Ring: RingSignature>:
//     SetSignature<Set = SubmoduleOfFinitelyFreeModule<Ring>>
// {
//     fn module(&self) -> &FinitelyFreeModuleStructure<Ring>;
//     fn ring(&self) -> &Ring {
//         self.module().ring()
//     }
// }

// pub trait SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<Ring: RingSignature>:
//     SubmoduleOfFinitelyFreeModuleSignature<Ring>
// {
//     fn unique_reduce(&self, basis: Self::Set) -> Self::Set;
// }

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct SubmoduleOfFinitelyFreeModuleStructure<Ring: BezoutDomainSignature> {
//     module: FinitelyFreeModuleStructure<Ring>,
// }

// impl<Ring: BezoutDomainSignature> Signature for SubmoduleOfFinitelyFreeModuleStructure<Ring> {}

// impl<Ring: BezoutDomainSignature> SetSignature for SubmoduleOfFinitelyFreeModuleStructure<Ring> {
//     type Set = SubmoduleOfFinitelyFreeModule<Ring>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         self.module.rank() == x.basis.cols()
//     }
// }

// impl<Ring: BezoutDomainSignature> SubmoduleOfFinitelyFreeModuleSignature<Ring>
//     for SubmoduleOfFinitelyFreeModuleStructure<Ring>
// {
//     fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
//         &self.module
//     }
// }

// impl<Field: FieldSignature> SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<Field>
//     for SubmoduleOfFinitelyFreeModuleStructure<Field>
// {
//     fn unique_reduce(&self, basis: Self::Set) -> Self::Set {
//         debug_assert!(self.is_element(&basis));
//         todo!()
//     }
// }

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct IntegralSubmoduleOfFinitelyFreeModuleStructure {
//     module: FinitelyFreeModuleStructure<IntegerCanonicalStructure>,
// }

// impl Signature for IntegralSubmoduleOfFinitelyFreeModuleStructure {}

// impl SetSignature for IntegralSubmoduleOfFinitelyFreeModuleStructure {
//     type Set = SubmoduleOfFinitelyFreeModule<IntegerCanonicalStructure>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         self.module.rank() == x.basis.cols()
//     }
// }

// impl SubmoduleOfFinitelyFreeModuleSignature<IntegerCanonicalStructure>
//     for IntegralSubmoduleOfFinitelyFreeModuleStructure
// {
//     fn module(&self) -> &FinitelyFreeModuleStructure<IntegerCanonicalStructure> {
//         &self.module
//     }
// }

// impl SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<IntegerCanonicalStructure>
//     for IntegralSubmoduleOfFinitelyFreeModuleStructure
// {
//     fn unique_reduce(&self, basis: Self::Set) -> Self::Set {
//         debug_assert!(self.is_element(&basis));
//         todo!()
//     }
// }

// #[derive(Debug, Clone)]
// pub struct FinitelyFreeSubmoduleStructure<
//     Ring: EqSignature + BezoutDomainSignature,
//     Module: FinitelyFreeModuleSignature<Ring>,
// > {
//     ring: PhantomData<Ring>,
//     module: Module,
//     // linearly independent rows
//     basis_matrix: Matrix<Ring::Set>,
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     pub fn ring(&self) -> &Ring {
//         self.module.ring()
//     }

//     pub fn module(&self) -> &Module {
//         &self.module
//     }

//     pub fn from_span(module: Module, span: Vec<Module::Set>) -> Self {
//         for v in span {
//             debug_assert!(module.is_element(&v));
//         }
//         module
//         todo!()
//     }

// pub fn from_basis(module: Module, basis: Vec<Module::Set>) -> Self {
//     for v in basis {
//         debug_assert!(module.is_element(&v));
//     }
//     todo!()
// }

// pub fn submodule(
//     &self,
// ) -> impl InjectiveFunction<FinitelyFreeModuleStructure<Ring>, Module>
// + LinearTransformation<Ring, FinitelyFreeModuleStructure<Ring>, Module> {
//     todo!()
// }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>> Signature
//     for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     SetSignature for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     type Set = Vec<Ring::Set>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         self.submodule().is
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     EqSignature for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         todo!()
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     ModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn ring(&self) -> &Ring {
//         todo!()
//     }

//     fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
//         todo!()
//     }

//     fn neg(&self, a: &Self::Set) -> Self::Set {
//         todo!()
//     }

//     fn scalar_mul(&self, x: &<Ring>::Set, a: &Self::Set) -> Self::Set {
//         todo!()
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FreeModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FinitelyFreeModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn rank(&self) -> usize {
//         self.basis_matrix.rows()
//     }

//     fn basis(&self) -> Vec<Self::Set> {
//         self.submodule().basis()
//     }
// }
