use std::marker::PhantomData;

use super::{
    finitely_free_affine::FinitelyFreeSubmoduleAffineSubsetStructure,
    finitely_free_coset::FinitelyFreeSubmoduleCosetStructure,
    finitely_free_submodule::{FinitelyFreeSubmodule, FinitelyFreeSubmoduleStructure},
    matrix::{Matrix, ReducedHermiteAlgorithmSignature},
};
use crate::{linear::matrix::MatrixStructure, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeModuleStructure<Ring: RingSignature, RingB: BorrowedStructure<Ring>> {
    _ring: PhantomData<Ring>,
    ring: RingB,
    rank: usize,
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> FinitelyFreeModuleStructure<Ring, RingB> {
    pub fn new(ring: RingB, rank: usize) -> Self {
        Self {
            _ring: PhantomData::default(),
            ring,
            rank,
        }
    }
}
pub trait RingToFinitelyFreeModuleStructure: RingSignature {
    fn free_module<'a>(&'a self, n: usize) -> FinitelyFreeModuleStructure<Self, &'a Self> {
        FinitelyFreeModuleStructure::new(self, n)
    }
    fn into_free_module(self, n: usize) -> FinitelyFreeModuleStructure<Self, Self> {
        FinitelyFreeModuleStructure::new(self, n)
    }
}
impl<Ring: RingSignature> RingToFinitelyFreeModuleStructure for Ring {}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> FinitelyFreeModuleStructure<Ring, RingB> {
    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    pub fn to_col(&self, v: &<Self as SetSignature>::Set) -> Matrix<Ring::Set> {
        debug_assert!(self.is_element(&v));
        Matrix::construct(self.rank, 1, |r, _| v[r].clone())
    }

    pub fn to_row(&self, v: &<Self as SetSignature>::Set) -> Matrix<Ring::Set> {
        debug_assert!(self.is_element(&v));
        Matrix::construct(1, self.rank, |_, c| v[c].clone())
    }

    pub fn from_row(&self, m: &Matrix<Ring::Set>) -> <Self as SetSignature>::Set {
        debug_assert_eq!(m.rows(), 1);
        debug_assert_eq!(m.cols(), self.rank);
        (0..self.rank)
            .map(|i| m.at(0, i).unwrap().clone())
            .collect()
    }

    pub fn from_col(&self, m: &Matrix<Ring::Set>) -> <Self as SetSignature>::Set {
        debug_assert_eq!(m.cols(), 1);
        debug_assert_eq!(m.rows(), self.rank);
        (0..self.rank)
            .map(|i| m.at(i, 0).unwrap().clone())
            .collect()
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Set {
        debug_assert!(i < self.rank);
        (0..self.rank)
            .map(|j| {
                if i == j {
                    self.ring().one()
                } else {
                    self.ring().zero()
                }
            })
            .collect()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, RingB: BorrowedStructure<Ring>>
    FinitelyFreeModuleStructure<Ring, RingB>
{
    pub fn submodules<'a>(&'a self) -> FinitelyFreeSubmoduleStructure<Ring, &'a Ring> {
        FinitelyFreeSubmoduleStructure::new(FinitelyFreeModuleStructure::new(
            self.ring(),
            self.rank(),
        ))
    }

    pub fn into_submodules(self) -> FinitelyFreeSubmoduleStructure<Ring, RingB> {
        FinitelyFreeSubmoduleStructure::new(self)
    }

    pub fn cosets<'a>(&'a self) -> FinitelyFreeSubmoduleCosetStructure<Ring, &'a Ring> {
        FinitelyFreeSubmoduleCosetStructure::new(FinitelyFreeModuleStructure::new(
            self.ring(),
            self.rank(),
        ))
    }

    pub fn into_cosets(self) -> FinitelyFreeSubmoduleCosetStructure<Ring, RingB> {
        FinitelyFreeSubmoduleCosetStructure::new(self)
    }

    pub fn affine_subsets<'a>(
        &'a self,
    ) -> FinitelyFreeSubmoduleAffineSubsetStructure<Ring, &'a Ring> {
        FinitelyFreeSubmoduleAffineSubsetStructure::new(FinitelyFreeModuleStructure::new(
            self.ring(),
            self.rank(),
        ))
    }

    pub fn into_affine_subsets(self) -> FinitelyFreeSubmoduleAffineSubsetStructure<Ring, RingB> {
        FinitelyFreeSubmoduleAffineSubsetStructure::new(self)
    }

    pub fn improper_submodule(&self) -> FinitelyFreeSubmodule<Ring::Set> {
        self.submodules()
            .matrix_row_span(MatrixStructure::new(self.ring().clone()).ident(self.rank()))
    }

    pub fn generated_submodule(
        &self,
        generators: Vec<&Vec<Ring::Set>>,
    ) -> FinitelyFreeSubmodule<Ring::Set> {
        for generator in &generators {
            debug_assert!(self.is_element(generator));
        }
        let row_span = Matrix::construct(generators.len(), self.rank(), |r, c| {
            generators[r][c].clone()
        });
        self.submodules().matrix_row_span(row_span)
    }
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> Signature
    for FinitelyFreeModuleStructure<Ring, RingB>
{
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> SetSignature
    for FinitelyFreeModuleStructure<Ring, RingB>
{
    type Set = Vec<Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        self.rank == v.len() && v.iter().all(|r| self.ring().is_element(r))
    }
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> EqSignature
    for FinitelyFreeModuleStructure<Ring, RingB>
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank).all(|i| self.ring().equal(&v[i], &w[i]))
    }
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> ModuleSignature<Ring>
    for FinitelyFreeModuleStructure<Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    fn zero(&self) -> Self::Set {
        (0..self.rank).map(|_| self.ring().zero()).collect()
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank)
            .map(|i| self.ring().add(&v[i], &w[i]))
            .collect()
    }

    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|r| self.ring().neg(r)).collect()
    }

    fn sub(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank)
            .map(|i| self.ring().sub(&v[i], &w[i]))
            .collect()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|s| self.ring().mul(r, s)).collect()
    }
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> FreeModuleSignature<Ring>
    for FinitelyFreeModuleStructure<Ring, RingB>
{
}

impl<Ring: RingSignature, RingB: BorrowedStructure<Ring>> FinitelyFreeModuleSignature<Ring>
    for FinitelyFreeModuleStructure<Ring, RingB>
{
    fn rank(&self) -> usize {
        self.rank
    }

    fn to_vec(&self, v: &Self::Set) -> Vec<Ring::Set> {
        v.clone()
    }

    fn from_vec(&self, v: &Vec<Ring::Set>) -> Self::Set {
        v.clone()
    }
}

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct FreeModuleFiniteNumberedBasisLinearTransformation<
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: RingB,
    domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
    range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
    matrix: Matrix<Ring::Set>, // v -> Mv
}

impl<
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    pub fn new(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        matrix: Matrix<Ring::Set>,
    ) -> Self {
        debug_assert_eq!(ring.borrow(), domain.ring());
        debug_assert_eq!(ring.borrow(), range.ring());
        debug_assert_eq!(domain.rank(), matrix.cols());
        debug_assert_eq!(range.rank(), matrix.rows());
        let rank = MatrixStructure::<Ring, _>::new(ring.borrow()).rank(matrix.clone());
        if INJECTIVE {
            debug_assert_eq!(rank, domain.rank());
        }
        if SURJECTIVE {
            debug_assert_eq!(rank, range.rank());
        }
        Self {
            ring,
            domain,
            range,
            matrix,
        }
    }

    fn construct_impl(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        let matrix = Matrix::from_cols(
            (0..domain.rank())
                .map(|i| {
                    let img_i = basis_image(i);
                    debug_assert!(range.is_element(&img_i));
                    img_i
                })
                .collect(),
        );
        Self::new(ring, domain, range, matrix)
    }
}

impl<
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        false,
        false,
    >
{
    pub fn construct(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        false,
    >
{
    pub fn construct_injective(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        false,
        true,
    >
{
    pub fn construct_surjective(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        true,
    >
{
    pub fn construct_bijective(
        ring: RingB,
        domain: FinitelyFreeModuleStructure<Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Morphism<
        FinitelyFreeModuleStructure<Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn domain(&self) -> &FinitelyFreeModuleStructure<Ring, RingDomainB> {
        &self.domain
    }

    fn range(&self) -> &FinitelyFreeModuleStructure<Ring, RingRangeB> {
        &self.range
    }
}

impl<
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Function<
        FinitelyFreeModuleStructure<Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn image(&self, x: &Vec<Ring::Set>) -> Vec<Ring::Set> {
        self.range.from_col(
            &MatrixStructure::new(self.ring.clone())
                .mul(&self.matrix, &self.domain.to_col(x))
                .unwrap(),
        )
    }
}

impl<
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const SURJECTIVE: bool,
>
    InjectiveFunction<
        FinitelyFreeModuleStructure<Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        SURJECTIVE,
    >
{
    fn try_preimage(&self, y: &Vec<Ring::Set>) -> Option<Vec<Ring::Set>> {
        MatrixStructure::new(self.ring.clone()).col_solve(self.matrix.clone(), &y)
    }
}

impl<
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    BijectiveFunction<
        FinitelyFreeModuleStructure<Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        true,
    >
{
    fn preimage(&self, y: &Vec<Ring::Set>) -> Vec<Ring::Set> {
        self.try_preimage(y).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::{FreeModuleFiniteNumberedBasisLinearTransformation, *};
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finite_rank_modules() {
        let m = FinitelyFreeModuleStructure::new(Integer::structure(), 3);

        let a = m.basis_element(0);
        let b = m.basis_element(1);
        let c = m.basis_element(2);

        assert_eq!(
            m.add(&m.neg(&b), &m.add(&a, &b)),
            vec![Integer::from(1), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            vec![Integer::from(1), Integer::from(2), Integer::from(1)]
        );

        assert_eq!(
            m.scalar_mul(&5.into(), &a),
            vec![Integer::from(5), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(m.basis(), vec![a, b, c]);
    }

    #[test]
    fn test_finite_rank_modules_linear_transformation() {
        let m = FinitelyFreeModuleStructure::new(Integer::structure(), 2);
        let n = FinitelyFreeModuleStructure::new(Integer::structure(), 5);

        let t = FreeModuleFiniteNumberedBasisLinearTransformation::construct_injective(
            Integer::structure(),
            m.clone(),
            n.clone(),
            |i| {
                if i == 0 {
                    vec![0, 2, 3, -4, 1]
                        .into_iter()
                        .map(|x| Integer::from(x))
                        .collect()
                } else if i == 1 {
                    vec![1, 2, 3, 2, 1]
                        .into_iter()
                        .map(|x| Integer::from(x))
                        .collect()
                } else {
                    unreachable!()
                }
            },
        );

        assert_eq!(
            t.image(&vec![Integer::from(1), Integer::from(2)]),
            vec![2, 6, 9, 0, 3]
                .into_iter()
                .map(|x| Integer::from(x))
                .collect::<Vec<_>>()
        );
    }
}
