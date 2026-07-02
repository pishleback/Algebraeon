use super::{
    finitely_free_affine::FinitelyFreeSubmoduleAffineSubsetStructure,
    finitely_free_coset::FinitelyFreeSubmoduleCosetStructure,
    finitely_free_submodule::{FinitelyFreeSubmodule, FinitelyFreeSubmoduleStructure},
};
use crate::{
    matrix::{Matrix, MatrixStructure, ReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_structures::*;
use std::{borrow::Cow, marker::PhantomData};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeModuleStructure<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _set: PhantomData<Set>,
    set: SetB,
    _ring: PhantomData<Ring>,
    ring: RingB,
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    pub fn new(set: SetB, ring: RingB) -> Self {
        Self {
            _set: PhantomData,
            set,
            _ring: PhantomData,
            ring,
        }
    }
}
pub trait RingToFinitelyFreeModuleSignature: RingSignature {
    fn free_module<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>(
        &self,
        set: SetB,
    ) -> FinitelyFreeModuleStructure<Set, SetB, Self, &Self> {
        FinitelyFreeModuleStructure::new(set, self)
    }

    fn into_free_module<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>(
        self,
        set: SetB,
    ) -> FinitelyFreeModuleStructure<Set, SetB, Self, Self> {
        FinitelyFreeModuleStructure::new(set, self)
    }
}
impl<Ring: RingSignature> RingToFinitelyFreeModuleSignature for Ring {}

pub trait SetToFinitelyFreeModuleSignature: EnumeratedOrdFiniteSetSignature {
    fn free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        &self,
        ring: RingB,
    ) -> FinitelyFreeModuleStructure<Self, &Self, Ring, RingB> {
        FinitelyFreeModuleStructure::new(self, ring)
    }

    fn into_free_module<Ring: RingSignature, RingB: BorrowedStructure<Ring>>(
        self,
        ring: RingB,
    ) -> FinitelyFreeModuleStructure<Self, Self, Ring, RingB> {
        FinitelyFreeModuleStructure::new(self, ring)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature> SetToFinitelyFreeModuleSignature for Set {}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    pub fn set(&self) -> &Set {
        self.set.borrow()
    }

    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    pub fn to_col(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(self.rank(), 1, |r, _| v[r].clone())
    }

    pub fn to_row(&self, v: &<Self as SetSignature>::Elem) -> Matrix<Ring::Elem> {
        debug_assert!(self.validate_element(v).is_ok());
        Matrix::construct(1, self.rank(), |_, c| v[c].clone())
    }

    pub fn from_row(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.rows(), 1);
        debug_assert_eq!(m.cols(), self.rank());
        (0..self.rank())
            .map(|i| m.at(0, i).unwrap().clone())
            .collect()
    }

    pub fn from_col(&self, m: &Matrix<Ring::Elem>) -> <Self as SetSignature>::Elem {
        debug_assert_eq!(m.cols(), 1);
        debug_assert_eq!(m.rows(), self.rank());
        (0..self.rank())
            .map(|i| m.at(i, 0).unwrap().clone())
            .collect()
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Elem {
        debug_assert!(i < self.rank());
        (0..self.rank())
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

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    pub fn submodules(&self) -> FinitelyFreeSubmoduleStructure<Set, &Set, Ring, &Ring> {
        FinitelyFreeSubmoduleStructure::new(FinitelyFreeModuleStructure::new(
            self.set(),
            self.ring(),
        ))
    }

    pub fn into_submodules(self) -> FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB> {
        FinitelyFreeSubmoduleStructure::new(self)
    }

    pub fn cosets(&self) -> FinitelyFreeSubmoduleCosetStructure<Set, &Set, Ring, &Ring> {
        FinitelyFreeSubmoduleCosetStructure::new(FinitelyFreeModuleStructure::new(
            self.set(),
            self.ring(),
        ))
    }

    pub fn into_cosets(self) -> FinitelyFreeSubmoduleCosetStructure<Set, SetB, Ring, RingB> {
        FinitelyFreeSubmoduleCosetStructure::new(self)
    }

    pub fn affine_subsets(
        &self,
    ) -> FinitelyFreeSubmoduleAffineSubsetStructure<Set, &Set, Ring, &Ring> {
        FinitelyFreeSubmoduleAffineSubsetStructure::new(FinitelyFreeModuleStructure::new(
            self.set(),
            self.ring(),
        ))
    }

    pub fn into_affine_subsets(
        self,
    ) -> FinitelyFreeSubmoduleAffineSubsetStructure<Set, SetB, Ring, RingB> {
        FinitelyFreeSubmoduleAffineSubsetStructure::new(self)
    }

    pub fn improper_submodule(&self) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.submodules()
            .matrix_row_span(MatrixStructure::new(self.ring().clone()).ident(self.rank()))
    }

    pub fn generated_submodule(
        &self,
        generators: Vec<&Vec<Ring::Elem>>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        for generator in &generators {
            debug_assert!(self.validate_element(generator).is_ok());
        }
        let row_span = Matrix::construct(generators.len(), self.rank(), |r, c| {
            generators[r][c].clone()
        });
        self.submodules().matrix_row_span(row_span)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> Signature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SetSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    type Elem = Vec<Ring::Elem>;

    fn validate_element(&self, v: &Self::Elem) -> Result<(), String> {
        if self.rank() != v.len() {
            return Err("wrong size".to_string());
        }
        for r in v {
            self.ring().validate_element(r)?;
        }
        Ok(())
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature + EqSignature,
    RingB: BorrowedStructure<Ring>,
> EqSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn equal(&self, v: &Self::Elem, w: &Self::Elem) -> bool {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank()).all(|i| self.ring().equal(&v[i], &w[i]))
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> RinglikeSpecializationSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> ZeroSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn zero(&self) -> Self::Elem {
        (0..self.rank()).map(|_| self.ring().zero()).collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditionSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn add(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().add(&v[i], &w[i]))
            .collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> CancellativeAdditionSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> TryNegateSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveMonoidSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveGroupSignature for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn neg(&self, v: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|r| self.ring().neg(r)).collect()
    }

    fn sub(&self, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        (0..self.rank())
            .map(|i| self.ring().sub(&v[i], &w[i]))
            .collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> SemiModuleSignature<Ring> for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    fn scalar_mul(&self, v: &Self::Elem, r: &Ring::Elem) -> Self::Elem {
        debug_assert!(self.validate_element(v).is_ok());
        v.iter().map(|s| self.ring().mul(r, s)).collect()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleSignature<Ring> for FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>
{
    type Basis = Set;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        self.set()
    }

    fn to_component<'a>(&self, b: &Set::Elem, v: &'a Self::Elem) -> Cow<'a, Ring::Elem> {
        let b: usize = self.set().element_to_enumeration(b).try_into().unwrap();
        debug_assert!(b < self.rank());
        Cow::Borrowed(&v[b])
    }

    fn from_component(&self, b: &Set::Elem, r: &<Ring>::Elem) -> Self::Elem {
        let b: usize = self.set().element_to_enumeration(b).try_into().unwrap();
        debug_assert!(b < self.rank());
        let mut element = self.zero();
        element[b] = r.clone();
        element
    }
}

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct FreeModuleFiniteNumberedBasisLinearTransformation<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: RingB,
    domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
    range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
    matrix: Matrix<Ring::Elem>, // v -> Mv
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        matrix: Matrix<Ring::Elem>,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        let matrix = Matrix::from_cols(
            (0..domain.rank())
                .map(|i| {
                    let img_i = basis_image(i);
                    debug_assert!(range.validate_element(&img_i).is_ok());
                    img_i
                })
                .collect(),
        );
        Self::new(ring, domain, range, matrix)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
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
        domain: FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Morphism<
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn domain(&self) -> &FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB> {
        &self.domain
    }

    fn range(&self) -> &FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB> {
        &self.range
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Function<
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn image(&self, x: &Vec<Ring::Elem>) -> Vec<Ring::Elem> {
        self.range.from_col(
            &MatrixStructure::new(self.ring.clone())
                .mul(&self.matrix, &self.domain.to_col(x))
                .unwrap(),
        )
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const SURJECTIVE: bool,
>
    InjectiveFunction<
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        SURJECTIVE,
    >
{
    fn try_preimage(&self, y: &Vec<Ring::Elem>) -> Option<Vec<Ring::Elem>> {
        MatrixStructure::new(self.ring.clone()).col_solve(self.matrix.clone(), y)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    BijectiveFunction<
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<Set, SetB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        Set,
        SetB,
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        true,
        true,
    >
{
    fn preimage(&self, y: &Vec<Ring::Elem>) -> Vec<Ring::Elem> {
        self.try_preimage(y).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::{FreeModuleFiniteNumberedBasisLinearTransformation, *};
    use algebraeon_sets::sets::EnumeratedFiniteSetStructure;

    #[test]
    fn test_finite_rank_modules() {
        let m = FinitelyFreeModuleStructure::new(
            EnumeratedFiniteSetStructure::new(3),
            Integer::structure(),
        );

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
            m.scalar_mul(&a, &5.into()),
            vec![Integer::from(5), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(m.basis_vecs(), vec![a, b, c]);
    }

    #[test]
    fn test_finite_rank_modules_linear_transformation() {
        let m = FinitelyFreeModuleStructure::new(
            EnumeratedFiniteSetStructure::new(2),
            Integer::structure(),
        );
        let n = FinitelyFreeModuleStructure::new(
            EnumeratedFiniteSetStructure::new(5),
            Integer::structure(),
        );

        let t = FreeModuleFiniteNumberedBasisLinearTransformation::construct_injective(
            Integer::structure(),
            m.clone(),
            n.clone(),
            |i| {
                if i == 0 {
                    vec![0, 2, 3, -4, 1]
                        .into_iter()
                        .map(Integer::from)
                        .collect()
                } else if i == 1 {
                    vec![1, 2, 3, 2, 1].into_iter().map(Integer::from).collect()
                } else {
                    unreachable!()
                }
            },
        );

        assert_eq!(
            t.image(&vec![Integer::from(1), Integer::from(2)]),
            vec![2, 6, 9, 0, 3]
                .into_iter()
                .map(Integer::from)
                .collect::<Vec<_>>()
        );
    }
}
