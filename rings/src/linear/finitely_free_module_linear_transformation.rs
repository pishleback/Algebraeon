use crate::{
    linear::finitely_free_module::FinitelyFreeModuleStructure,
    matrix::{Matrix, MatrixStructure, ReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_structures::*;

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct FreeModuleFiniteNumberedBasisLinearTransformation<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: RingB,
    domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
    range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
    matrix: Matrix<Ring::Elem>, // v -> Mv
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
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
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: BezoutDomainSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
        domain: FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        range: FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Morphism<
        FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
        Ring,
        RingB,
        RingDomainB,
        RingRangeB,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn domain(&self) -> &FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB> {
        &self.domain
    }

    fn range(&self) -> &FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB> {
        &self.range
    }
}

impl<
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
    const SURJECTIVE: bool,
>
    InjectiveFunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
    SetDomain: EnumeratedOrdFiniteSetSignature,
    SetDomainB: BorrowedStructure<SetDomain>,
    SetRange: EnumeratedOrdFiniteSetSignature,
    SetRangeB: BorrowedStructure<SetRange>,
    Ring: ReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    RingDomainB: BorrowedStructure<Ring>,
    RingRangeB: BorrowedStructure<Ring>,
>
    BijectiveFunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, SetDomainB, Ring, RingDomainB>,
        FinitelyFreeModuleStructure<SetRange, SetRangeB, Ring, RingRangeB>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetDomainB,
        SetRange,
        SetRangeB,
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
    use super::*;
    use algebraeon_sets::sets::EnumeratedFiniteSetStructure;

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
