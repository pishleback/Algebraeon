use std::rc::Rc;

use crate::{
    linear::finitely_free_module::FinitelyFreeModuleStructure,
    matrix::{Matrix, MatrixStructure, ReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_structures::*;

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct FreeModuleFiniteNumberedBasisLinearTransformation<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: RingSignature,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: Rc<Ring>,
    domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
    range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
    matrix: Matrix<Ring::Elem>, // v -> Mv
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: BezoutDomainSignature,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetRange,
        Ring,
        INJECTIVE,
        SURJECTIVE,
    >
{
    pub fn new(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        matrix: Matrix<Ring::Elem>,
    ) -> Rc<Self> {
        debug_assert_eq!(ring, domain.ring());
        debug_assert_eq!(ring, range.ring());
        debug_assert_eq!(domain.rank(), matrix.cols());
        debug_assert_eq!(range.rank(), matrix.rows());
        #[cfg(debug_assertions)]
        {
            let rank = MatrixStructure::<Ring>::new(ring.clone()).rank(matrix.clone());
            if INJECTIVE {
                assert_eq!(rank, domain.rank());
            }
            if SURJECTIVE {
                assert_eq!(rank, range.rank());
            }
        }
        Self {
            ring,
            domain,
            range,
            matrix,
        }
        .into()
    }

    fn construct_impl(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Rc<Self> {
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
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: BezoutDomainSignature,
> FreeModuleFiniteNumberedBasisLinearTransformation<SetDomain, SetRange, Ring, false, false>
{
    pub fn construct(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Rc<Self> {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: BezoutDomainSignature,
> FreeModuleFiniteNumberedBasisLinearTransformation<SetDomain, SetRange, Ring, true, false>
{
    pub fn construct_injective(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Rc<Self> {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: BezoutDomainSignature,
> FreeModuleFiniteNumberedBasisLinearTransformation<SetDomain, SetRange, Ring, false, true>
{
    pub fn construct_surjective(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Rc<Self> {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: BezoutDomainSignature,
> FreeModuleFiniteNumberedBasisLinearTransformation<SetDomain, SetRange, Ring, true, true>
{
    pub fn construct_bijective(
        ring: Rc<Ring>,
        domain: Rc<FinitelyFreeModuleStructure<SetDomain, Ring>>,
        range: Rc<FinitelyFreeModuleStructure<SetRange, Ring>>,
        basis_image: impl Fn(usize) -> Vec<Ring::Elem>,
    ) -> Rc<Self> {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: RingSignature,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    Morphism<
        FinitelyFreeModuleStructure<SetDomain, Ring>,
        FinitelyFreeModuleStructure<SetRange, Ring>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetRange,
        Ring,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn domain(self: &Rc<Self>) -> Rc<FinitelyFreeModuleStructure<SetDomain, Ring>> {
        self.domain.clone()
    }

    fn range(self: &Rc<Self>) -> Rc<FinitelyFreeModuleStructure<SetRange, Ring>> {
        self.range.clone()
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: RingSignature,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
>
    FunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, Ring>,
        FinitelyFreeModuleStructure<SetRange, Ring>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetRange,
        Ring,
        INJECTIVE,
        SURJECTIVE,
    >
{
    fn image(self: &Rc<Self>, x: &Vec<Ring::Elem>) -> Vec<Ring::Elem> {
        self.range.from_col(
            &MatrixStructure::new(self.ring.clone())
                .mul(&self.matrix, &self.domain.to_col(x))
                .unwrap(),
        )
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    const SURJECTIVE: bool,
>
    InjectiveFunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, Ring>,
        FinitelyFreeModuleStructure<SetRange, Ring>,
    >
    for FreeModuleFiniteNumberedBasisLinearTransformation<
        SetDomain,
        SetRange,
        Ring,
        true,
        SURJECTIVE,
    >
{
    fn try_preimage(self: &Rc<Self>, y: &Vec<Ring::Elem>) -> Option<Vec<Ring::Elem>> {
        MatrixStructure::new(self.ring.clone()).col_solve(self.matrix.clone(), y)
    }
}

impl<
    SetDomain: OrderedFiniteSetSignature,
    SetRange: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
>
    BijectiveFunctionMorphism<
        FinitelyFreeModuleStructure<SetDomain, Ring>,
        FinitelyFreeModuleStructure<SetRange, Ring>,
    > for FreeModuleFiniteNumberedBasisLinearTransformation<SetDomain, SetRange, Ring, true, true>
{
    fn preimage(self: &Rc<Self>, y: &Vec<Ring::Elem>) -> Vec<Ring::Elem> {
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
