use super::matrix::Matrix;
use crate::{linear::matrix::MatrixStructure, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleFiniteNumberedBasisStructure<Ring: RingSignature> {
    ring: Ring,
    rank: usize,
}

impl<Ring: RingSignature> FreeModuleFiniteNumberedBasisStructure<Ring> {
    pub fn new(ring: Ring, rank: usize) -> Self {
        Self { ring, rank }
    }

    pub fn ring(&self) -> &Ring {
        &self.ring
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
                    self.ring.one()
                } else {
                    self.ring.zero()
                }
            })
            .collect()
    }
}

impl<Ring: RingSignature> Signature for FreeModuleFiniteNumberedBasisStructure<Ring> {}

impl<Ring: RingSignature> SetSignature for FreeModuleFiniteNumberedBasisStructure<Ring> {
    type Set = Vec<Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        self.rank == v.len() && v.iter().all(|r| self.ring.is_element(r))
    }
}

impl<Ring: RingSignature> EqSignature for FreeModuleFiniteNumberedBasisStructure<Ring> {
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank).all(|i| self.ring.equal(&v[i], &w[i]))
    }
}

impl<Ring: RingSignature> ModuleSignature<Ring> for FreeModuleFiniteNumberedBasisStructure<Ring> {
    fn ring(&self) -> &Ring {
        &self.ring
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank)
            .map(|i| self.ring.add(&v[i], &w[i]))
            .collect()
    }

    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|r| self.ring.neg(r)).collect()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|s| self.ring.mul(r, s)).collect()
    }
}

impl<Ring: RingSignature> FreeModuleSignature<Ring>
    for FreeModuleFiniteNumberedBasisStructure<Ring>
{
}

impl<Ring: RingSignature> FiniteRankModuleSignature<Ring>
    for FreeModuleFiniteNumberedBasisStructure<Ring>
{
    fn basis(&self) -> Vec<Self::Set> {
        (0..self.rank)
            .into_iter()
            .map(|i| self.basis_element(i))
            .collect()
    }
}

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct FreeModuleFiniteNumberedBasisLinearTransformation<
    Ring: RingSignature,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: Ring,
    domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
    range: FreeModuleFiniteNumberedBasisStructure<Ring>,
    matrix: Matrix<Ring::Set>, // v -> Mv
}

impl<Ring: BezoutDomainSignature, const INJECTIVE: bool, const SURJECTIVE: bool>
    FreeModuleFiniteNumberedBasisLinearTransformation<Ring, INJECTIVE, SURJECTIVE>
{
    pub fn new(
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
        matrix: Matrix<Ring::Set>,
    ) -> Self {
        debug_assert_eq!(&ring, domain.ring());
        debug_assert_eq!(&ring, range.ring());
        debug_assert_eq!(domain.rank(), matrix.cols());
        debug_assert_eq!(range.rank(), matrix.rows());
        let rank = MatrixStructure::new(ring.clone()).rank(matrix.clone());
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
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
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

impl<Ring: BezoutDomainSignature>
    FreeModuleFiniteNumberedBasisLinearTransformation<Ring, false, false>
{
    pub fn construct(
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<Ring: BezoutDomainSignature>
    FreeModuleFiniteNumberedBasisLinearTransformation<Ring, true, false>
{
    pub fn construct_injective(
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<Ring: BezoutDomainSignature>
    FreeModuleFiniteNumberedBasisLinearTransformation<Ring, false, true>
{
    pub fn construct_surjective(
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<Ring: BezoutDomainSignature>
    FreeModuleFiniteNumberedBasisLinearTransformation<Ring, true, true>
{
    pub fn construct_bijective(
        ring: Ring,
        domain: FreeModuleFiniteNumberedBasisStructure<Ring>,
        range: FreeModuleFiniteNumberedBasisStructure<Ring>,
        basis_image: impl Fn(usize) -> Vec<Ring::Set>,
    ) -> Self {
        Self::construct_impl(ring, domain, range, basis_image)
    }
}

impl<Ring: RingSignature, const INJECTIVE: bool, const SURJECTIVE: bool>
    Morphism<
        FreeModuleFiniteNumberedBasisStructure<Ring>,
        FreeModuleFiniteNumberedBasisStructure<Ring>,
    > for FreeModuleFiniteNumberedBasisLinearTransformation<Ring, INJECTIVE, SURJECTIVE>
{
    fn domain(&self) -> &FreeModuleFiniteNumberedBasisStructure<Ring> {
        &self.domain
    }

    fn range(&self) -> &FreeModuleFiniteNumberedBasisStructure<Ring> {
        &self.range
    }
}

impl<Ring: RingSignature, const INJECTIVE: bool, const SURJECTIVE: bool>
    Function<
        FreeModuleFiniteNumberedBasisStructure<Ring>,
        FreeModuleFiniteNumberedBasisStructure<Ring>,
    > for FreeModuleFiniteNumberedBasisLinearTransformation<Ring, INJECTIVE, SURJECTIVE>
{
    fn image(&self, x: &Vec<Ring::Set>) -> Vec<Ring::Set> {
        self.range.from_col(
            &MatrixStructure::new(self.ring.clone())
                .mul(&self.matrix, &self.domain.to_col(x))
                .unwrap(),
        )
    }
}

impl<Ring: BezoutDomainSignature, const SURJECTIVE: bool>
    InjectiveFunction<
        FreeModuleFiniteNumberedBasisStructure<Ring>,
        FreeModuleFiniteNumberedBasisStructure<Ring>,
    > for FreeModuleFiniteNumberedBasisLinearTransformation<Ring, true, SURJECTIVE>
{
    fn try_preimage(&self, y: &Vec<Ring::Set>) -> Option<Vec<Ring::Set>> {
        Some(
            self.domain.from_col(
                &MatrixStructure::new(self.ring.clone())
                    .col_solve(&self.matrix, &self.range.to_col(&y))?,
            ),
        )
    }
}

impl<Ring: BezoutDomainSignature>
    BijectiveFunction<
        FreeModuleFiniteNumberedBasisStructure<Ring>,
        FreeModuleFiniteNumberedBasisStructure<Ring>,
    > for FreeModuleFiniteNumberedBasisLinearTransformation<Ring, true, true>
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
        let m = FreeModuleFiniteNumberedBasisStructure::new(Integer::structure(), 3);

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
        let m = FreeModuleFiniteNumberedBasisStructure::new(Integer::structure(), 2);
        let n = FreeModuleFiniteNumberedBasisStructure::new(Integer::structure(), 5);

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
