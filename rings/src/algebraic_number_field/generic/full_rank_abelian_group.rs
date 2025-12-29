use crate::{
    algebraic_number_field::AlgebraicNumberFieldSignature,
    matrix::Matrix,
    module::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, FiniteDimensionalFieldExtension,
        FiniteRankFreeRingExtension, MetaCharZeroRing,
    },
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Rational};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, Function, InjectiveFunction, MetaType, Morphism, SetSignature,
    Signature,
};
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
> {
    anf: KB,
    // length = degree of k
    basis: Vec<K::Set>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    pub fn new(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let n = anf.borrow().n();
        if n != basis.len() {
            return Err("Basis has wrong length".to_string());
        }
        for v in &basis {
            if let Err(e) = anf.borrow().is_element(v) {
                return Err(format!(
                    "Vector is not a valid element of the number field: {}",
                    e
                ));
            }
        }
        let mat = Matrix::join_cols(
            n,
            (0..n)
                .map(|i| {
                    anf.borrow()
                        .finite_dimensional_rational_extension()
                        .to_col(&basis[i])
                })
                .collect(),
        );
        if mat.rank() != n {
            return Err("Vectors do not form a basis".to_string());
        }
        Ok(Self { anf, basis })
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let n = anf.borrow().n();
        #[cfg(debug_assertions)]
        {
            assert_eq!(n, basis.len());
            for v in &basis {
                assert!(anf.borrow().is_element(v).is_ok());
            }
            let mat = Matrix::join_cols(
                n,
                (0..n)
                    .map(|i| {
                        anf.borrow()
                            .finite_dimensional_rational_extension()
                            .to_col(&basis[i])
                    })
                    .collect(),
            );
            assert_eq!(mat.rank(), n);
        }
        Self { anf, basis }
    }

    pub fn anf(&self) -> &K {
        self.anf.borrow()
    }

    pub fn n(&self) -> usize {
        debug_assert_eq!(self.anf().n(), self.basis.len());
        self.basis.len()
    }

    pub fn z_module(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, IntegerCanonicalStructure> {
        Integer::structure().into_free_module(self.n())
    }

    pub fn discriminant(&self) -> Rational {
        self.anf()
            .finite_dimensional_rational_extension()
            .discriminant(&self.basis)
    }

    pub fn into_anf_extension(
        self,
    ) -> AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, Self> {
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion {
            _k: PhantomData,
            _kb: PhantomData,
            abelian_subgroup: self,
        }
    }

    pub fn anf_extension<'a>(
        &'a self,
    ) -> AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, &'a Self> {
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion {
            _k: PhantomData,
            _kb: PhantomData,
            abelian_subgroup: self,
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> PartialEq
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    fn eq(&self, other: &Self) -> bool {
        let n = self.anf().n();
        if self.anf == other.anf {
            debug_assert_eq!(other.anf().n(), n);
            debug_assert_eq!(self.basis.len(), n);
            debug_assert_eq!(other.basis.len(), n);
            (0..n).all(|i| {
                let b = self.anf().equal(&self.basis[i], &other.basis[i]);
                debug_assert_eq!(b, other.anf().equal(&self.basis[i], &other.basis[i]));
                b
            })
        } else {
            false
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Eq
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Signature
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> SetSignature
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if x.len() != self.n() {
            return Err("wrong length".to_string());
        }
        Ok(())
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> EqSignature
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.z_module().equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveMonoidSignature
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    fn zero(&self) -> Self::Set {
        self.z_module().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.z_module().add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveGroupSignature
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.z_module().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.z_module().sub(a, b)
    }
}

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    AB: BorrowedStructure<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>>,
> {
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    abelian_subgroup: AB,
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    AB: BorrowedStructure<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>>,
> Morphism<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>, K>
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, AB>
{
    fn domain(&self) -> &AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB> {
        self.abelian_subgroup.borrow()
    }

    fn range(&self) -> &K {
        self.abelian_subgroup.borrow().anf()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    AB: BorrowedStructure<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>>,
> Function<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>, K>
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, AB>
{
    fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
        debug_assert!(self.abelian_subgroup.borrow().is_element(x).is_ok());
        let k = self.abelian_subgroup.borrow().anf();
        let n = k.n();
        debug_assert_eq!(n, x.len());
        k.sum(
            (0..n)
                .map(|i| k.mul(&k.from_int(&x[i]), &self.abelian_subgroup.borrow().basis[i]))
                .collect(),
        )
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    AB: BorrowedStructure<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>>,
> InjectiveFunction<AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>, K>
    for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, AB>
{
    fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
        let k = self.abelian_subgroup.borrow().anf();
        let n = k.n();
        debug_assert!(k.is_element(y).is_ok());
        let mat = Matrix::join_cols(
            n,
            (0..n)
                .map(|i| {
                    k.finite_dimensional_rational_extension()
                        .to_col(&self.abelian_subgroup.borrow().basis[i])
                })
                .collect(),
        );
        let y = k.finite_dimensional_rational_extension().to_vec(y);
        let x_rat = mat.col_solve(&y)?;
        let mut x_int = Vec::with_capacity(n);
        for c_rat in x_rat {
            x_int.push(c_rat.try_to_int()?);
        }
        Some(x_int)
    }
}

#[cfg(test)]
mod tests {
    use crate::parsing::parse_rational_polynomial;

    use super::*;

    #[test]
    fn full_rank_abelian_group_to_algebraic_number_field_image_and_preimage() {
        let anf = parse_rational_polynomial("x^2+1", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();

        let abgrp = AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis::new(
            anf.clone(),
            vec![
                parse_rational_polynomial("2 + 3*x", "x").unwrap(),
                parse_rational_polynomial("1 - x", "x").unwrap(),
            ],
        )
        .unwrap();

        assert!(
            anf.equal(
                &abgrp
                    .anf_extension()
                    .image(&vec![Integer::from(1), Integer::from(2)]),
                &parse_rational_polynomial("4 + x", "x").unwrap()
            )
        );

        assert_eq!(
            abgrp
                .anf_extension()
                .try_preimage(&parse_rational_polynomial("4 + x", "x").unwrap()),
            Some(vec![Integer::from(1), Integer::from(2)])
        );

        assert!(
            abgrp
                .anf_extension()
                .try_preimage(&parse_rational_polynomial("x", "x").unwrap())
                .is_none()
        );
    }
}
