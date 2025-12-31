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
    BorrowedStructure, EqSignature, Function, InjectiveFunction, Morphism, SetSignature, Signature,
    ToStringSignature,
};
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct FullRankZSubmoduleWithBasis<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> {
    anf: KB,
    // length = degree of k
    basis: Vec<K::Set>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    FullRankZSubmoduleWithBasis<K, KB>
{
    fn check(&self) -> Result<(), String> {
        let n = self.anf.borrow().n();
        if n != self.basis.len() {
            return Err("Basis has wrong length".to_string());
        }
        for v in &self.basis {
            if let Err(e) = self.anf.borrow().is_element(v) {
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
                    self.anf
                        .borrow()
                        .inbound_finite_dimensional_rational_extension()
                        .to_col(&self.basis[i])
                })
                .collect(),
        );
        if mat.rank() != n {
            return Err("Vectors do not form a basis".to_string());
        }
        Ok(())
    }

    fn new_impl(anf: KB, basis: Vec<K::Set>) -> Self {
        Self { anf, basis }
    }

    pub fn new(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let s = Self::new_impl(anf, basis);
        s.check()?;
        Ok(s)
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let s = Self::new_impl(anf, basis);
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }

    pub fn anf(&self) -> &K {
        self.anf.borrow()
    }

    pub fn basis(&self) -> &Vec<K::Set> {
        &self.basis
    }

    pub fn basis_vector(&self, i: usize) -> &K::Set {
        debug_assert!(i < self.n());
        &self.basis[i]
    }

    pub fn n(&self) -> usize {
        debug_assert_eq!(self.anf().n(), self.basis.len());
        self.basis.len()
    }

    pub fn free_z_module_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, &'static IntegerCanonicalStructure>
    {
        Integer::structure_ref().free_module(self.n())
    }

    pub fn discriminant(&self) -> Rational {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(&self.basis)
    }

    pub fn into_outbound_order_to_anf_inclusion(
        self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion::new(self)
    }

    pub fn outbound_order_to_anf_inclusion<'a>(
        &'a self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, &'a Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion::new(self)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> PartialEq
    for FullRankZSubmoduleWithBasis<K, KB>
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
    for FullRankZSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Signature
    for FullRankZSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> SetSignature
    for FullRankZSubmoduleWithBasis<K, KB>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if x.len() != self.n() {
            return Err("wrong length".to_string());
        }
        Ok(())
    }
}

impl<K: AlgebraicNumberFieldSignature + ToStringSignature, KB: BorrowedStructure<K>>
    ToStringSignature for FullRankZSubmoduleWithBasis<K, KB>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.anf()
            .to_string(&self.outbound_order_to_anf_inclusion().image(elem))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> EqSignature
    for FullRankZSubmoduleWithBasis<K, KB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.free_z_module_restructure().equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveMonoidSignature
    for FullRankZSubmoduleWithBasis<K, KB>
{
    fn zero(&self) -> Self::Set {
        self.free_z_module_restructure().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.free_z_module_restructure().add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveGroupSignature
    for FullRankZSubmoduleWithBasis<K, KB>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.free_z_module_restructure().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.free_z_module_restructure().sub(a, b)
    }
}

mod anf_inclusion {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        AB: BorrowedStructure<FullRankZSubmoduleWithBasis<K, KB>>,
    > {
        _k: PhantomData<K>,
        _kb: PhantomData<KB>,
        abelian_subgroup: AB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        AB: BorrowedStructure<FullRankZSubmoduleWithBasis<K, KB>>,
    > AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, AB>
    {
        pub fn new(abelian_subgroup: AB) -> Self {
            Self {
                _k: PhantomData,
                _kb: PhantomData,
                abelian_subgroup,
            }
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        AB: BorrowedStructure<FullRankZSubmoduleWithBasis<K, KB>>,
    > Morphism<FullRankZSubmoduleWithBasis<K, KB>, K>
        for AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion<K, KB, AB>
    {
        fn domain(&self) -> &FullRankZSubmoduleWithBasis<K, KB> {
            self.abelian_subgroup.borrow()
        }

        fn range(&self) -> &K {
            self.abelian_subgroup.borrow().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        AB: BorrowedStructure<FullRankZSubmoduleWithBasis<K, KB>>,
    > Function<FullRankZSubmoduleWithBasis<K, KB>, K>
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
        AB: BorrowedStructure<FullRankZSubmoduleWithBasis<K, KB>>,
    > InjectiveFunction<FullRankZSubmoduleWithBasis<K, KB>, K>
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
                        k.inbound_finite_dimensional_rational_extension()
                            .to_col(&self.abelian_subgroup.borrow().basis[i])
                    })
                    .collect(),
            );
            let y = k.inbound_finite_dimensional_rational_extension().to_vec(y);
            let x_rat = mat.col_solve(&y)?;
            let mut x_int = Vec::with_capacity(n);
            for c_rat in x_rat {
                x_int.push(c_rat.try_to_int()?);
            }
            Some(x_int)
        }
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

        let abgrp = FullRankZSubmoduleWithBasis::new(
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
                    .outbound_order_to_anf_inclusion()
                    .image(&vec![Integer::from(1), Integer::from(2)]),
                &parse_rational_polynomial("4 + x", "x").unwrap()
            )
        );

        assert_eq!(
            abgrp
                .outbound_order_to_anf_inclusion()
                .try_preimage(&parse_rational_polynomial("4 + x", "x").unwrap()),
            Some(vec![Integer::from(1), Integer::from(2)])
        );

        assert!(
            abgrp
                .outbound_order_to_anf_inclusion()
                .try_preimage(&parse_rational_polynomial("x", "x").unwrap())
                .is_none()
        );
    }
}
