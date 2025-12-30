use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis, AlgebraicNumberFieldSignature,
    },
    matrix::SymmetricMatrix,
    module::finitely_free_module::FinitelyFreeModuleStructure,
    structure::{FiniteDimensionalFieldExtension, FiniteRankFreeRingExtension, MetaCharZeroRing},
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, Function, InjectiveFunction, Morphism, SetSignature, Signature,
};
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldOrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> {
    full_rank_abelian_group: AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>,
    products: SymmetricMatrix<Vec<Integer>>,
}

fn make_products<K: AlgebraicNumberFieldSignature>(
    anf: &K,
    basis: &Vec<K::Set>,
) -> Result<SymmetricMatrix<Vec<Integer>>, String> {
    let n = anf.n();
    let mut s = SymmetricMatrix::filled(n, vec![]);
    for r in 0..n {
        for c in r..n {
            let mut x = vec![];
            for y in anf
                .inbound_finite_dimensional_rational_extension()
                .to_vec(&anf.mul(&basis[r], &basis[c]))
            {
                if let Some(y) = y.try_to_int() {
                    x.push(y)
                } else {
                    return Err(format!("{} is not an integer", y));
                }
            }
            s.set(r, c, x).unwrap();
        }
    }
    Ok(s)
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn check(&self) -> Result<(), String> {
        for b in self.full_rank_abelian_group.basis() {
            if !self.full_rank_abelian_group.anf().is_algebraic_integer(b) {
                return Err("Basis vectors must be algebraic integers for an order".to_string());
            }
        }
        Ok(())
    }

    fn new_impl(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let products = make_products(anf.borrow(), &basis)
            .map_err(|err| format!("Not an order: Not closed under multiplication: {}", err))
            .unwrap();
        let abelian_group =
            AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis::new_unchecked(anf, basis);
        let s = Self {
            full_rank_abelian_group: abelian_group,
            products,
        };
        Ok(s)
    }

    pub fn new(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let s = Self::new_impl(anf, basis)?;
        s.check()?;
        Ok(s)
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let s = Self::new_impl(anf, basis).unwrap();
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }

    pub fn anf(&self) -> &K {
        self.full_rank_abelian_group.anf()
    }

    pub fn basis(&self) -> &Vec<K::Set> {
        self.full_rank_abelian_group.basis()
    }

    pub fn n(&self) -> usize {
        self.full_rank_abelian_group.n()
    }

    pub fn full_rank_abelian_group_restructure(
        &self,
    ) -> &AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB> {
        &self.full_rank_abelian_group
    }

  pub  fn abelian_group_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, IntegerCanonicalStructure> {
        self.full_rank_abelian_group.abelian_group_restructure()
    }

    pub fn discriminant(&self) -> Integer {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(self.full_rank_abelian_group.basis())
            .try_to_int()
            .unwrap()
    }

    pub fn into_outbound_anf_inclusion(
        self,
    ) -> anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, Self> {
        anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion::new(self)
    }

    pub fn outbound_anf_inclusion<'a>(
        &'a self,
    ) -> anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, &'a Self> {
        anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion::new(self)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> PartialEq
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn eq(&self, other: &Self) -> bool {
        self.full_rank_abelian_group == other.full_rank_abelian_group
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Eq
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Signature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> SetSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.full_rank_abelian_group.is_element(x)
    }
}

mod anf_inclusion {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldOrderWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > {
        _k: PhantomData<K>,
        _kb: PhantomData<KB>,
        order: OB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        pub fn new(order: OB) -> Self {
            Self {
                _k: PhantomData,
                _kb: PhantomData,
                order,
            }
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > Morphism<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn domain(&self) -> &AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL> {
            self.order.borrow()
        }

        fn range(&self) -> &K {
            self.order.borrow().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > Function<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
            todo!()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > InjectiveFunction<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
            todo!()
        }
    }
}
