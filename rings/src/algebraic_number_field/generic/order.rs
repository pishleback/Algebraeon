use std::marker::PhantomData;

use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis,
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasisInclusion,
        AlgebraicNumberFieldSignature,
    },
    matrix::SymmetricMatrix,
    module::finitely_free_module::FinitelyFreeModuleStructure,
    structure::{FiniteDimensionalFieldExtension, MetaCharZeroRing},
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, Function, InjectiveFunction, Morphism, SetSignature, Signature,
};

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldOrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
> {
    abelian_group: AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>,
    products: SymmetricMatrix<K::Set>,
}

fn make_products<K: AlgebraicNumberFieldSignature>(
    anf: &K,
    basis: &Vec<K::Set>,
) -> SymmetricMatrix<K::Set> {
    SymmetricMatrix::construct_bottom_left(anf.n(), |r, c| anf.mul(&basis[r], &basis[c]))
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    AlgebraicNumberFieldOrderWithBasis<K, KB>
{
    fn check(&self) -> Result<(), String> {
        for b in self.abelian_group.basis() {
            if !self.abelian_group.anf().is_algebraic_integer(b) {
                return Err("Basis vectors must be algebraic integers for an order".to_string());
            }
        }
        compile_error!("Also need to check closure under products");
        Ok(())
    }

    pub fn new(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let s = Self::new_unchecked(anf, basis);
        s.check()?;
        Ok(s)
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let products = make_products(anf.borrow(), &basis);
        let abelian_group =
            AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis::new_unchecked(anf, basis);
        let s = Self {
            abelian_group,
            products,
        };
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }

    pub fn anf(&self) -> &K {
        self.abelian_group.anf()
    }

    pub fn basis(&self) -> &Vec<K::Set> {
        self.abelian_group.basis()
    }

    pub fn n(&self) -> usize {
        self.abelian_group.n()
    }

    pub fn z_module(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, IntegerCanonicalStructure> {
        self.abelian_group.z_module()
    }

    pub fn discriminant(&self) -> Integer {
        self.anf()
            .finite_dimensional_rational_extension()
            .discriminant(self.abelian_group.basis())
            .try_to_int()
            .unwrap()
    }

    pub fn into_anf_extension(self) -> AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, Self> {
        AlgebraicNumberFieldOrderWithBasisInclusion {
            _k: PhantomData,
            _kb: PhantomData,
            order: self,
        }
    }

    pub fn anf_extension<'a>(
        &'a self,
    ) -> AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, &'a Self> {
        AlgebraicNumberFieldOrderWithBasisInclusion {
            _k: PhantomData,
            _kb: PhantomData,
            order: self,
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> PartialEq
    for AlgebraicNumberFieldOrderWithBasis<K, KB>
{
    fn eq(&self, other: &Self) -> bool {
        self.abelian_group == other.abelian_group
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Eq
    for AlgebraicNumberFieldOrderWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Signature
    for AlgebraicNumberFieldOrderWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> SetSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.abelian_group.is_element(x)
    }
}

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldOrderWithBasisInclusion<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> {
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    order: OB,
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> Morphism<AlgebraicNumberFieldOrderWithBasis<K, KB>, K>
    for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, OB>
{
    fn domain(&self) -> &AlgebraicNumberFieldOrderWithBasis<K, KB> {
        self.order.borrow()
    }

    fn range(&self) -> &K {
        self.order.borrow().anf()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> Function<AlgebraicNumberFieldOrderWithBasis<K, KB>, K>
    for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, OB>
{
    fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
        todo!()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> InjectiveFunction<AlgebraicNumberFieldOrderWithBasis<K, KB>, K>
    for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, OB>
{
    fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
        todo!()
    }
}
