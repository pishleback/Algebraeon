use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis, AlgebraicNumberFieldSignature,
    },
    matrix::SymmetricMatrix,
};
use algebraeon_sets::structure::{BorrowedStructure, Signature};

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldOrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
> {
    abelian_group: AlgebraicNumberFieldFullRankAbelianSubgroupWithBasis<K, KB>,
    products: SymmetricMatrix<K::Set>,
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
