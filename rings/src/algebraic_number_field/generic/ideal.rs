use crate::algebraic_number_field::{
    AlgebraicNumberFieldOrderWithBasis, AlgebraicNumberFieldSignature,
};
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::BorrowedStructure;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldIdeal<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> {
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    order: OB,
    ideal: Vec<Integer>,
}

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldFractionalIdeal<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB>>,
> {
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    order: OB,
    ideal: Vec<Integer>,
    denominator: Integer,
}
