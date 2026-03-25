use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldSignature, FullRankIntegerSubmoduleWithBasisSignature,
    },
    matrix::Matrix,
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, FiniteRankFreeRingExtension,
        RinglikeSpecializationSignature, TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, Function, SetSignature, Signature, ToStringSignature,
};

#[derive(Debug, Clone)]
pub struct FullRankIntegerSubmoduleWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
> {
    anf: KB,
    // length = degree of k
    basis: Vec<K::Elem>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn check(&self) -> Result<(), String> {
        let n = self.anf.borrow().n();
        if n != self.basis.len() {
            return Err("Basis has wrong length".to_string());
        }
        for v in &self.basis {
            if let Err(e) = self.anf.borrow().validate_element(v) {
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

    fn new_impl(anf: KB, basis: Vec<K::Elem>) -> Self {
        Self { anf, basis }
    }

    pub fn new(anf: KB, basis: Vec<K::Elem>) -> Result<Self, String> {
        let s = Self::new_impl(anf, basis);
        s.check()?;
        Ok(s)
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Elem>) -> Self {
        let s = Self::new_impl(anf, basis);
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> PartialEq
    for FullRankIntegerSubmoduleWithBasis<K, KB>
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
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> Signature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> SetSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    type Elem = Vec<Integer>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if x.len() != self.n() {
            return Err("wrong length".to_string());
        }
        Ok(())
    }
}

impl<K: AlgebraicNumberFieldSignature + ToStringSignature, KB: BorrowedStructure<K>>
    ToStringSignature for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn to_string(&self, elem: &Self::Elem) -> String {
        self.anf()
            .to_string(&self.outbound_order_to_anf_inclusion().image(elem))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> EqSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.free_integer_submodule_restructure().equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> RinglikeSpecializationSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> ZeroSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn zero(&self) -> Self::Elem {
        self.free_integer_submodule_restructure().zero()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditionSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.free_integer_submodule_restructure().add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> CancellativeAdditionSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> TryNegateSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveMonoidSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AdditiveGroupSignature
    for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.free_integer_submodule_restructure().neg(a)
    }

    fn sub(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.free_integer_submodule_restructure().sub(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    FullRankIntegerSubmoduleWithBasisSignature<K> for FullRankIntegerSubmoduleWithBasis<K, KB>
{
    fn anf(&self) -> &K {
        self.anf.borrow()
    }

    fn basis(&self) -> &Vec<<K>::Elem> {
        &self.basis
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parsing::parse_rational_polynomial;
    use algebraeon_sets::structure::InjectiveFunction;

    #[test]
    fn full_rank_abelian_group_to_algebraic_number_field_image_and_preimage() {
        let anf = parse_rational_polynomial("x^2+1", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();

        let abgrp = FullRankIntegerSubmoduleWithBasis::new(
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
