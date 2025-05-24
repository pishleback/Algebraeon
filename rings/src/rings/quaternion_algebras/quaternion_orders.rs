use crate::{
    rings::algebraic_number_fields::number_field::AlgebraicNumberFieldStructure,
    structure::SemiRingSignature,
};

use super::{QuaternionAlgebraElement, QuaternionAlgebraStructure};

#[derive(Debug)]
pub struct QuaternionOrderZBasis {
    algebra: QuaternionAlgebraStructure<AlgebraicNumberFieldStructure>,
    basis: Vec<QuaternionAlgebraElement<AlgebraicNumberFieldStructure>>, // 4n elements
}

impl QuaternionOrderZBasis {
    fn contains(&self, q: &QuaternionAlgebraElement<AlgebraicNumberFieldStructure>) -> bool {
        unimplemented!("linear algebra")
    }

    fn check_basis(self) -> bool {
        // 1. check that 1 belongs to the order
        if !&self.contains(&self.algebra.clone().one()) {
            return false;
        }

        // 2. Check that the basis is closed under multiplication
        for (i, bi) in self.basis.iter().enumerate() {
            for (j, bj) in self.basis.iter().enumerate() {
                let product = self.algebra.mul(&bi, &bj);
                if !&self.contains(&product) {
                    println!(
                        "Basis not closed under multiplication: b[{}] * b[{}] = {:?} not in order",
                        i, j, product
                    );
                    return false;
                }
            }
        }

        true
    }
}
