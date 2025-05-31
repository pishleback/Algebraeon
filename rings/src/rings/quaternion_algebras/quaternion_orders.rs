use algebraeon_nzq::{Integer, IntegerCanonicalStructure};
use algebraeon_sets::structure::{EqSignature, MetaType, SetSignature, Signature};

use crate::{
    rings::algebraic_number_fields::number_field::AlgebraicNumberFieldStructure,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, FreeModuleSignature, ModuleSignature,
        SemiRingSignature,
    },
};

use super::{QuaternionAlgebraElement, QuaternionAlgebraStructure};

#[derive(Debug, Clone)]
pub struct QuaternionOrderZBasis {
    algebra: QuaternionAlgebraStructure<AlgebraicNumberFieldStructure>,
    basis: Vec<QuaternionAlgebraElement<AlgebraicNumberFieldStructure>>, // 4n elements
}

impl PartialEq for QuaternionOrderZBasis {
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && false
    }
}

impl Eq for QuaternionOrderZBasis {}

impl EqSignature for QuaternionOrderZBasis {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.algebra.equal(&a, &b)
    }
}

impl Signature for QuaternionOrderZBasis {}

impl SetSignature for QuaternionOrderZBasis {
    type Set = QuaternionAlgebraElement<AlgebraicNumberFieldStructure>;

    fn is_element(&self, x: &Self::Set) -> bool {
        // let submodules = self.algebra.submodules();
        unimplemented!("linear algebra")
    }
}

impl AdditiveMonoidSignature for QuaternionOrderZBasis {
    fn is_zero(&self, a: &Self::Set) -> bool {
        self.algebra.is_zero(a)
    }

    fn zero(&self) -> Self::Set {
        self.algebra.zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.algebra.add(a, b)
    }
}

impl AdditiveGroupSignature for QuaternionOrderZBasis {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.algebra.neg(a)
    }
}

impl FreeModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {}

impl ModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {
    fn ring(&self) -> &IntegerCanonicalStructure {
        &Integer::structure().clone()
    }

    fn scalar_mul(&self, x: &<IntegerCanonicalStructure>::Set, a: &Self::Set) -> Self::Set {
        self.algebra.scalar_mul(x, a)
    }
}

impl QuaternionOrderZBasis {
    fn check_basis(self) -> bool {
        // 1. check that 1 belongs to the order
        if !self.is_element(&self.algebra.clone().one()) {
            return false;
        }

        // 2. Check that the basis is closed under multiplication
        for (i, bi) in self.basis.iter().enumerate() {
            for (j, bj) in self.basis.iter().enumerate() {
                let product = self.algebra.mul(&bi, &bj);
                if !&self.is_element(&product) {
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
