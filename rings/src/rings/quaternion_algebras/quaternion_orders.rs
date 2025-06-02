use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Rational, RationalCanonicalStructure};
use algebraeon_sets::structure::{EqSignature, Function, MetaType, SetSignature, Signature};

use crate::{
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    rings::algebraic_number_fields::number_field::AlgebraicNumberFieldStructure,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, FinitelyFreeModuleSignature,
        FreeModuleSignature, ModuleSignature, RingSignature, SemiModuleSignature,
        SemiRingSignature,
    },
};

use super::{QuaternionAlgebraElement, QuaternionAlgebraStructure};

#[derive(Debug, Clone)]
pub struct QuaternionOrderZBasis {
    integers: IntegerCanonicalStructure, // so we can return a reference to it in .ring()
    algebra: QuaternionAlgebraStructure<AlgebraicNumberFieldStructure>,
    basis: Vec<QuaternionAlgebraElement<AlgebraicNumberFieldStructure>>, // 4n elements
}

impl PartialEq for QuaternionOrderZBasis {
    #[allow(clippy::overly_complex_bool_expr)]
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

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let submodules = self
            .algebra
            .ring()
            .free_module(self.algebra.rank())
            .submodules();

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

// impl FreeModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {
//     type Basis;

//     fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
//         todo!()
//     }

//     fn to_component<'a>(
//         &'a self,
//         b: &<Self::Basis as SetSignature>::Set,
//         v: &'a Self::Set,
//     ) -> &'a Integer {
//         todo!()
//     }

//     fn from_component(&self, b: &<Self::Basis as SetSignature>::Set, r: &Integer) -> Self::Set {
//         todo!()
//     }
// }

impl SemiModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {
    fn ring(&self) -> &IntegerCanonicalStructure {
        &self.integers
    }

    fn scalar_mul(&self, x: &Integer, a: &Self::Set) -> Self::Set {
        self.algebra.scalar_mul(
            &self
                .algebra
                .base_field()
                .principal_subring_inclusion()
                .image(x),
            a,
        )
    }
}

// impl ModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {}

impl QuaternionOrderZBasis {
    fn check_basis(self) -> bool {
        // 1. check that 1 belongs to the order
        if self.is_element(&self.algebra.clone().one()).is_err() {
            return false;
        }

        // 2. Check that the basis is closed under multiplication
        for (i, bi) in self.basis.iter().enumerate() {
            for (j, bj) in self.basis.iter().enumerate() {
                let product = self.algebra.mul(&bi, &bj);
                if self.is_element(&product).is_err() {
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
