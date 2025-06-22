use super::{QuaternionAlgebraElement, QuaternionAlgebraStructure};
use crate::{
    algebraic_number_field::structure::AlgebraicNumberFieldSignature,
    module::finitely_free_module::RingToFinitelyFreeModuleSignature,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, FiniteRankFreeRingExtension,
        FinitelyFreeModuleSignature, SemiModuleSignature, SemiRingSignature,
    },
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Rational};
use algebraeon_sets::structure::{EqSignature, Function, MetaType, SetSignature, Signature};
use itertools::Itertools;

#[derive(Debug, Clone)]
pub struct QuaternionOrderZBasis<ANF: AlgebraicNumberFieldSignature> {
    integers: IntegerCanonicalStructure, // so we can return a reference to it in .ring()
    algebra: QuaternionAlgebraStructure<ANF>,
    basis: Vec<QuaternionAlgebraElement<ANF::Set>>, // 4n elements
}

impl<ANF: AlgebraicNumberFieldSignature> PartialEq for QuaternionOrderZBasis<ANF> {
    #[allow(clippy::overly_complex_bool_expr)]
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && false
    }
}

impl<ANF: AlgebraicNumberFieldSignature> Eq for QuaternionOrderZBasis<ANF> {}

impl<ANF: AlgebraicNumberFieldSignature> EqSignature for QuaternionOrderZBasis<ANF> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.algebra.equal(a, b)
    }
}

impl<ANF: AlgebraicNumberFieldSignature> Signature for QuaternionOrderZBasis<ANF> {}

impl<ANF: AlgebraicNumberFieldSignature> SetSignature for QuaternionOrderZBasis<ANF> {
    type Set = QuaternionAlgebraElement<ANF::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let algebra = &self.algebra;
        let submodules = Rational::structure()
            .into_free_module(self.basis.len())
            .into_submodules();

        let basis_vecs: Vec<Vec<Rational>> = self
            .basis
            .iter()
            .map(|b| {
                algebra
                    .to_vec(b)
                    .into_iter()
                    .flat_map(|p| {
                        algebra
                            .base_field()
                            .finite_dimensional_rational_extension()
                            .to_vec(&p)
                    })
                    .collect_vec()
            })
            .collect_vec();

        let basis_refs: Vec<&Vec<Rational>> = basis_vecs.iter().collect();

        let v = submodules.span(basis_refs);

        let x_vec = algebra
            .to_vec(x)
            .into_iter()
            .flat_map(|p| {
                algebra
                    .base_field()
                    .finite_dimensional_rational_extension()
                    .to_vec(&p)
            })
            .collect_vec();

        let (coset, element_reduced) = submodules.reduce_element(&v, &x_vec);
        debug_assert!(element_reduced.iter().all(|coeff| *coeff == Rational::ZERO));

        if coset.iter().all(|coeff| coeff.is_integer()) {
            Ok(())
        } else {
            Err("Element not in order".to_string())
        }
    }
}

impl<ANF: AlgebraicNumberFieldSignature> AdditiveMonoidSignature for QuaternionOrderZBasis<ANF> {
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

impl<ANF: AlgebraicNumberFieldSignature> AdditiveGroupSignature for QuaternionOrderZBasis<ANF> {
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

impl<ANF: AlgebraicNumberFieldSignature> SemiModuleSignature<IntegerCanonicalStructure>
    for QuaternionOrderZBasis<ANF>
{
    fn ring(&self) -> &IntegerCanonicalStructure {
        &self.integers
    }

    fn scalar_mul(&self, a: &Self::Set, x: &Integer) -> Self::Set {
        self.algebra.scalar_mul(
            a,
            &self
                .algebra
                .base_field()
                .principal_subring_inclusion()
                .image(x),
        )
    }
}

// impl ModuleSignature<IntegerCanonicalStructure> for QuaternionOrderZBasis {}

impl<ANF: AlgebraicNumberFieldSignature> QuaternionOrderZBasis<ANF> {
    fn check_basis(self) -> bool {
        // 1. check that 1 belongs to the order
        if self.is_element(&self.algebra.clone().one()).is_err() {
            return false;
        }

        // 2. Check that the basis is closed under multiplication
        for (i, bi) in self.basis.iter().enumerate() {
            for (j, bj) in self.basis.iter().enumerate() {
                let product = self.algebra.mul(bi, bj);
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

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;

    #[test]
    fn test_lipschitz_order() {
        let rat = Rational::structure();
        let a = -Rational::ONE;
        let b = -Rational::ONE;

        // Hamilton quaternion algebra: H = (-1, -1 / QQ)
        let h = QuaternionAlgebraStructure::new(rat.clone(), a.clone(), b.clone());
        let one = h.one();
        let i = h.i();
        let j = h.j();
        let k = h.k();

        // Lipschitz order <1, i, j, k>_QQ
        let order = QuaternionOrderZBasis {
            integers: Integer::structure(),
            algebra: h,
            basis: vec![one.clone(), i.clone(), j.clone(), k.clone()],
        };

        for (idx, b) in order.basis.iter().enumerate() {
            assert!(order.is_element(b).is_ok(),);
        }

        assert!(order.check_basis());
    }
}
