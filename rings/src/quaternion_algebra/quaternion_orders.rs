use super::{QuaternionAlgebraElement, QuaternionAlgebraStructure};
use crate::{
    algebraic_number_field::AlgebraicNumberFieldSignature,
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, FiniteRankFreeRingExtension, FinitelyFreeModuleSignature,
        MultiplicationSignature, OneSignature, RinglikeSpecializationSignature,
        SemiModuleSignature, TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use itertools::Itertools;

// an integer submodule and subring of a quaternion algebra over an algebraic number field
#[derive(Debug, Clone)]
pub struct QuaternionOrderZBasis<ANF: AlgebraicNumberFieldSignature> {
    integers: IntegerCanonicalStructure, // so we can return a reference to it in .ring()
    algebra: QuaternionAlgebraStructure<ANF>,
    basis: Vec<QuaternionAlgebraElement<ANF::Elem>>, // 4n elements
}

impl<ANF: AlgebraicNumberFieldSignature> PartialEq for QuaternionOrderZBasis<ANF> {
    #[allow(clippy::overly_complex_bool_expr)]
    fn eq(&self, other: &Self) -> bool {
        self.algebra == other.algebra && false
    }
}

impl<ANF: AlgebraicNumberFieldSignature> Eq for QuaternionOrderZBasis<ANF> {}

impl<ANF: AlgebraicNumberFieldSignature> EqSignature for QuaternionOrderZBasis<ANF> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.algebra.equal(a, b)
    }
}

impl<ANF: AlgebraicNumberFieldSignature> Signature for QuaternionOrderZBasis<ANF> {}

impl<ANF: AlgebraicNumberFieldSignature> SetSignature for QuaternionOrderZBasis<ANF> {
    type Elem = QuaternionAlgebraElement<ANF::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        let algebra = &self.algebra;
        let submodules = Rational::structure()
            .into_free_module(EnumeratedFiniteSetStructure::new(self.basis.len()))
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
                            .inbound_finite_dimensional_rational_extension()
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
                    .inbound_finite_dimensional_rational_extension()
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

impl<ANF: AlgebraicNumberFieldSignature> RinglikeSpecializationSignature
    for QuaternionOrderZBasis<ANF>
{
}

impl<ANF: AlgebraicNumberFieldSignature> ZeroSignature for QuaternionOrderZBasis<ANF> {
    fn zero(&self) -> Self::Elem {
        self.algebra.zero()
    }
}

impl<ANF: AlgebraicNumberFieldSignature> AdditionSignature for QuaternionOrderZBasis<ANF> {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.algebra.add(a, b)
    }
}

impl<ANF: AlgebraicNumberFieldSignature> CancellativeAdditionSignature
    for QuaternionOrderZBasis<ANF>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<ANF: AlgebraicNumberFieldSignature> TryNegateSignature for QuaternionOrderZBasis<ANF> {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<ANF: AlgebraicNumberFieldSignature> AdditiveMonoidSignature for QuaternionOrderZBasis<ANF> {}

impl<ANF: AlgebraicNumberFieldSignature> AdditiveGroupSignature for QuaternionOrderZBasis<ANF> {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.algebra.neg(a)
    }
}

impl<ANF: AlgebraicNumberFieldSignature> SemiModuleSignature<IntegerCanonicalStructure>
    for QuaternionOrderZBasis<ANF>
{
    fn ring(&self) -> &IntegerCanonicalStructure {
        &self.integers
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &Integer) -> Self::Elem {
        self.algebra.scalar_mul(
            a,
            &self
                .algebra
                .base_field()
                .inbound_principal_integer_map()
                .image(x),
        )
    }
}

impl<ANF: AlgebraicNumberFieldSignature> QuaternionOrderZBasis<ANF> {
    #[allow(unused)]
    fn check_basis(self) -> bool {
        // 1. Check that the basis has 4n elements
        let expected_len = 4 * self.algebra.base_field().n();

        if self.basis.len() != expected_len {
            println!(
                "Incorrect basis size: got {}, expected {}",
                self.basis.len(),
                expected_len
            );
            return false;
        }

        // 2. check that 1 belongs to the order
        if self.validate_element(&self.algebra.clone().one()).is_err() {
            return false;
        }

        // 3. Check that the basis is closed under multiplication
        for (i, bi) in self.basis.iter().enumerate() {
            for (j, bj) in self.basis.iter().enumerate() {
                let product = self.algebra.mul(bi, bj);
                if self.validate_element(&product).is_err() {
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

        for b in order.basis.iter() {
            assert!(order.validate_element(b).is_ok(),);
        }

        assert!(order.check_basis());
    }
}
