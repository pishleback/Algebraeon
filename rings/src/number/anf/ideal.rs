use super::ring_of_integers::*;
use crate::{linear::subspace::*, structure::*};
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::MetaType;

#[derive(Debug, Clone)]
pub struct RingOfIntegersIdeal {
    // 1 column and n rows
    lattice: LinearLattice<Integer>,
}

impl RingOfIntegersWithIntegralBasisStructure {
    #[cfg(debug_assertions)]
    fn check_ideal(&self, ideal: &RingOfIntegersIdeal) {
        assert_eq!(ideal.lattice.rows(), self.degree());
        assert_eq!(ideal.lattice.cols(), 1);
        //TODO: check that it's actually an ideal
    }
}

impl RingOfIntegersIdeal {
    /// A basis of this ideal as a Z-module.
    pub fn integer_basis(&self) -> Vec<RingOfIntegersWithIntegralBasisElement> {
        self.lattice
            .basis_matrices()
            .into_iter()
            .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
            .collect()
    }

    /// Construct an ideal from a Z-linear span
    fn from_integer_span(n: usize, span: Vec<RingOfIntegersWithIntegralBasisElement>) -> Self {
        Self {
            lattice: LinearLattice::from_span(
                n,
                1,
                span.into_iter().map(|elem| elem.into_col()).collect(),
            ),
        }
    }

    /// Construct an ideal from a Z-linear basis
    fn from_integer_basis(n: usize, basis: Vec<RingOfIntegersWithIntegralBasisElement>) -> Self {
        Self {
            lattice: LinearLattice::from_basis(
                n,
                1,
                basis.into_iter().map(|elem| elem.into_col()).collect(),
            ),
        }
    }
}

impl IdealStructure for RingOfIntegersWithIntegralBasisStructure {
    type Ideal = RingOfIntegersIdeal;
}

impl IdealArithmeticStructure for RingOfIntegersWithIntegralBasisStructure {
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal {
        let n = self.degree();
        Self::Ideal::from_integer_basis(
            n,
            (0..n)
                .map(|i| {
                    self.anf_to_roi(self.anf().mul(self.basis_element(i), &self.roi_to_anf(a)))
                        .unwrap()
                })
                .collect(),
        )
    }

    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        LinearLatticeStructure::new(Integer::structure()).equal(&a.lattice, &b.lattice)
    }

    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        LinearLatticeStructure::new(Integer::structure())
            .contains_sublattice(&a.lattice, &b.lattice)
    }

    fn ideal_intersection(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::*;
    use algebraeon_nzq::*;

    #[test]
    fn ring_of_integers_ideals() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        let a = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE, Rational::ZERO]);
        let b = Polynomial::<Rational>::from_coeffs(vec![Rational::ZERO, Rational::ONE]);

        // Q[sqrt(2)]
        let anf = (x.pow(2) - 2).into_verbose().algebraic_number_field();
        let roi = RingOfIntegersWithIntegralBasisStructure::new(
            anf.clone(),
            vec![a.clone(), b.clone()],
            Integer::from(8),
        );

        // 1 + sqrt(2)
        let alpha = roi.anf_to_roi((&x + 1).into_verbose()).unwrap();

        // (a + b sqrt(2)) * (1 + sqrt(2)) = a(1 + sqrt(2)) + b(2 + sqrt(2))
        assert!(roi.ideal_equal(
            &roi.principal_ideal(&alpha),
            &RingOfIntegersIdeal::from_integer_basis(
                2,
                vec![
                    roi.anf_to_roi((1 + &x).into_verbose()).unwrap(),
                    roi.anf_to_roi((2 + &x).into_verbose()).unwrap()
                ]
            )
        ));
    }
}
