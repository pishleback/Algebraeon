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
                    self.try_anf_to_roi(&self.anf().mul(self.basis_element(i), &self.roi_to_anf(a)))
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

    fn ideal_intersect(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        Self::Ideal {
            lattice: LinearLatticeStructure::new(Integer::structure()).intersect_pair(
                self.degree(),
                1,
                &a.lattice,
                &b.lattice,
            ),
        }
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        Self::Ideal {
            lattice: LinearLatticeStructure::new(Integer::structure()).sum_pair(
                self.degree(),
                1,
                &a.lattice,
                &b.lattice,
            ),
        }
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        let n = self.degree();

        let a_basis = a
            .lattice
            .basis_matrices()
            .into_iter()
            .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
            .collect::<Vec<_>>();
        let b_basis = b
            .lattice
            .basis_matrices()
            .into_iter()
            .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
            .collect::<Vec<_>>();

        let mut span = vec![];
        for i in 0..n {
            for j in 0..n {
                span.push(self.mul(&a_basis[i], &b_basis[j]));
            }
        }
        Self::Ideal::from_integer_span(n, span)
    }
}

impl DedekindDomainStructure for RingOfIntegersWithIntegralBasisStructure {}

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

        {
            // 1 + sqrt(2)
            let alpha = roi.try_anf_to_roi(&(&x + 1).into_verbose()).unwrap();

            // (a + b sqrt(2)) * (1 + sqrt(2)) = a(1 + sqrt(2)) + b(2 + sqrt(2))
            assert!(roi.ideal_equal(
                &roi.principal_ideal(&alpha),
                &RingOfIntegersIdeal::from_integer_basis(
                    2,
                    vec![
                        roi.try_anf_to_roi(&(1 + &x).into_verbose()).unwrap(),
                        roi.try_anf_to_roi(&(2 + &x).into_verbose()).unwrap()
                    ]
                )
            ));
        }

        {
            // 6
            let alpha = roi.try_anf_to_roi(&(6 * x.pow(0)).into_verbose()).unwrap();
            // 15
            let beta = roi.try_anf_to_roi(&(15 * x.pow(0)).into_verbose()).unwrap();

            let alpha_ideal = roi.principal_ideal(&alpha);
            let beta_ideal = roi.principal_ideal(&beta);

            let alpha_beta_add = roi.ideal_add(&alpha_ideal, &beta_ideal);
            let alpha_beta_intersect = roi.ideal_intersect(&alpha_ideal, &beta_ideal);
            let alpha_beta_mul = roi.ideal_mul(&alpha_ideal, &beta_ideal);

            // sum is 3
            assert!(roi.ideal_equal(
                &alpha_beta_add,
                &RingOfIntegersIdeal::from_integer_basis(
                    2,
                    vec![
                        roi.try_anf_to_roi(&(3 * x.pow(0)).into_verbose()).unwrap(),
                        roi.try_anf_to_roi(&(3 * x.pow(1)).into_verbose()).unwrap()
                    ]
                )
            ));

            // intersection is 30
            assert!(roi.ideal_equal(
                &alpha_beta_intersect,
                &RingOfIntegersIdeal::from_integer_basis(
                    2,
                    vec![
                        roi.try_anf_to_roi(&(30 * x.pow(0)).into_verbose()).unwrap(),
                        roi.try_anf_to_roi(&(30 * x.pow(1)).into_verbose()).unwrap()
                    ]
                )
            ));

            // product is 90
            assert!(roi.ideal_equal(
                &alpha_beta_mul,
                &RingOfIntegersIdeal::from_integer_basis(
                    2,
                    vec![
                        roi.try_anf_to_roi(&(90 * x.pow(0)).into_verbose()).unwrap(),
                        roi.try_anf_to_roi(&(90 * x.pow(1)).into_verbose()).unwrap()
                    ]
                )
            ));
        }
    }
}
