use super::{ring_of_integer_extensions::RingOfIntegersExtension, ring_of_integers::*};
use crate::{linear::subspace::*, structure::*};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::{
    combinatorics::num_partitions_part_pool,
    structure::{MetaType, SetStructure},
};
use itertools::Itertools;

#[derive(Debug, Clone)]
pub enum RingOfIntegersIdeal {
    Zero,
    NonZero {
        // 1 column and n rows
        lattice: LinearLattice<Integer>,
    },
}

impl RingOfIntegersWithIntegralBasisStructure {
    #[cfg(debug_assertions)]
    pub fn check_ideal(&self, ideal: &RingOfIntegersIdeal) {
        match ideal {
            RingOfIntegersIdeal::Zero => {}
            RingOfIntegersIdeal::NonZero { lattice } => {
                assert_eq!(lattice.rows(), self.degree());
                assert_eq!(lattice.cols(), 1);
                // check it's actually an ideal
                for ideal_basis_elem in lattice
                    .basis_matrices()
                    .into_iter()
                    .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
                {
                    for integral_basis_elem in (0..self.degree()).map(|i| {
                        RingOfIntegersWithIntegralBasisElement::basis_element(self.degree(), i)
                    }) {
                        let x = self.mul(&ideal_basis_elem, &integral_basis_elem).into_col();
                        assert!(lattice.contains_point(&x));
                    }
                }
            }
        }
    }
}

impl RingOfIntegersIdeal {
    /// A basis of this ideal as a Z-module.
    pub fn integer_basis(&self) -> Option<Vec<RingOfIntegersWithIntegralBasisElement>> {
        match self {
            RingOfIntegersIdeal::Zero => None,
            RingOfIntegersIdeal::NonZero { lattice } => Some(
                lattice
                    .basis_matrices()
                    .into_iter()
                    .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
                    .collect(),
            ),
        }
    }
}

impl RingOfIntegersWithIntegralBasisStructure {
    /// Construct an ideal from a Z-linear span
    pub fn ideal_from_integer_span(
        &self,
        span: Vec<RingOfIntegersWithIntegralBasisElement>,
    ) -> RingOfIntegersIdeal {
        for elem in &span {
            debug_assert!(self.is_element(elem));
        }
        let n = self.degree();
        RingOfIntegersIdeal::NonZero {
            lattice: LinearLattice::from_span(
                n,
                1,
                span.into_iter().map(|elem| elem.into_col()).collect(),
            ),
        }
    }

    /// Construct an ideal from a Z-linear basis
    fn ideal_from_integer_basis(
        &self,
        basis: Vec<RingOfIntegersWithIntegralBasisElement>,
    ) -> RingOfIntegersIdeal {
        for elem in &basis {
            debug_assert!(self.is_element(elem));
        }
        let n = self.degree();
        RingOfIntegersIdeal::NonZero {
            lattice: LinearLattice::from_basis(
                n,
                1,
                basis.into_iter().map(|elem| elem.into_col()).collect(),
            ),
        }
    }

    pub fn ideal_norm(&self, ideal: &RingOfIntegersIdeal) -> Natural {
        RingOfIntegersExtension::new_integer_extension(self.clone()).ideal_norm(ideal)
    }
}

impl RingOfIntegersWithIntegralBasisStructure {
    /// generate all ideals of norm equal to n
    pub fn all_ideals_norm_eq<'a>(
        &'a self,
        n: &Natural,
    ) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        match Integer::factor_ideal(n) {
            Some(n) => {
                let sq = RingOfIntegersExtension::new_integer_extension(self.clone());
                Box::new(
                    n.into_factor_powers()
                        .into_iter()
                        .map(|(p, k)| {
                            let k: usize = k.try_into().unwrap();
                            let primes_over_p = sq.factor_prime_ideal(p).into_factors();
                            num_partitions_part_pool(
                                k,
                                primes_over_p
                                    .iter()
                                    .map(|f| f.residue_class_degree)
                                    .collect(),
                            )
                            .into_iter()
                            .map(|idxs| {
                                self.ideal_product(
                                    idxs.into_iter()
                                        .map(|i| primes_over_p[i].prime_ideal.ideal().clone())
                                        .collect(),
                                )
                            })
                            .collect::<Vec<RingOfIntegersIdeal>>()
                        })
                        .multi_cartesian_product()
                        .map(|ideals| self.ideal_product(ideals)),
                )
            }
            None => Box::new(vec![self.zero_ideal()].into_iter()),
        }
    }

    /// generate all non-zero ideals of norm at most n
    pub fn all_nonzero_ideals_norm_le<'a>(
        &'a self,
        n: &'a Natural,
    ) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        Box::new(
            (1usize..)
                .map(|m| Natural::from(m))
                .take_while(|m| m <= n)
                .map(|m| self.all_ideals_norm_eq(&m))
                .flatten(),
        )
    }

    /// generate all ideals
    pub fn all_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        Box::new(
            (0usize..)
                .map(|m| Natural::from(m))
                .map(|m| self.all_ideals_norm_eq(&m))
                .flatten(),
        )
    }

    /// generate all non-zero ideals
    pub fn all_nonzero_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        Box::new(
            (1usize..)
                .map(|m| Natural::from(m))
                .map(|m| self.all_ideals_norm_eq(&m))
                .flatten(),
        )
    }
}

impl IdealStructure for RingOfIntegersWithIntegralBasisStructure {
    type Ideal = RingOfIntegersIdeal;
}

impl IdealArithmeticStructure for RingOfIntegersWithIntegralBasisStructure {
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal {
        if self.is_zero(a) {
            Self::Ideal::Zero
        } else {
            let n = self.degree();
            let ideal = self.ideal_from_integer_basis(
                (0..n)
                    .map(|i| {
                        self.try_anf_to_roi(
                            &self.anf().mul(self.basis_element(i), &self.roi_to_anf(a)),
                        )
                        .unwrap()
                    })
                    .collect(),
            );
            #[cfg(debug_assertions)]
            self.check_ideal(&ideal);
            ideal
        }
    }

    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        match (a, b) {
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::Zero) => true,
            (
                RingOfIntegersIdeal::NonZero { lattice: a_lattice },
                RingOfIntegersIdeal::NonZero { lattice: b_lattice },
            ) => LinearLatticeStructure::new(Integer::structure()).equal(a_lattice, b_lattice),
            _ => false,
        }
    }

    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }

        match (a, b) {
            (_, RingOfIntegersIdeal::Zero) => true,
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::NonZero { .. }) => {
                debug_assert_ne!(self.degree(), 0);
                false
            }
            (
                RingOfIntegersIdeal::NonZero { lattice: a_lattice },
                RingOfIntegersIdeal::NonZero { lattice: b_lattice },
            ) => LinearLatticeStructure::new(Integer::structure())
                .contains_sublattice(a_lattice, b_lattice),
        }
    }

    fn ideal_intersect(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        match (a, b) {
            (
                RingOfIntegersIdeal::NonZero { lattice: a_lattice },
                RingOfIntegersIdeal::NonZero { lattice: b_lattice },
            ) => Self::Ideal::NonZero {
                lattice: LinearLatticeStructure::new(Integer::structure()).intersect_pair(
                    self.degree(),
                    1,
                    a_lattice,
                    b_lattice,
                ),
            },
            _ => Self::Ideal::Zero,
        }
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        match (a, b) {
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::Zero) => RingOfIntegersIdeal::Zero,
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::NonZero { .. }) => b.clone(),
            (RingOfIntegersIdeal::NonZero { .. }, RingOfIntegersIdeal::Zero) => a.clone(),
            (
                RingOfIntegersIdeal::NonZero { lattice: a_lattice },
                RingOfIntegersIdeal::NonZero { lattice: b_lattice },
            ) => Self::Ideal::NonZero {
                lattice: LinearLatticeStructure::new(Integer::structure()).sum_pair(
                    self.degree(),
                    1,
                    a_lattice,
                    b_lattice,
                ),
            },
        }
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }

        match (a, b) {
            (
                RingOfIntegersIdeal::NonZero { lattice: a_lattice },
                RingOfIntegersIdeal::NonZero { lattice: b_lattice },
            ) => {
                let n = self.degree();

                let a_basis = a_lattice
                    .basis_matrices()
                    .into_iter()
                    .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
                    .collect::<Vec<_>>();
                let b_basis = b_lattice
                    .basis_matrices()
                    .into_iter()
                    .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
                    .collect::<Vec<_>>();

                debug_assert_eq!(a_basis.len(), n);
                debug_assert_eq!(b_basis.len(), n);

                let mut span = vec![];
                for i in 0..n {
                    for j in 0..n {
                        span.push(self.mul(&a_basis[i], &b_basis[j]));
                    }
                }
                self.ideal_from_integer_span(span)
            }
            _ => Self::Ideal::Zero,
        }
    }
}

impl DedekindDomainStructure for RingOfIntegersWithIntegralBasisStructure {}

impl FactorableIdealsStructure for RingOfIntegersWithIntegralBasisStructure {
    fn factor_ideal(&self, ideal: &Self::Ideal) -> Option<DedekindDomainIdealFactorization<Self>> {
        Some(
            RingOfIntegersExtension::new_integer_extension(self.clone())
                .factor_ideal(ideal)?
                .into_full_factorization(),
        )
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

        {
            // 1 + sqrt(2)
            let alpha = roi.try_anf_to_roi(&(&x + 1).into_verbose()).unwrap();

            // (a + b sqrt(2)) * (1 + sqrt(2)) = a(1 + sqrt(2)) + b(2 + sqrt(2))
            assert!(roi.ideal_equal(
                &roi.principal_ideal(&alpha),
                &roi.ideal_from_integer_basis(vec![
                    roi.try_anf_to_roi(&(1 + &x).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(2 + &x).into_verbose()).unwrap()
                ])
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
                &roi.ideal_from_integer_basis(vec![
                    roi.try_anf_to_roi(&(3 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(3 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));

            // intersection is 30
            assert!(roi.ideal_equal(
                &alpha_beta_intersect,
                &roi.ideal_from_integer_basis(vec![
                    roi.try_anf_to_roi(&(30 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(30 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));

            // product is 90
            assert!(roi.ideal_equal(
                &alpha_beta_mul,
                &roi.ideal_from_integer_basis(vec![
                    roi.try_anf_to_roi(&(90 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(90 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));
        }
    }

    #[test]
    fn test_count_all_ideals_norm_eq() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1).into_verbose().algebraic_number_field();
        let roi = anf.ring_of_integers();

        assert_eq!(
            roi.all_ideals_norm_eq(&Natural::from(5040 as u32))
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            roi.all_ideals_norm_eq(&Natural::from(5040 * 7 as u32))
                .collect::<Vec<_>>()
                .len(),
            2
        );
    }
}
