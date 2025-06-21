use super::integer_lattice_ring_of_integers::*;
use crate::{
    algebraic_number_field::structure::AlgebraicIntegerRingInAlgebraicNumberField,
    matrix::Matrix,
    module::{
        finitely_free_affine::FinitelyFreeSubmoduleAffineSubset,
        finitely_free_submodule::FinitelyFreeSubmodule,
    },
    structure::*,
};
use algebraeon_nzq::{Integer, Natural, traits::Abs};
use algebraeon_sets::{
    combinatorics::num_partitions_part_pool,
    structure::{BorrowedStructure, EqSignature, MetaType, SetSignature, Signature},
};
use itertools::Itertools;

#[derive(Debug, Clone)]
pub enum RingOfIntegersIdeal {
    Zero,
    NonZero(FinitelyFreeSubmodule<Integer>),
}

impl RingOfIntegersIdeal {
    /// A basis of this ideal as a Z-module.
    pub fn basis(&self) -> Option<Vec<Vec<Integer>>> {
        match self {
            RingOfIntegersIdeal::Zero => None,
            RingOfIntegersIdeal::NonZero(lattice) => Some(lattice.basis()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingOfIntegersIdealsStructure<
    RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>,
> {
    roi: RingB,
}

impl RingToIdealsSignature for RingOfIntegersWithIntegralBasisStructure {
    type Ideals<SelfB: BorrowedStructure<Self>> = RingOfIntegersIdealsStructure<SelfB>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        RingOfIntegersIdealsStructure { roi: self }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        RingOfIntegersIdealsStructure { roi: self }
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>> Signature
    for RingOfIntegersIdealsStructure<RingB>
{
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>> SetSignature
    for RingOfIntegersIdealsStructure<RingB>
{
    type Set = RingOfIntegersIdeal;

    fn is_element(&self, ideal: &Self::Set) -> Result<(), String> {
        match ideal {
            RingOfIntegersIdeal::Zero => Ok(()),
            RingOfIntegersIdeal::NonZero(lattice) => {
                // check it's a submodule
                self.ring().z_module().submodules().is_element(lattice)?;
                // check it's an ideal
                for ideal_basis_elem in lattice.basis() {
                    for integral_basis_elem in
                        (0..self.ring().degree()).map(|i| self.ring().z_module().basis_element(i))
                    {
                        let x = self.ring().mul(&ideal_basis_elem, &integral_basis_elem);
                        if !self
                            .ring()
                            .z_module()
                            .submodules()
                            .contains_element(lattice, &x)
                        {
                            return Err("submodule is not an ideal".to_string());
                        }
                    }
                }
                Ok(())
            }
        }
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    IdealsSignature<RingOfIntegersWithIntegralBasisStructure, RingB>
    for RingOfIntegersIdealsStructure<RingB>
{
    fn ring(&self) -> &RingOfIntegersWithIntegralBasisStructure {
        self.roi.borrow()
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    RingOfIntegersIdealsStructure<RingB>
{
    /// Construct an ideal from a Z-linear span
    pub fn ideal_from_integer_span(&self, span: Vec<Vec<Integer>>) -> RingOfIntegersIdeal {
        for elem in &span {
            debug_assert!(self.ring().is_element(elem).is_ok());
        }
        let n = self.ring().degree();
        RingOfIntegersIdeal::NonZero(
            Matrix::join_cols(
                n,
                span.into_iter()
                    .map(|elem| Matrix::from_cols(vec![elem]))
                    .collect(),
            )
            .col_span(),
        )
    }

    pub fn ideal_norm(&self, ideal: &RingOfIntegersIdeal) -> Natural {
        debug_assert!(self.is_element(ideal).is_ok());
        match ideal {
            RingOfIntegersIdeal::Zero => Natural::ZERO,
            RingOfIntegersIdeal::NonZero(lattice) => {
                let n = self.ring().degree();
                let cols = lattice.basis();
                #[cfg(debug_assertions)]
                for col in &cols {
                    assert_eq!(col.len(), n);
                }
                let mat = Matrix::construct(n, n, |i, j| cols[i][j].clone());
                mat.det().unwrap().abs()
            }
        }
    }

    // Order of the multiplicative group of the quotient modulo the ideal.
    pub fn euler_phi(&self, ideal: &RingOfIntegersIdeal) -> Option<Natural> {
        match ideal {
            RingOfIntegersIdeal::Zero => None,
            RingOfIntegersIdeal::NonZero { .. } => Some(
                self.factorizations()
                    .into_powers(self.factor_ideal(ideal).unwrap())
                    .iter()
                    .map(|(prime_ideal, exponent)| {
                        let norm = self.ideal_norm(prime_ideal.ideal());
                        let e_minus_1 = exponent - Natural::ONE;
                        (&norm - Natural::ONE) * norm.pow(&e_minus_1)
                    })
                    .fold(Natural::ONE, |acc, x| acc * x),
            ),
        }
    }

    /// generate all ideals of norm equal to n
    pub fn all_ideals_norm_eq<'a>(
        &'a self,
        n: &Natural,
    ) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        match Integer::ideals().factor_ideal(n) {
            Some(n) => {
                let roi_to_anf =
                    RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(
                        self.ring().clone(),
                    );
                let sq = roi_to_anf.zq_extension();
                Box::new(
                    Integer::structure()
                        .ideals()
                        .factorizations()
                        .into_powers(n)
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
                .map(Natural::from)
                .take_while(|m| m <= n)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }

    /// generate all ideals
    pub fn all_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        Box::new(
            (0usize..)
                .map(Natural::from)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }

    /// generate all non-zero ideals
    pub fn all_nonzero_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = RingOfIntegersIdeal>> {
        Box::new(
            (1usize..)
                .map(Natural::from)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    IdealsArithmeticSignature<RingOfIntegersWithIntegralBasisStructure, RingB>
    for RingOfIntegersIdealsStructure<RingB>
{
    fn principal_ideal(&self, a: &Vec<Integer>) -> Self::Set {
        if self.ring().is_zero(a) {
            Self::Set::Zero
        } else {
            let n = self.ring().degree();
            let ideal = self.ideal_from_integer_span(
                (0..n)
                    .map(|i| {
                        self.ring()
                            .try_anf_to_roi(
                                &self
                                    .ring()
                                    .anf()
                                    .mul(self.ring().basis_element(i), &self.ring().roi_to_anf(a)),
                            )
                            .unwrap()
                    })
                    .collect(),
            );
            debug_assert!(self.is_element(&ideal).is_ok());
            ideal
        }
    }

    fn ideal_equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::Zero) => true,
            (RingOfIntegersIdeal::NonZero(a_lattice), RingOfIntegersIdeal::NonZero(b_lattice)) => {
                self.ring()
                    .z_module()
                    .submodules()
                    .equal(a_lattice, b_lattice)
            }
            _ => false,
        }
    }

    fn ideal_contains(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (_, RingOfIntegersIdeal::Zero) => true,
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::NonZero { .. }) => {
                debug_assert_ne!(self.ring().degree(), 0);
                false
            }
            (RingOfIntegersIdeal::NonZero(a_lattice), RingOfIntegersIdeal::NonZero(b_lattice)) => {
                self.ring()
                    .z_module()
                    .submodules()
                    .contains(a_lattice, b_lattice)
            }
        }
    }

    fn ideal_contains_element(&self, a: &Self::Set, x: &Vec<Integer>) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.ring().is_element(x).is_ok());
        match a {
            RingOfIntegersIdeal::Zero => self.ring().is_zero(x),
            RingOfIntegersIdeal::NonZero(lattice) => self
                .ring()
                .z_module()
                .submodules()
                .contains_element(lattice, x),
        }
    }

    fn ideal_intersect(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (RingOfIntegersIdeal::NonZero(a_lattice), RingOfIntegersIdeal::NonZero(b_lattice)) => {
                Self::Set::NonZero(
                    self.ring()
                        .z_module()
                        .submodules()
                        .intersect(a_lattice, b_lattice),
                )
            }
            _ => Self::Set::Zero,
        }
    }

    fn ideal_add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::Zero) => RingOfIntegersIdeal::Zero,
            (RingOfIntegersIdeal::Zero, RingOfIntegersIdeal::NonZero { .. }) => b.clone(),
            (RingOfIntegersIdeal::NonZero { .. }, RingOfIntegersIdeal::Zero) => a.clone(),
            (RingOfIntegersIdeal::NonZero(a_lattice), RingOfIntegersIdeal::NonZero(b_lattice)) => {
                Self::Set::NonZero(
                    self.ring()
                        .z_module()
                        .submodules()
                        .sum(a_lattice, b_lattice),
                )
            }
        }
    }

    fn ideal_mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (RingOfIntegersIdeal::NonZero(a_lattice), RingOfIntegersIdeal::NonZero(b_lattice)) => {
                let n = self.ring().degree();
                let a_basis = a_lattice.basis();
                let b_basis = b_lattice.basis();
                debug_assert_eq!(a_basis.len(), n);
                debug_assert_eq!(b_basis.len(), n);

                let mut span = vec![];
                for i in 0..n {
                    for j in 0..n {
                        span.push(self.ring().mul(&a_basis[i], &b_basis[j]));
                    }
                }
                self.ideal_from_integer_span(span)
            }
            _ => Self::Set::Zero,
        }
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    DedekindDomainIdealsSignature<RingOfIntegersWithIntegralBasisStructure, RingB>
    for RingOfIntegersIdealsStructure<RingB>
{
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    FactorableIdealsSignature<RingOfIntegersWithIntegralBasisStructure, RingB>
    for RingOfIntegersIdealsStructure<RingB>
{
    fn factor_ideal(
        &self,
        ideal: &Self::Set,
    ) -> Option<DedekindDomainIdealFactorization<Self::Set>> {
        Some(
            RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(
                self.ring().clone(),
            )
            .zq_extension()
            .factor_ideal(ideal)?
            .into_full_factorization(),
        )
    }
}

impl<RingB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>>
    RingOfIntegersIdealsStructure<RingB>
{
    /// given an ideal I and element a find an element b such that I = (a, b)
    pub fn ideal_other_generator(
        &self,
        g: &Vec<Integer>,
        ideal: &RingOfIntegersIdeal,
    ) -> Vec<Integer> {
        debug_assert!(self.ideal_contains_element(ideal, g));
        debug_assert!(!self.ring().is_zero(g));
        // prod_i p^{e_i}
        let ideal_factored = self.factor_ideal(ideal).unwrap();
        // prod_i p^{f_i} * prod_j q^{g_j}
        let g_factored = self.factor_ideal(&self.principal_ideal(g)).unwrap();
        // want b not in any q and in all p^{e_i} and not in any p^{e_i+1}

        // this is all b not in any q and in all p^{e_i}
        let b_set = self.ring().z_module().affine_subsets().intersect_list(
            self.factorizations()
                .to_powers(&ideal_factored)
                .into_iter()
                .map(|(p, k)| match self.ideal_nat_pow(p.ideal(), k) {
                    RingOfIntegersIdeal::Zero => unreachable!(),
                    RingOfIntegersIdeal::NonZero(pk_lattice) => {
                        FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                            self.ring().z_module().cosets().from_submodule(pk_lattice),
                        )
                    }
                })
                .chain(
                    self.factorizations()
                        .into_prime_support(g_factored)
                        .into_iter()
                        .filter(|prime_ideal| {
                            !self
                                .factorizations()
                                .to_prime_support(&ideal_factored)
                                .into_iter()
                                .any(|p| self.ideal_equal(p.ideal(), prime_ideal.ideal()))
                        })
                        .map(|q| match q.into_ideal() {
                            RingOfIntegersIdeal::Zero => unreachable!(),
                            RingOfIntegersIdeal::NonZero(q_lattice) => {
                                FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                                    self.ring()
                                        .z_module()
                                        .cosets()
                                        .from_offset_and_submodule(&self.ring().one(), q_lattice),
                                )
                            }
                        }),
                )
                .collect(),
        );

        //need to filter out the b in some p^{e_i+1}
        let rm_b_set = self.ring().z_module().affine_subsets().intersect_list(
            self.factorizations()
                .to_powers(&ideal_factored)
                .into_iter()
                .map(
                    |(p, k)| match self.ideal_nat_pow(p.ideal(), &(k + Natural::ONE)) {
                        RingOfIntegersIdeal::Zero => unreachable!(),
                        RingOfIntegersIdeal::NonZero(pk_lattice) => {
                            FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                                self.ring().z_module().cosets().from_submodule(pk_lattice),
                            )
                        }
                    },
                )
                .collect(),
        );

        //if all basis elements of b_set were contain in rm_b_set then we'd have b_set contained in rm_b_set
        //but this is not the case, so some basis of b_set is not in rm_b_set

        self.ring()
            .z_module()
            .affine_subsets()
            .affine_basis(&b_set)
            .into_iter()
            .find(|b| {
                !self
                    .ring()
                    .z_module()
                    .affine_subsets()
                    .contains_element(&rm_b_set, b)
            })
            .unwrap()
    }

    /// return two elements which generate the ideal
    pub fn ideal_two_generators(
        &self,
        ideal: &RingOfIntegersIdeal,
    ) -> (Vec<Integer>, Vec<Integer>) {
        let (a, b) = match ideal {
            RingOfIntegersIdeal::Zero => (self.ring().zero(), self.ring().zero()),
            RingOfIntegersIdeal::NonZero(lattice) => {
                let a = lattice.basis().into_iter().next().unwrap();
                let b = self.ideal_other_generator(&a, ideal);
                (a, b)
            }
        };
        debug_assert!(self.ideal_equal(ideal, &self.generated_ideal(vec![a.clone(), b.clone()])));
        (a, b)
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
        let anf = (x.pow(2) - 2)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = RingOfIntegersWithIntegralBasisStructure::new(
            anf.clone(),
            vec![a.clone(), b.clone()],
            Integer::from(8),
        );
        let roi_ideals = roi.ideals();

        {
            // 1 + sqrt(2)
            let alpha = roi.try_anf_to_roi(&(&x + 1).into_verbose()).unwrap();

            // (a + b sqrt(2)) * (1 + sqrt(2)) = a(1 + sqrt(2)) + b(2 + sqrt(2))
            assert!(roi_ideals.ideal_equal(
                &roi_ideals.principal_ideal(&alpha),
                &roi_ideals.ideal_from_integer_span(vec![
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

            let alpha_ideal = roi_ideals.principal_ideal(&alpha);
            let beta_ideal = roi_ideals.principal_ideal(&beta);

            let alpha_beta_add = roi_ideals.ideal_add(&alpha_ideal, &beta_ideal);
            let alpha_beta_intersect = roi_ideals.ideal_intersect(&alpha_ideal, &beta_ideal);
            let alpha_beta_mul = roi_ideals.ideal_mul(&alpha_ideal, &beta_ideal);

            // sum is 3
            assert!(roi_ideals.ideal_equal(
                &alpha_beta_add,
                &roi_ideals.ideal_from_integer_span(vec![
                    roi.try_anf_to_roi(&(3 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(3 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));

            // intersection is 30
            assert!(roi_ideals.ideal_equal(
                &alpha_beta_intersect,
                &roi_ideals.ideal_from_integer_span(vec![
                    roi.try_anf_to_roi(&(30 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(30 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));

            // product is 90
            assert!(roi_ideals.ideal_equal(
                &alpha_beta_mul,
                &roi_ideals.ideal_from_integer_span(vec![
                    roi.try_anf_to_roi(&(90 * x.pow(0)).into_verbose()).unwrap(),
                    roi.try_anf_to_roi(&(90 * x.pow(1)).into_verbose()).unwrap()
                ])
            ));
        }
    }

    #[test]
    fn test_count_all_ideals_norm_eq() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.compute_ring_of_integers();
        let roi_ideals = roi.ideals();

        assert_eq!(
            roi_ideals
                .all_ideals_norm_eq(&Natural::from(5040 as u32))
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            roi_ideals
                .all_ideals_norm_eq(&Natural::from(5040 * 7 as u32))
                .collect::<Vec<_>>()
                .len(),
            2
        );
    }

    #[test]
    fn test_euler_phi_of_principal_ideal() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Construct the number field Q(i), which has ring of integers Z[i]
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.compute_ring_of_integers();
        let roi_ideals = roi.ideals();

        // Consider the ideal (5)
        let ideal = roi_ideals.principal_ideal(&roi.from_int(5));

        let phi = roi_ideals.euler_phi(&ideal).unwrap();
        assert_eq!(phi, Natural::from(16u32));
    }
}
