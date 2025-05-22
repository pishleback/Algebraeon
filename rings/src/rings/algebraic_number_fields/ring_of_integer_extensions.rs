use super::ideal::RingOfIntegersIdeal;
use super::number_field::*;
use super::ring_of_integers::*;
use crate::linear::matrix::Matrix;
use crate::polynomial::Polynomial;
use crate::polynomial::PolynomialStructure;
use crate::rings::quotient::QuotientStructure;
use crate::rings::valuation::Valuation;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct RingOfIntegersExtension {
    z: IntegerCanonicalStructure,
    q: RationalCanonicalStructure,
    r: RingOfIntegersWithIntegralBasisStructure,
    k: AlgebraicNumberFieldStructure,
    z_to_q: PrincipalSubringInclusion<RationalCanonicalStructure>,
    z_to_r: PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
    q_to_k: PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
    r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion,
}

impl RingOfIntegersExtension {
    pub fn new_integer_extension(roi: RingOfIntegersWithIntegralBasisStructure) -> Self {
        let anf = roi.anf();
        RingOfIntegersExtension {
            z: Integer::structure(),
            q: Rational::structure(),
            r: roi.clone(),
            k: anf.clone(),
            z_to_q: PrincipalSubringInclusion::new(Rational::structure()),
            z_to_r: PrincipalSubringInclusion::new(roi.clone()),
            q_to_k: PrincipalRationalSubfieldInclusion::new(anf.clone()),
            r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(roi),
        }
    }

    pub fn padic_anf_valuation_element(
        &self,
        prime_ideal: RingOfIntegersIdeal,
        a: &Polynomial<Rational>,
    ) -> Valuation {
        let d = self.integralize_multiplier(a);
        let m = self.k_field().mul(a, &self.z_to_k().image(&d));
        self.r
            .padic_roi_element_valuation(prime_ideal.clone(), self.r.try_anf_to_roi(&m).unwrap())
            - self
                .r
                .padic_roi_element_valuation(prime_ideal, self.r.from_int(d))
    }
}

impl
    IntegralClosureExtension<
        IntegerCanonicalStructure,
        RationalCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldStructure,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > for RingOfIntegersExtension
{
    fn z_ring(&self) -> &IntegerCanonicalStructure {
        &self.z
    }
    fn r_ring(&self) -> &RingOfIntegersWithIntegralBasisStructure {
        &self.r
    }
    fn q_field(&self) -> &RationalCanonicalStructure {
        &self.q
    }
    fn k_field(&self) -> &AlgebraicNumberFieldStructure {
        &self.k
    }
    fn z_to_q(&self) -> &PrincipalSubringInclusion<RationalCanonicalStructure> {
        &self.z_to_q
    }
    fn z_to_r(&self) -> &PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure> {
        &self.z_to_r
    }
    fn q_to_k(&self) -> &PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure> {
        &self.q_to_k
    }
    fn r_to_k(&self) -> &RingOfIntegersToAlgebraicNumberFieldInclusion {
        &self.r_to_k
    }

    fn integralize_multiplier(&self, alpha: &Polynomial<Rational>) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            self.k_field().denominator(alpha)
        }
    }
}

impl
    DedekindDomainExtension<
        IntegerCanonicalStructure,
        RationalCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldStructure,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > for RingOfIntegersExtension
{
    fn ideal_norm(&self, ideal: &RingOfIntegersIdeal) -> Natural {
        #[cfg(debug_assertions)]
        self.r_ring().check_ideal(ideal);
        match ideal {
            RingOfIntegersIdeal::Zero => Natural::ZERO,
            RingOfIntegersIdeal::NonZero { lattice } => {
                let n = self.r_ring().degree();
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

    fn factor_prime_ideal(
        &self,
        prime_ideal: DedekindDomainPrimeIdeal<IntegerCanonicalStructure>,
    ) -> DedekindExtensionIdealFactorsAbovePrime<
        IntegerCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
    > {
        // https://en.wikipedia.org/wiki/Dedekind%E2%80%93Kummer_theorem
        let p = Integer::ideal_generator(prime_ideal.ideal());
        let anf = self.k_field();
        let roi = self.r_ring();
        let mod_p = QuotientStructure::new_field_unchecked(Integer::structure(), p.clone());
        let poly_mod_p = PolynomialStructure::new(mod_p);
        let poly_roi = PolynomialStructure::new(roi.clone());

        // alpha generates the algebraic number field but it is not necessarily an algebraic integer
        let alpha = anf.generator();
        // beta generates the algebraic number field and belongs to the ring of integers
        let beta = self.integral_scalar_multiple_r(&alpha);
        // factor the minimal polynomial of beta over the integers modulo p
        let beta_min_poly = self.min_poly_r_over_z(&beta);
        let beta_min_poly_factored = poly_mod_p.factor(&beta_min_poly).unwrap();
        // there is one prime ideal factor for each irreducible factor of beta's minimal polynomial modulo p
        // the prime ideal coresponding to an irreducible factor g(x) is generated by (p, g(beta))
        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
            roi.clone(),
            prime_ideal,
            beta_min_poly_factored
                .into_factor_powers()
                .into_iter()
                .map(|(g, power)| {
                    debug_assert!(g.is_monic());
                    let prime_ideal = roi.generated_ideal(vec![
                        self.z_to_r().image(&p),
                        poly_roi.evaluate(&g.apply_map(|c| self.z_to_r().image(&c)), &beta),
                    ]);
                    // norm(I) = p^deg(g)
                    debug_assert_eq!(
                        roi.ideal_norm(&prime_ideal),
                        p.clone().abs().nat_pow(&g.degree().unwrap().into())
                    );
                    DedekindExtensionIdealFactorsAbovePrimeFactor {
                        prime_ideal: DedekindDomainPrimeIdeal::from_ideal_unchecked(prime_ideal),
                        residue_class_degree: g.degree().unwrap(),
                        power,
                    }
                })
                .collect(),
        )
    }

    fn factor_ideal(
        &self,
        ideal: &RingOfIntegersIdeal,
    ) -> Option<
        DedekindExtensionIdealFactorization<
            IntegerCanonicalStructure,
            RingOfIntegersWithIntegralBasisStructure,
        >,
    > {
        let roi = self.r_ring();
        let extension_square = RingOfIntegersExtension::new_integer_extension(roi.clone());
        let norm = extension_square.ideal_norm(ideal);
        let norm_prime_factors = Integer::factor_ideal(&norm)?;
        Some(
            DedekindExtensionIdealFactorization::from_ideal_factors_above_primes(
                roi.clone(),
                norm_prime_factors
                    .into_squarefree_factor_list()
                    .into_iter()
                    .map(|prime| {
                        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
                            roi.clone(),
                            prime.clone(),
                            extension_square
                                .factor_prime_ideal(prime.clone())
                                .into_factors()
                                .into_iter()
                                .filter_map(|factor_above_prime| {
                                    let k = roi.largest_prime_ideal_factor_power(
                                        &factor_above_prime.prime_ideal,
                                        ideal,
                                    );
                                    if k == Natural::ZERO {
                                        None
                                    } else {
                                        Some(DedekindExtensionIdealFactorsAbovePrimeFactor {
                                            prime_ideal: factor_above_prime.prime_ideal,
                                            residue_class_degree: factor_above_prime
                                                .residue_class_degree,
                                            power: k,
                                        })
                                    }
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            ),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::Polynomial;

    #[test]
    fn integral_multiplier() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(3) + x + 1).into_verbose().algebraic_number_field();
        let roi = anf.ring_of_integers();
        let sq = RingOfIntegersExtension::new_integer_extension(roi.clone());
        let r_to_k_fof = sq.r_to_k_field_of_fractions();

        let sample_rats = Rational::exhaustive_rationals()
            .take(10)
            .collect::<Vec<_>>();

        for x in &sample_rats {
            for y in &sample_rats {
                for z in &sample_rats {
                    println!();
                    let alpha =
                        Polynomial::<Rational>::from_coeffs(vec![x.clone(), y.clone(), z.clone()]);
                    println!("alpha = {:?}", alpha);
                    let m = sq.integralize_multiplier(&alpha);
                    println!("m = {:?}", m);
                    let m_times_alpha = anf.mul(&anf.from_int(&m), &alpha);
                    println!("m * alpha = {:?}", m_times_alpha);
                    assert!(anf.is_algebraic_integer(&m_times_alpha));
                    let (n, d) = r_to_k_fof.numerator_and_denominator(&alpha);
                    println!("n={:?} d={:?}", n, d);
                    let (n, d) = (r_to_k_fof.image(&n), r_to_k_fof.image(&d));
                    println!("n={:?} d={:?}", n, d);
                    assert!(anf.equal(&alpha, &anf.div(&n, &d).unwrap()));
                }
            }
        }
    }

    #[test]
    fn ideals_opps_roi_extension_for_gaussian_integers() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1).into_verbose().algebraic_number_field();
        let roi = anf.ring_of_integers();
        let sq = RingOfIntegersExtension::new_integer_extension(roi.clone());

        let f2 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(2u32.into()).unwrap());
        assert!(!f2.is_inert());
        assert!(f2.is_ramified());
        println!("f2 = {:?}", f2);
        for ideal in f2.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 2u32.into())
        }

        let f3 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(3u32.into()).unwrap());
        assert!(f3.is_inert());
        assert!(!f3.is_ramified());
        println!("f3 = {:?}", f3);
        for ideal in f3.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 9u32.into())
        }

        let f5 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(5u32.into()).unwrap());
        assert!(!f5.is_inert());
        assert!(!f5.is_ramified());
        println!("f5 = {:?}", f5);
        for ideal in f5.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 5u32.into())
        }
    }
}
