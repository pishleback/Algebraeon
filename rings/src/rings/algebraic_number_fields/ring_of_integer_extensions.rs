use super::ideal::RingOfIntegersIdeal;
use super::ideal::RingOfIntegersIdealsStructure;
use super::number_field::*;
use super::ring_of_integers::*;
use crate::linear::matrix::Matrix;
use crate::polynomial::Polynomial;
use crate::polynomial::PolynomialStructure;
use crate::rings::integer::ideal::IntegerIdealsStructure;
use crate::rings::quotient::QuotientStructure;
use crate::rings::valuation::Valuation;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct RingOfIntegersExtension<
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<
            RingOfIntegersWithIntegralBasisStructure,
            RingOfIntegersWithIntegralBasisStructure,
        >,
> {
    z: IntegerCanonicalStructure,
    q: RationalCanonicalStructure,
    r: RingOfIntegersWithIntegralBasisStructure,
    k: AlgebraicNumberFieldStructure,
    z_to_q: PrincipalSubringInclusion<RationalCanonicalStructure>,
    z_to_r: PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
    q_to_k: PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
    r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion,
    ideals_z: IdealsZ,
    ideals_r: IdealsR,
}

impl
    RingOfIntegersExtension<
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        RingOfIntegersIdealsStructure<RingOfIntegersWithIntegralBasisStructure>,
    >
{
    pub fn new_integer_extension(roi: RingOfIntegersWithIntegralBasisStructure) -> Self {
        let ideals_r = roi.clone().into_ideals();
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
            ideals_z: Integer::ideals(),
            ideals_r,
        }
    }

    pub fn padic_anf_valuation_element(
        &self,
        prime_ideal: RingOfIntegersIdeal,
        a: &Polynomial<Rational>,
    ) -> Valuation {
        let d = self.integralize_multiplier(a);
        let m = self.k_field().mul(a, &self.z_to_k().image(&d));
        self.ideals_r
            .padic_roi_element_valuation(prime_ideal.clone(), self.r.try_anf_to_roi(&m).unwrap())
            - self
                .ideals_r
                .padic_roi_element_valuation(prime_ideal, self.r.from_int(d))
    }

    // An element is S-integral, if its valuations at all primes not in S are nonnegative.
    // If S is the empty set, this coincides with the usual integrality.
    #[allow(non_snake_case)]
    pub fn is_S_integral(
        &self,
        S: Vec<&DedekindDomainPrimeIdeal<RingOfIntegersIdeal>>,
        a: &Polynomial<Rational>,
    ) -> bool {
        let d = self.integralize_multiplier(a);
        let m = self.k_field().mul(a, &self.z_to_k().image(&d));
        // for each prime factor P of d not in S, check if valuation_P(m) ≥ valuation_P(d)

        let d_as_roi = self.r.from_int(d.clone());
        let principal_ideal_d = self.ideals_r.generated_ideal(vec![d_as_roi.clone()]);

        let d_factorization = self.factor_ideal(&principal_ideal_d);
        if d_factorization.is_none() {
            return true;
        }

        for prime in d_factorization.unwrap().into_powers() {
            let prime_ideal = prime.0.into_ideal();
            // Skip primes in S
            if S.iter().any(|s_ideal| {
                self.ideals_r
                    .ideal_equal(&(*s_ideal).clone().into_ideal(), &prime_ideal)
            }) {
                continue;
            }

            let m_val = self.ideals_r.padic_roi_element_valuation(
                prime_ideal.clone(),
                self.r.try_anf_to_roi(&m).unwrap(),
            );
            let d_val = self
                .ideals_r
                .padic_roi_element_valuation(prime_ideal, d_as_roi.clone());

            if m_val < d_val {
                return false;
            }
        }

        true
    }
}

impl<
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<
            RingOfIntegersWithIntegralBasisStructure,
            RingOfIntegersWithIntegralBasisStructure,
        >,
>
    IntegralClosureExtension<
        IntegerCanonicalStructure,
        RationalCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldStructure,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > for RingOfIntegersExtension<IdealsZ, IdealsR>
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
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        RingOfIntegersIdealsStructure<RingOfIntegersWithIntegralBasisStructure>,
    >
    for RingOfIntegersExtension<
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        RingOfIntegersIdealsStructure<RingOfIntegersWithIntegralBasisStructure>,
    >
{
    fn z_ideals(&self) -> &IntegerIdealsStructure<IntegerCanonicalStructure> {
        &self.ideals_z
    }

    fn r_ideals(&self) -> &RingOfIntegersIdealsStructure<RingOfIntegersWithIntegralBasisStructure> {
        &self.ideals_r
    }

    fn ideal_norm(&self, ideal: &RingOfIntegersIdeal) -> Natural {
        debug_assert!(self.r_ideals().is_element(ideal));
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
        prime_ideal: DedekindDomainPrimeIdeal<Natural>,
    ) -> DedekindExtensionIdealFactorsAbovePrime<Natural, RingOfIntegersIdeal> {
        // https://en.wikipedia.org/wiki/Dedekind%E2%80%93Kummer_theorem
        let p = Integer::ideals().ideal_generator(prime_ideal.ideal());
        let anf = self.k_field();
        let roi = self.r_ring();
        let roi_ideals = roi.ideals();
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
            prime_ideal,
            poly_mod_p
                .factorizations()
                .into_factor_powers(beta_min_poly_factored)
                .into_iter()
                .map(|(g, power)| {
                    debug_assert!(g.is_monic());
                    let prime_ideal = roi_ideals.generated_ideal(vec![
                        self.z_to_r().image(&p),
                        poly_roi.evaluate(&g.apply_map(|c| self.z_to_r().image(&c)), &beta),
                    ]);
                    // norm(I) = p^deg(g)
                    debug_assert_eq!(
                        roi_ideals.ideal_norm(&prime_ideal),
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
    ) -> Option<DedekindExtensionIdealFactorization<Natural, RingOfIntegersIdeal>> {
        let roi = self.r_ring();
        let roi_ideals = roi.ideals();
        let extension_square = RingOfIntegersExtension::new_integer_extension(roi.clone());
        let norm = extension_square.ideal_norm(ideal);
        let norm_prime_factors = Integer::ideals().factor_ideal(&norm)?;
        Some(
            DedekindExtensionIdealFactorization::from_ideal_factors_above_primes(
                Integer::ideals()
                    .factorizations()
                    .into_squarefree_factor_list(norm_prime_factors)
                    .into_iter()
                    .map(|prime| {
                        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
                            prime.clone(),
                            extension_square
                                .factor_prime_ideal(prime.clone())
                                .into_factors()
                                .into_iter()
                                .filter_map(|factor_above_prime| {
                                    let k = roi_ideals.largest_prime_ideal_factor_power(
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

    #[test]
    #[allow(non_snake_case)]
    fn test_is_S_integral() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 5).into_verbose().algebraic_number_field();
        let roi = anf.ring_of_integers();
        let ext = RingOfIntegersExtension::new_integer_extension(roi.clone());

        // Element: (1/2) + sqrt(-5)
        let poly = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE_HALF, Rational::ONE]);

        // primes above (2)
        let f2 =
            ext.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(2u32.into()).unwrap());
        let prime_ideals_above_2 = f2.unique_prime_factors();

        let d = ext.integralize_multiplier(&poly);
        let m = ext.k_field().mul(&poly, &ext.z_to_k().image(&d));
        println!("poly = {:?}", poly);
        println!("d = {:?}", d);
        println!("m = {:?}", m);
        assert!(ext.r.try_anf_to_roi(&m).is_some(), "m not in ROI: {:?}", m);

        // Case 1: S = empty set → should NOT be S-integral, denominator 2 not inverted
        assert!(!ext.is_S_integral(vec![], &poly));
        // Case 2: S = prime ideals above 2 → now we allow inversion of 2
        assert!(ext.is_S_integral(prime_ideals_above_2, &poly));
    }
}
