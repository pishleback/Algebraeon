use super::number_field::*;
use super::ring_of_integers::*;
use crate::polynomial::PolynomialStructure;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct RingOfIntegersSquare {
    z: IntegerCanonicalStructure,
    q: RationalCanonicalStructure,
    r: RingOfIntegersWithIntegralBasisStructure,
    k: AlgebraicNumberFieldStructure,
    z_to_q: PrincipalSubringInclusion<RationalCanonicalStructure>,
    z_to_r: PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
    q_to_k: PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
    r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion,
}

impl RingOfIntegersSquare {
    pub fn new(roi: &RingOfIntegersWithIntegralBasisStructure) -> Self {
        let anf = roi.anf();
        RingOfIntegersSquare {
            z: Integer::structure(),
            q: Rational::structure(),
            r: roi.clone(),
            k: anf.clone(),
            z_to_q: PrincipalSubringInclusion::new(Rational::structure()),
            z_to_r: PrincipalSubringInclusion::new(roi.clone()),
            q_to_k: PrincipalRationalSubfieldInclusion::new(anf.clone()),
            r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(
                roi.clone(),
            ),
        }
    }
}

impl
    IntegralClosureSquare<
        IntegerCanonicalStructure,
        RationalCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldStructure,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > for RingOfIntegersSquare
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

    fn integralize_multiplier(
        &self,
        alpha: &<AlgebraicNumberFieldStructure as SetStructure>::Set,
    ) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            let q_poly = PolynomialStructure::new(self.q_field().clone());
            let alpha_min_poly_monic = self.q_to_k().min_poly(alpha);
            debug_assert!(q_poly.is_monic(&alpha_min_poly_monic));
            let alpha_min_poly_monic_coeffs = alpha_min_poly_monic.into_coeffs();
            let alpha_min_poly_monic_coeffs_denominators = alpha_min_poly_monic_coeffs
                .into_iter()
                .map(|c| self.z_to_q().denominator(&c))
                .collect();
            Integer::lcm_list(alpha_min_poly_monic_coeffs_denominators)
        }
    }
}

impl
    FactorablePrimeIdealsSquare<
        IntegerCanonicalStructure,
        RationalCanonicalStructure,
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldStructure,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > for RingOfIntegersSquare
{
    fn factor_prime_ideal(
        &self,
        p: &DedekindDomainPrimeIdeal<IntegerCanonicalStructure>,
    ) -> DedekindDomainIdealFactorization<RingOfIntegersWithIntegralBasisStructure> {
        // https://en.wikipedia.org/wiki/Dedekind%E2%80%93Kummer_theorem
        let p = Integer::ideal_generator(p.ideal());
        let anf = self.k_field();
        let roi = self.r_ring();
        let mod_p = QuotientStructure::new_field_unchecked(Integer::structure(), p.clone());
        let poly_mod_p = PolynomialStructure::new(mod_p);
        let poly_roi = PolynomialStructure::new(roi.clone());

        // alpha generates the algebraic number field but it is not necessarily an algebraic integer
        let alpha = anf.generator();
        // beta generates the algebraic number field and belongs to the ring of integers
        let beta = self.integral_scalar_multiple(&alpha);
        // factor the minimal polynomial of beta over the integers modulo p
        let beta_min_poly = self.min_poly_r_over_z(&beta);
        let beta_min_poly_factored = poly_mod_p.factor(&beta_min_poly).unwrap();
        // there is one prime ideal factor for each irreducible factor of beta's minimal polynomial modulo p
        // the prime ideal coresponding to an irreducible factor g(x) is generated by (p, g(beta))
        DedekindDomainIdealFactorization::from_powers_unchecked(
            beta_min_poly_factored
                .into_factors()
                .into_iter()
                .map(|(g, power)| {
                    debug_assert!(g.is_monic());
                    (
                        DedekindDomainPrimeIdeal::from_ideal_unchecked(roi.generated_ideal(vec![
                            self.z_to_r().image(&p),
                            poly_roi.evaluate(&g.apply_map(|c| self.z_to_r().image(&c)), &beta),
                        ])),
                        power,
                    )
                })
                .collect(),
        )
    }

    fn residue_class_degree(
        &self,
        p: &DedekindDomainPrimeIdeal<IntegerCanonicalStructure>,
        q: &DedekindDomainPrimeIdeal<RingOfIntegersWithIntegralBasisStructure>,
    ) -> Natural {
        todo!()
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
        let sq = RingOfIntegersSquare::new(&roi);
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
    fn factor_integer_prime_in_ring_of_ints() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1).into_verbose().algebraic_number_field();
        let roi = anf.ring_of_integers();
        let sq = RingOfIntegersSquare::new(&roi);

        let f2 =
            sq.factor_prime_ideal(&DedekindDomainPrimeIdeal::try_from_nat(2u32.into()).unwrap());
        assert!(!f2.is_inert());
        assert!(f2.is_ramified());

        let f3 =
            sq.factor_prime_ideal(&DedekindDomainPrimeIdeal::try_from_nat(3u32.into()).unwrap());
        assert!(f3.is_inert());
        assert!(!f3.is_ramified());

        let f5 =
            sq.factor_prime_ideal(&DedekindDomainPrimeIdeal::try_from_nat(5u32.into()).unwrap());
        assert!(!f5.is_inert());
        assert!(!f5.is_ramified());
    }
}
