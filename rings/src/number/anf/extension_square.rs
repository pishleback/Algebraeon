use super::number_field::*;
use super::ring_of_integers::*;
use crate::polynomial::PolynomialStructure;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
struct RingOfIntegersSquare {
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
    fn new(roi: &RingOfIntegersWithIntegralBasisStructure) -> Self {
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
        _p: &DedekindDomainPrimeIdeal<IntegerCanonicalStructure>,
    ) -> DedekindDomainIdealFactorization<RingOfIntegersWithIntegralBasisStructure> {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::Polynomial;

    #[test]
    fn ring_of_integers_integral_multiplier_test() {
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
}
