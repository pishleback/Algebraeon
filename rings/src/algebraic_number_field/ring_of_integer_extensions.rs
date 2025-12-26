use super::ideal::RingOfIntegersIdeal;
use super::ideal::RingOfIntegersIdealsStructure;
use super::integer_lattice_ring_of_integers::*;
use super::polynomial_quotient_number_field::*;
use super::structure::AlgebraicNumberFieldSignature;
use crate::integer::ideal::IntegerIdealsStructure;
use crate::polynomial::RingToPolynomialSignature;
use crate::structure::*;
use algebraeon_nzq::traits::Abs;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Cow;
use std::marker::PhantomData;

/// Q -> K
/// ↑    ↑
/// Z -> R
///
/// Where Q is the rationals, Z is the integers, K is an algebraic number field, R is its ring of integers
///
#[derive(Debug, Clone)]
pub struct RingOfIntegersExtension<
    ANF: AlgebraicNumberFieldSignature,
    ROItoANF: BorrowedMorphism<ANF::RingOfIntegers, ANF, ANF::RingOfIntegersInclusion>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RB: BorrowedStructure<ANF::RingOfIntegers>,
    IdealsR: DedekindDomainIdealsSignature<ANF::RingOfIntegers, RB>,
> {
    _anf: PhantomData<ANF>,
    _roi: PhantomData<RB>,
    r_to_k: ROItoANF,
    ideals_z: IdealsZ,
    ideals_r: IdealsR,
}

impl<
    ANF: AlgebraicNumberFieldSignature,
    ROItoANF: BorrowedMorphism<ANF::RingOfIntegers, ANF, ANF::RingOfIntegersInclusion>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RB: BorrowedStructure<ANF::RingOfIntegers>,
    IdealsR: DedekindDomainIdealsSignature<ANF::RingOfIntegers, RB>,
> RingOfIntegersExtension<ANF, ROItoANF, IdealsZ, RB, IdealsR>
{
    pub fn new_integer_extension(r_to_k: ROItoANF, ideals_z: IdealsZ, ideals_r: IdealsR) -> Self {
        Self {
            _anf: PhantomData,
            _roi: PhantomData,
            r_to_k,
            ideals_z,
            ideals_r,
        }
    }
}

impl<
    ANF: AlgebraicNumberFieldSignature,
    ROItoANF: BorrowedMorphism<ANF::RingOfIntegers, ANF, ANF::RingOfIntegersInclusion>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RB: BorrowedStructure<ANF::RingOfIntegers>,
    IdealsR: DedekindDomainIdealsSignature<ANF::RingOfIntegers, RB>,
> IntegralClosureExtension for RingOfIntegersExtension<ANF, ROItoANF, IdealsZ, RB, IdealsR>
{
    type QKBasis = ANF::Basis;
    type Z = IntegerCanonicalStructure;
    type Q = RationalCanonicalStructure;
    type R = ANF::RingOfIntegers;
    type K = ANF;
    type ZQ<BZ: BorrowedStructure<Self::Z>, BQ: BorrowedStructure<Self::Q>> =
        PrincipalSubringInclusion<RationalCanonicalStructure, RationalCanonicalStructure>;
    type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>> =
        PrincipalSubringInclusion<ANF::RingOfIntegers, BR>;
    type QK<BQ: BorrowedStructure<Self::Q>, BK: BorrowedStructure<Self::K>> =
        ANF::RationalInclusion<BK>;
    type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>> =
        ANF::RingOfIntegersInclusion;

    fn z_ring(&self) -> &Self::Z {
        Integer::structure_ref()
    }
    fn r_ring(&self) -> &Self::R {
        self.r_to_k.borrow().domain()
    }
    fn q_field(&self) -> &Self::Q {
        Rational::structure_ref()
    }
    fn k_field(&self) -> &Self::K {
        self.r_to_k.borrow().range()
    }

    fn z_to_q<'a>(&'a self) -> Cow<'a, Self::ZQ<&'a Self::Z, &'a Self::Q>> {
        Cow::Owned(Rational::structure().into_principal_subring_inclusion())
    }
    fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>> {
        Cow::Owned(self.r_ring().principal_subring_inclusion())
    }
    fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>> {
        Cow::Owned(self.k_field().finite_dimensional_rational_extension())
    }
    fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>> {
        Cow::Borrowed(self.r_to_k.borrow())
    }

    fn integralize_multiplier(&self, alpha: &ANF::Set) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            self.k_field().min_poly_denominator_lcm(alpha)
        }
    }
}

impl<
    ROItoANF: BorrowedMorphism<
            RingOfIntegersWithIntegralBasisStructure,
            AlgebraicNumberFieldPolynomialQuotientStructure,
            RingOfIntegersToAlgebraicNumberFieldInclusion,
        >,
    RB: BorrowedStructure<RingOfIntegersWithIntegralBasisStructure>,
> DedekindDomainExtension<IntegerCanonicalStructure, RB>
    for RingOfIntegersExtension<
        AlgebraicNumberFieldPolynomialQuotientStructure,
        ROItoANF,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        RB,
        RingOfIntegersIdealsStructure<RB>,
    >
{
    type IdealsZ = IntegerIdealsStructure<IntegerCanonicalStructure>;
    type IdealsR = RingOfIntegersIdealsStructure<RB>;

    fn z_ideals(&self) -> &Self::IdealsZ {
        &self.ideals_z
    }

    fn r_ideals(&self) -> &Self::IdealsR {
        &self.ideals_r
    }

    fn ideal_norm(&self, ideal: &RingOfIntegersIdeal) -> Natural {
        self.r_ideals().ideal_norm(ideal)
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
        let mod_p = Integer::structure().into_quotient_field_unchecked(p.clone());
        let poly_mod_p = mod_p.polynomial_ring();
        let poly_roi = roi.polynomial_ring();

        // alpha generates the algebraic number field but it is not necessarily an algebraic integer
        let alpha = anf.generator();
        // beta generates the algebraic number field and belongs to the ring of integers
        let beta = self.integral_scalar_multiple_r(&alpha);
        // factor the minimal polynomial of beta over the integers modulo p
        let beta_min_poly = self.min_poly_r_over_z(&beta);
        let beta_min_poly_factored = poly_mod_p.factor(&beta_min_poly).unwrap();
        // there is one prime ideal factor for each irreducible factor of beta's minimal polynomial modulo p
        // the prime ideal corresponding to an irreducible factor g(x) is generated by (p, g(beta))
        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
            prime_ideal,
            poly_mod_p
                .factorizations()
                .into_powers(beta_min_poly_factored)
                .into_iter()
                .map(|(g, power)| {
                    debug_assert!(g.is_monic());
                    let prime_ideal = roi_ideals.generated_ideal(vec![
                        self.z_to_r().image(&p),
                        poly_roi.evaluate(&g.apply_map(|c| self.z_to_r().image(c)), &beta),
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
        let norm = self.ideal_norm(ideal);
        let norm_prime_factors = Integer::ideals().factor_ideal(&norm)?;
        Some(
            DedekindExtensionIdealFactorization::from_ideal_factors_above_primes(
                Integer::ideals()
                    .factorizations()
                    .into_prime_support(norm_prime_factors)
                    .into_iter()
                    .map(|prime| {
                        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
                            prime.clone(),
                            self.factor_prime_ideal(prime.clone())
                                .into_factors()
                                .into_iter()
                                .filter_map(|factor_above_prime| {
                                    let k = self.r_ideals().largest_prime_ideal_factor_power(
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
    use crate::{
        algebraic_number_field::structure::AlgebraicIntegerRingInAlgebraicNumberField,
        polynomial::Polynomial,
    };

    #[test]
    fn integral_multiplier() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(3) + x + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi_to_anf = anf.clone().into_ring_of_integers_extension();
        let sq = roi_to_anf.zq_extension();

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
    fn ideals_operations_roi_extension_for_gaussian_integers() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi_to_anf = anf.clone().into_ring_of_integers_extension();
        let sq = roi_to_anf.zq_extension();

        let f2 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(2u32.into()).unwrap());
        assert!(!f2.is_inert());
        assert!(f2.is_ramified());
        println!("f2 = {:?}", f2);
        for ideal in f2.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 2u32.into());
        }

        let f3 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(3u32.into()).unwrap());
        assert!(f3.is_inert());
        assert!(!f3.is_ramified());
        println!("f3 = {:?}", f3);
        for ideal in f3.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 9u32.into());
        }

        let f5 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(5u32.into()).unwrap());
        assert!(!f5.is_inert());
        assert!(!f5.is_ramified());
        println!("f5 = {:?}", f5);
        for ideal in f5.unique_prime_factors() {
            assert_eq!(sq.ideal_norm(&ideal.clone().into_ideal()), 5u32.into());
        }
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_is_S_integral() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 5)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi_to_anf = anf.clone().into_ring_of_integers_extension();
        let sq = roi_to_anf.zq_extension();

        // Element: (1/2) + sqrt(-5)
        let poly = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE_HALF, Rational::ONE]);

        // primes above (2)
        let f2 =
            sq.factor_prime_ideal(DedekindDomainPrimeIdeal::try_from_nat(2u32.into()).unwrap());
        let prime_ideals_above_2 = f2.unique_prime_factors();

        let d = sq.integralize_multiplier(&poly);
        let m = sq.k_field().mul(&poly, &sq.z_to_k().image(&d));
        println!("poly = {:?}", poly);
        println!("d = {:?}", d);
        println!("m = {:?}", m);
        assert!(
            sq.r_ring().try_anf_to_roi(&m).is_some(),
            "m not in ROI: {:?}",
            m
        );

        // Case 1: S = empty set → should NOT be S-integral, denominator 2 not inverted
        assert!(!sq.is_S_integral(vec![], &poly));
        // Case 2: S = prime ideals above 2 → now we allow inversion of 2
        assert!(sq.is_S_integral(prime_ideals_above_2, &poly));
    }
}
