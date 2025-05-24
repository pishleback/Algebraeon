use super::*;
use crate::polynomial::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
struct FieldOfFractionsInclusionForIntegralClosure<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> {
    z: PhantomData<Z>,
    q: PhantomData<Q>,
    r: PhantomData<R>,
    k: PhantomData<K>,
    zq: PhantomData<ZQ>,
    zr: PhantomData<ZR>,
    qk: PhantomData<QK>,
    rk: PhantomData<RK>,
    square: ICS,
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn new(square: ICS) -> Self {
        Self {
            z: PhantomData,
            q: PhantomData,
            r: PhantomData,
            k: PhantomData,
            zq: PhantomData,
            zr: PhantomData,
            qk: PhantomData,
            rk: PhantomData,
            square,
        }
    }
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> Morphism<R, K> for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn domain(&self) -> &R {
        self.square.r_ring()
    }

    fn range(&self) -> &K {
        self.square.k_field()
    }
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> Function<R, K> for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn image(&self, x: &R::Set) -> K::Set {
        self.square.r_to_k().image(x)
    }
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> InjectiveFunction<R, K>
    for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn try_preimage(&self, x: &K::Set) -> Option<R::Set> {
        self.square.r_to_k().try_preimage(x)
    }
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> RingHomomorphism<R, K>
    for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
}

impl<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>,
> FieldOfFractionsInclusion<R, K>
    for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn numerator_and_denominator(&self, a: &K::Set) -> (R::Set, R::Set) {
        // let d in Z such that d*a is in R
        let d = self.square.integralize_multiplier(a);
        // take d in R
        let d = self.square.z_to_r().image(&d);
        // now a = (d*a) / d
        let n = self
            .try_preimage(&self.range().mul(&self.image(&d), a))
            .unwrap();
        debug_assert!(self.range().equal(
            a,
            &self.range().div(&self.image(&n), &self.image(&d)).unwrap()
        ));
        (n, d)
    }
}

/// Given a commuting square of injective ring homomorphisms
///
/// Q → K
/// ↑   ↑
/// Z → R
///
/// such that
///  - Q is the field of fractions of Z
///  - Q → K is a finite dimensional field extension
///
/// This trait expresses that R is the integral closure of Z in K
pub trait IntegralClosureExtension<
    Z: IntegralDomainSignature,
    Q: FieldSignature,
    R: IntegralDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
>: Debug + Clone
{
    fn z_ring(&self) -> &Z;
    fn q_field(&self) -> &Q;
    fn r_ring(&self) -> &R;
    fn k_field(&self) -> &K;

    fn z_to_q(&self) -> &ZQ;
    fn z_to_r(&self) -> &ZR;
    fn q_to_k(&self) -> &QK;
    fn r_to_k(&self) -> &RK;

    /// The square should commute, so this should be both
    /// - z_to_q followed by q_to_k
    /// - z_to_r followed by r_to_k
    fn z_to_k(&self) -> impl RingHomomorphism<Z, K> + InjectiveFunction<Z, K> {
        CompositionMorphism::new(self.z_to_q().clone(), self.q_to_k().clone())
    }

    /// The monic minimal polynomial of alpha in K over Q
    fn min_poly_k_over_q(&self, alpha: &K::Set) -> Polynomial<Q::Set> {
        let alpha_min_poly_monic = self.q_to_k().min_poly(&alpha);
        #[cfg(debug_assertions)]
        {
            let q_poly = PolynomialStructure::new(self.q_field().clone());
            assert!(q_poly.is_monic(&alpha_min_poly_monic));
        }
        alpha_min_poly_monic
    }

    /// By definition of R as the integral closure of Z in K every element of R, when considered as an element of K, has minimal polynomial over Q which is monic with coefficients in Z
    fn min_poly_r_over_z(&self, alpha: &R::Set) -> Polynomial<Z::Set> {
        let alpha_min_poly_monic = self
            .min_poly_k_over_q(&self.r_to_k().image(alpha))
            .apply_map_into(|c| self.z_to_q().try_preimage(&c).unwrap());
        #[cfg(debug_assertions)]
        {
            let z_poly = PolynomialStructure::new(self.z_ring().clone());
            assert!(z_poly.is_monic(&alpha_min_poly_monic));
        }
        alpha_min_poly_monic
    }

    /// For alpha in K return non-zero d in Z such that d*alpha is in R
    fn integralize_multiplier(&self, alpha: &K::Set) -> Z::Set;

    fn integral_scalar_multiple_r(&self, alpha: &K::Set) -> R::Set {
        self.r_to_k()
            .try_preimage(&self.integral_scalar_multiple_k(alpha))
            .unwrap()
    }

    fn integral_scalar_multiple_k(&self, alpha: &K::Set) -> K::Set {
        let d = self.integralize_multiplier(&alpha);
        // This is the LCM of the denominators of the coefficients of a,
        // and thus it may well be > 1 even when the element is an algebraic integer.
        self.k_field()
            .mul(&self.r_to_k().image(&self.z_to_r().image(&d)), alpha)
    }

    /// Every element of K is a fraction of elements of R
    fn r_to_k_field_of_fractions(&self) -> impl FieldOfFractionsInclusion<R, K> {
        FieldOfFractionsInclusionForIntegralClosure::new(self.clone())
    }
}

/// A commuting square of injective ring homomorphisms
///
/// Q → K
/// ↑   ↑
/// Z → R
///
/// such that
///  - Q is the field of fractions of Z
///  - Q → K is a finite dimensional field extension
///  - R is the integral closure of Z in K
///  - Z and R are Dedekind domains
///
/// This trait allows, for each prime ideal p of Z, the ideal pR of R to be factored into prime ideals in R
pub trait DedekindDomainExtension<
    Z: DedekindDomainSignature,
    Q: FieldSignature,
    R: DedekindDomainSignature,
    K: FieldSignature,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    IdealsZ: DedekindDomainIdealsSignature<Z>,
    IdealsR: DedekindDomainIdealsSignature<R>,
>: IntegralClosureExtension<Z, Q, R, K, ZQ, ZR, QK, RK>
{
    fn ideal_norm(&self, ideal: &IdealsR::Set) -> IdealsZ::Set;

    fn factor_prime_ideal(
        &self,
        prime_ideal: DedekindDomainPrimeIdeal<Z, IdealsZ>,
    ) -> DedekindExtensionIdealFactorsAbovePrime<Z, R, IdealsZ, IdealsR>;

    fn factor_ideal(
        &self,
        ideal: &IdealsR::Set,
    ) -> Option<DedekindExtensionIdealFactorization<Z, R, IdealsZ, IdealsR>>;
}

#[derive(Debug, Clone)]
pub struct DedekindExtensionIdealFactorsAbovePrimeFactor<
    R: DedekindDomainSignature,
    IdealsR: DedekindDomainIdealsSignature<R>,
> {
    pub prime_ideal: DedekindDomainPrimeIdeal<R, IdealsR>,
    pub residue_class_degree: usize,
    pub power: Natural,
}

#[derive(Debug, Clone)]
pub struct DedekindExtensionIdealFactorsAbovePrime<
    Z: DedekindDomainSignature,
    R: DedekindDomainSignature,
    IdealsZ: DedekindDomainIdealsSignature<Z>,
    IdealsR: DedekindDomainIdealsSignature<R>,
> {
    ideals_r: IdealsR,
    base_prime: DedekindDomainPrimeIdeal<Z, IdealsZ>,
    // All factors lie above base_prime
    // All powers are >= 1
    factors: Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<R, IdealsR>>,
}

impl<
    Z: DedekindDomainSignature,
    R: DedekindDomainSignature,
    IdealsZ: DedekindDomainIdealsSignature<Z>,
    IdealsR: DedekindDomainIdealsSignature<R>,
> DedekindExtensionIdealFactorsAbovePrime<Z, R, IdealsZ, IdealsR>
{
    pub fn from_powers_unchecked(
        ideals_r: IdealsR,
        base_prime: DedekindDomainPrimeIdeal<Z, IdealsZ>,
        factors: Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<R, IdealsR>>,
    ) -> Self {
        for f in &factors {
            debug_assert_ne!(f.power, Natural::ZERO);
        }
        Self {
            ideals_r,
            base_prime,
            factors,
        }
    }

    pub fn into_factors(self) -> Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<R, IdealsR>> {
        self.factors
    }

    pub fn into_powers(self) -> Vec<(DedekindDomainPrimeIdeal<R, IdealsR>, Natural)> {
        self.factors
            .into_iter()
            .map(|f| (f.prime_ideal, f.power))
            .collect()
    }

    pub fn unique_prime_factors(&self) -> Vec<&DedekindDomainPrimeIdeal<R, IdealsR>> {
        self.factors.iter().map(|f| &f.prime_ideal).collect()
    }

    pub fn into_full_factorization(self) -> DedekindDomainIdealFactorization<R, IdealsR> {
        DedekindDomainIdealFactorization::from_factor_powers(
            self.ideals_r,
            self.factors
                .into_iter()
                .map(|f| (f.prime_ideal, f.power))
                .collect(),
        )
    }

    /// Do any prime factors appear with multiplicity greater than 1?
    pub fn is_ramified(&self) -> bool {
        self.factors.iter().any(|f| f.power > Natural::ONE)
    }

    /// Is there exactly one prime factor with multiplicity equal to 1?
    pub fn is_inert(&self) -> bool {
        match self.factors.len() {
            1 => {
                let f = &self.factors[0];
                debug_assert_ne!(f.power, Natural::ZERO);
                f.power == Natural::ONE
            }
            _ => false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DedekindExtensionIdealFactorization<
    Z: DedekindDomainSignature,
    R: DedekindDomainSignature,
    IdealsZ: DedekindDomainIdealsSignature<Z>,
    IdealsR: DedekindDomainIdealsSignature<R>,
> {
    ideals_r: IdealsR,
    // Each should be above a different prime
    factors_above_primes: Vec<DedekindExtensionIdealFactorsAbovePrime<Z, R, IdealsZ, IdealsR>>,
}

impl<
    Z: DedekindDomainSignature,
    R: DedekindDomainSignature,
    IdealsZ: DedekindDomainIdealsSignature<Z>,
    IdealsR: DedekindDomainIdealsSignature<R>,
> DedekindExtensionIdealFactorization<Z, R, IdealsZ, IdealsR>
{
    pub fn from_ideal_factors_above_primes(
        ideals_r: IdealsR,
        factors_above_primes: Vec<DedekindExtensionIdealFactorsAbovePrime<Z, R, IdealsZ, IdealsR>>,
    ) -> Self {
        for f in &factors_above_primes {
            debug_assert_eq!(f.ideals_r, ideals_r);
        }
        Self {
            ideals_r,
            factors_above_primes,
        }
    }

    pub fn into_powers(self) -> Vec<(DedekindDomainPrimeIdeal<R, IdealsR>, Natural)> {
        self.into_ring_and_powers().1
    }

    pub fn into_ring_and_powers(self) -> (R, Vec<(DedekindDomainPrimeIdeal<R, IdealsR>, Natural)>) {
        (
            self.ideals_r.ring().clone(),
            self.factors_above_primes
                .into_iter()
                .map(|factors| factors.into_powers())
                .flatten()
                .collect(),
        )
    }

    pub fn into_full_factorization(self) -> DedekindDomainIdealFactorization<R, IdealsR> {
        let ideals_r = self.ideals_r.clone();
        let (_ring, powers) = self.into_ring_and_powers();
        DedekindDomainIdealFactorization::from_factor_powers(ideals_r, powers)
    }
}
