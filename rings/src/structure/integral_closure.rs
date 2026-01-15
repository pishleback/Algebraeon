use super::*;
use crate::polynomial::*;
use crate::valuation::Valuation;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::borrow::Cow;
use std::fmt::Debug;

/// Given a commuting square for an integral closure
///
/// Q → K
/// ↑   ↑
/// Z → R
///
/// Provide an implementation of K as the field of fractions of R
#[derive(Debug, Clone, PartialEq, Eq)]
struct FieldOfFractionsInclusionForIntegralClosure<ICS: IntegralClosureExtension> {
    square: ICS,
}

impl<ICS: IntegralClosureExtension> FieldOfFractionsInclusionForIntegralClosure<ICS> {
    fn new(square: ICS) -> Self {
        Self { square }
    }
}

impl<ICS: IntegralClosureExtension> Morphism<ICS::R, ICS::K>
    for FieldOfFractionsInclusionForIntegralClosure<ICS>
{
    fn domain(&self) -> &ICS::R {
        self.square.r_ring()
    }

    fn range(&self) -> &ICS::K {
        self.square.k_field()
    }
}

impl<ICS: IntegralClosureExtension> Function<ICS::R, ICS::K>
    for FieldOfFractionsInclusionForIntegralClosure<ICS>
{
    fn image(&self, x: &<ICS::R as SetSignature>::Set) -> <ICS::K as SetSignature>::Set {
        self.square.r_to_k().image(x)
    }
}

impl<ICS: IntegralClosureExtension> InjectiveFunction<ICS::R, ICS::K>
    for FieldOfFractionsInclusionForIntegralClosure<ICS>
{
    fn try_preimage(
        &self,
        x: &<ICS::K as SetSignature>::Set,
    ) -> Option<<ICS::R as SetSignature>::Set> {
        self.square.r_to_k().try_preimage(x)
    }
}

impl<ICS: IntegralClosureExtension> RingHomomorphism<ICS::R, ICS::K>
    for FieldOfFractionsInclusionForIntegralClosure<ICS>
{
}

impl<ICS: IntegralClosureExtension> FieldOfFractionsInclusion<ICS::R, ICS::K>
    for FieldOfFractionsInclusionForIntegralClosure<ICS>
{
    fn numerator_and_denominator(
        &self,
        a: &<ICS::K as SetSignature>::Set,
    ) -> (<ICS::R as SetSignature>::Set, <ICS::R as SetSignature>::Set) {
        // let d in Z such that d*a is in R
        let d = self.square.integralize_multiplier(a);
        // take d in R
        let d = self.square.z_to_r().image(&d);
        // now a = (d*a) / d
        let n = self
            .try_preimage(&self.range().mul(&self.image(&d), a))
            .unwrap();
        debug_assert!(
            self.range().equal(
                a,
                &self
                    .range()
                    .try_divide(&self.image(&n), &self.image(&d))
                    .unwrap()
            )
        );
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
pub trait IntegralClosureExtension: Debug + Clone + Send + Sync {
    type QKBasis: FiniteSetSignature;
    type Z: IntegralDomainSignature;
    type Q: FieldSignature;
    type R: IntegralDomainSignature;
    type K: FieldSignature;
    type ZQ<BZ : BorrowedStructure<Self::Z>, BQ : BorrowedStructure<Self::Q>>: FieldOfFractionsInclusion<Self::Z, Self::Q>;
    type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>>: RingHomomorphism<Self::Z, Self::R>
        + InjectiveFunction<Self::Z, Self::R>;
    type QK<BQ : BorrowedStructure<Self::Q>, BK : BorrowedStructure<Self::K>>: FiniteDimensionalFieldExtension<Self::Q, Self::K>;
    type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>>: RingHomomorphism<Self::R, Self::K>
        + InjectiveFunction<Self::R, Self::K>;

    fn z_ring(&self) -> &Self::Z;
    fn q_field(&self) -> &Self::Q;
    fn r_ring(&self) -> &Self::R;
    fn k_field(&self) -> &Self::K;

    fn z_to_q<'a>(&'a self) -> Cow<'a, Self::ZQ<&'a Self::Z, &'a Self::Q>>;
    fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>>;
    fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>>;
    fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>>;

    /// The square should commute, so this should be both
    /// - `z_to_q` followed by `q_to_k`
    /// - `z_to_r` followed by `r_to_k`
    fn z_to_k(
        &self,
    ) -> impl RingHomomorphism<Self::Z, Self::K> + InjectiveFunction<Self::Z, Self::K> {
        CompositionMorphism::new(self.z_to_q().into_owned(), self.q_to_k().into_owned())
    }

    /// The monic minimal polynomial of alpha in K over Q
    fn min_poly_k_over_q(
        &self,
        alpha: &<Self::K as SetSignature>::Set,
    ) -> Polynomial<<Self::Q as SetSignature>::Set> {
        let alpha_min_poly_monic = self.q_to_k().min_poly(alpha);
        #[cfg(debug_assertions)]
        {
            let q_field = self.q_field();
            let q_poly = q_field.polynomials();
            assert!(q_poly.is_monic(&alpha_min_poly_monic));
        }
        alpha_min_poly_monic
    }

    /// By definition of R as the integral closure of Z in K every element of R, when considered as an element of K, has minimal polynomial over Q which is monic with coefficients in Z
    fn min_poly_r_over_z(
        &self,
        alpha: &<Self::R as SetSignature>::Set,
    ) -> Polynomial<<Self::Z as SetSignature>::Set> {
        let alpha_min_poly_monic = self
            .min_poly_k_over_q(&self.r_to_k().image(alpha))
            .apply_map_into(|c| self.z_to_q().try_preimage(&c).unwrap());
        #[cfg(debug_assertions)]
        {
            let z_ring = self.z_ring();
            let z_poly = z_ring.polynomials();
            assert!(z_poly.is_monic(&alpha_min_poly_monic));
        }
        alpha_min_poly_monic
    }

    /// For alpha in K return non-zero d in Z such that d*alpha is in R
    fn integralize_multiplier(
        &self,
        alpha: &<Self::K as SetSignature>::Set,
    ) -> <Self::Z as SetSignature>::Set;

    fn integral_scalar_multiple_r(
        &self,
        alpha: &<Self::K as SetSignature>::Set,
    ) -> <Self::R as SetSignature>::Set {
        self.r_to_k()
            .try_preimage(&self.integral_scalar_multiple_k(alpha))
            .unwrap()
    }

    fn integral_scalar_multiple_k(
        &self,
        alpha: &<Self::K as SetSignature>::Set,
    ) -> <Self::K as SetSignature>::Set {
        let d = self.integralize_multiplier(alpha);
        // This is the LCM of the denominators of the coefficients of a,
        // and thus it may well be > 1 even when the element is an algebraic integer.
        self.k_field()
            .mul(&self.r_to_k().image(&self.z_to_r().image(&d)), alpha)
    }

    /// Every element of K is a fraction of elements of R
    fn r_to_k_field_of_fractions(&self) -> impl FieldOfFractionsInclusion<Self::R, Self::K> {
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
/// This trait allows the ideal pR of R to be factored into prime ideals in R for each prime ideal p of Z
pub trait DedekindDomainExtension<ZB: BorrowedStructure<Self::Z>, RB: BorrowedStructure<Self::R>>:
    IntegralClosureExtension
where
    Self::Z: DedekindDomainSignature,
    Self::R: DedekindDomainSignature,
{
    type IdealsZ: DedekindDomainIdealsSignature<Self::Z, ZB>;
    type IdealsR: DedekindDomainIdealsSignature<Self::R, RB>;

    fn z_ideals(&self) -> &Self::IdealsZ;

    fn r_ideals(&self) -> &Self::IdealsR;

    fn ideal_norm(
        &self,
        ideal: &<Self::IdealsR as SetSignature>::Set,
    ) -> <Self::IdealsZ as SetSignature>::Set;

    fn factor_prime_ideal(
        &self,
        prime_ideal: <Self::IdealsZ as SetSignature>::Set,
    ) -> DedekindExtensionIdealFactorsAbovePrime<
        <Self::IdealsZ as SetSignature>::Set,
        <Self::IdealsR as SetSignature>::Set,
    >;

    fn factor_ideal(
        &self,
        ideal: &<Self::IdealsR as SetSignature>::Set,
    ) -> Option<
        DedekindExtensionIdealFactorization<
            <Self::IdealsZ as SetSignature>::Set,
            <Self::IdealsR as SetSignature>::Set,
        >,
    >;

    fn padic_k_element_valuation(
        &self,
        prime: &<Self::IdealsR as SetSignature>::Set,
        a: &<Self::K as SetSignature>::Set,
    ) -> Valuation {
        #[cfg(debug_assertions)]
        assert_ne!(self.r_ideals().try_is_irreducible(prime), Some(false));
        let d = self.integralize_multiplier(a);
        let m = self.k_field().mul(a, &self.z_to_k().image(&d));
        self.r_ideals()
            .padic_r_element_valuation(prime, &self.r_to_k().try_preimage(&m).unwrap())
            - self
                .r_ideals()
                .padic_r_element_valuation(prime, &self.z_to_r().image(&d))
    }

    // An element is S-integral, if its valuations at all primes not in S are nonnegative.
    // If S is the empty set, this coincides with the usual integrality.
    #[allow(non_snake_case)]
    fn is_S_integral(
        &self,
        S: Vec<&<Self::IdealsR as SetSignature>::Set>,
        a: &<Self::K as SetSignature>::Set,
    ) -> bool {
        #[cfg(debug_assertions)]
        for s in &S {
            assert_ne!(self.r_ideals().try_is_irreducible(s), Some(false));
        }

        let d = self.integralize_multiplier(a);
        let m = self.k_field().mul(a, &self.z_to_k().image(&d));
        // for each prime factor P of d not in S, check if valuation_P(m) ≥ valuation_P(d)

        let d_as_roi = self.z_to_r().image(&d);
        let principal_ideal_d = self.r_ideals().generated_ideal(vec![d_as_roi.clone()]);

        let d_factorization = self.factor_ideal(&principal_ideal_d);
        if d_factorization.is_none() {
            return true;
        }

        for (prime, _) in d_factorization.unwrap().into_powers() {
            // Skip primes in S
            if S.iter()
                .any(|s_ideal| self.r_ideals().equal(s_ideal, &prime))
            {
                continue;
            }

            let m_val = self
                .r_ideals()
                .padic_r_element_valuation(&prime, &self.r_to_k().try_preimage(&m).unwrap());
            let d_val = self.r_ideals().padic_r_element_valuation(&prime, &d_as_roi);

            if m_val < d_val {
                return false;
            }
        }

        true
    }

    fn expand_extension_ideal_factorization(
        &self,
        f: &DedekindExtensionIdealFactorization<
            <Self::IdealsZ as SetSignature>::Set,
            <Self::IdealsR as SetSignature>::Set,
        >,
    ) -> Factored<<Self::IdealsR as SetSignature>::Set, Natural> {
        self.r_ideals()
            .factorizations()
            .new_unit_and_powers_unchecked(self.r_ideals().one(), f.clone().into_powers())
    }

    /*
      pub fn into_full_factorization(self) -> IdealR {
        DedekindDomainIdealFactorization::from_factor_powers(
            self.factors
                .into_iter()
                .map(|f| (f.prime_ideal, f.power))
                .collect(),
        )
    }

     pub fn into_full_factorization(self) -> Factored<<IdealR as SetSignature>::Set, Natural> {

        Factored::NonZero(NonZeroFactored { unit: , powers: () })

        Factored:: ::from_factor_powers(self.into_powers())
    }
    */
}

#[derive(Debug, Clone)]
pub struct DedekindExtensionIdealFactorsAbovePrimeFactor<Ideal> {
    pub prime_ideal: Ideal,
    pub residue_class_degree: usize,
    pub power: Natural,
}

#[derive(Debug, Clone)]
pub struct DedekindExtensionIdealFactorsAbovePrime<IdealZ, IdealR> {
    #[allow(unused)]
    base_prime: IdealZ,
    // All factors lie above base_prime
    // All powers are >= 1
    factors: Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<IdealR>>,
}

impl<IdealZ, IdealR> DedekindExtensionIdealFactorsAbovePrime<IdealZ, IdealR> {
    pub fn from_powers_unchecked(
        base_prime: IdealZ,
        factors: Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<IdealR>>,
    ) -> Self {
        for f in &factors {
            debug_assert_ne!(f.power, Natural::ZERO);
        }
        Self {
            base_prime,
            factors,
        }
    }

    pub fn into_factors(self) -> Vec<DedekindExtensionIdealFactorsAbovePrimeFactor<IdealR>> {
        self.factors
    }

    pub fn into_powers(self) -> Vec<(IdealR, Natural)> {
        self.factors
            .into_iter()
            .map(|f| (f.prime_ideal, f.power))
            .collect()
    }

    pub fn unique_prime_factors(&self) -> Vec<&IdealR> {
        self.factors.iter().map(|f| &f.prime_ideal).collect()
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
pub struct DedekindExtensionIdealFactorization<IdealZ, IdealR> {
    // Each should be above a different prime
    factors_above_primes: Vec<DedekindExtensionIdealFactorsAbovePrime<IdealZ, IdealR>>,
}

impl<IdealZ, IdealR> DedekindExtensionIdealFactorization<IdealZ, IdealR> {
    pub fn from_ideal_factors_above_primes(
        factors_above_primes: Vec<DedekindExtensionIdealFactorsAbovePrime<IdealZ, IdealR>>,
    ) -> Self {
        Self {
            factors_above_primes,
        }
    }

    pub fn into_powers(self) -> Vec<(IdealR, Natural)> {
        #[allow(clippy::redundant_closure_for_method_calls)]
        self.factors_above_primes
            .into_iter()
            .flat_map(|factors| factors.into_powers())
            .collect()
    }
}
