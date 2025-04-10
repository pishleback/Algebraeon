use super::*;
use crate::polynomial::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
struct FieldOfFractionsInclusionForIntegralClosure<
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
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
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
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
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
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
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
> Function<R, K> for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn image(&self, x: &R::Set) -> K::Set {
        self.square.r_to_k().image(x)
    }
}

impl<
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
> InjectiveFunction<R, K>
    for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
    fn try_preimage(&self, x: &K::Set) -> Option<R::Set> {
        self.square.r_to_k().try_preimage(x)
    }
}

impl<
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
> RingHomomorphism<R, K>
    for FieldOfFractionsInclusionForIntegralClosure<Z, Q, R, K, ZQ, ZR, QK, RK, ICS>
{
}

impl<
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
    ICS: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>,
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
        debug_assert!(
            self.range()
                .equal(a, &self.range().div(&self.image(&n), &self.image(&d)).unwrap())
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
pub trait IntegralClosureSquare<
    Z: IntegralDomainStructure,
    Q: FieldStructure,
    R: IntegralDomainStructure,
    K: FieldStructure,
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
pub trait FactorableIdealsSquare<
    Z: DedekindDomainStructure,
    Q: FieldStructure,
    R: DedekindDomainStructure,
    K: FieldStructure,
    ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
    QK: FiniteDimensionalFieldExtension<Q, K>,
    ZQ: FieldOfFractionsInclusion<Z, Q>,
    RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
>: IntegralClosureSquare<Z, Q, R, K, ZQ, ZR, QK, RK>
{
    fn factor_prime_ideal(p: &DedekindDomainPrimeIdeal<Z>) -> DedekindDomainIdealFactorization<R>;
}


// #[derive(Clone)]
// pub struct IntegralClosureSquare<
//     Z: IntegralDomainStructure,
//     R: IntegralDomainStructure,
//     Q: FieldStructure,
//     K: FieldStructure,
//     ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
//     QK: FiniteDimensionalFieldExtension<Q, K>,
//     ZQ: FieldOfFractionsInclusion<Z, Q>,
//     RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
// > {
//     z: PhantomData<Z>,
//     r: PhantomData<R>,
//     q: PhantomData<Q>,
//     k: PhantomData<K>,
//     z_to_r: ZR,
//     q_to_k: QK,
//     z_to_q: ZQ,
//     r_to_k: RK,
// }

// impl<
//     Z: IntegralDomainStructure,
//     R: IntegralDomainStructure,
//     Q: FieldStructure,
//     K: FieldStructure,
//     ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
//     QK: FiniteDimensionalFieldExtension<Q, K>,
//     ZQ: FieldOfFractionsInclusion<Z, Q>,
//     RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
// > IntegralClosureSquare<Z, R, Q, K, ZQ, ZR, QK, RK>
// {
//     pub fn new(z_to_r: ZR, q_to_k: QK, z_to_q: ZQ, r_to_k: RK) -> Self {
//         assert_eq!(z_to_r.domain(), z_to_q.domain());
//         assert_eq!(q_to_k.domain(), z_to_q.range());
//         assert_eq!(z_to_r.range(), r_to_k.domain());
//         assert_eq!(q_to_k.range(), r_to_k.range());
//         Self {
//             z: PhantomData,
//             r: PhantomData,
//             q: PhantomData,
//             k: PhantomData,
//             z_to_r,
//             q_to_k,
//             z_to_q,
//             r_to_k,
//         }
//     }
// }
