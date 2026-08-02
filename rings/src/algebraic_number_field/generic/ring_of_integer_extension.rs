use crate::algebraic_number_field::{
    AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature,
    RingOfIntegersToAlgebraicNumberFieldInclusion,
};
use crate::num_theory::integer_ideal::IntegerIdealsStructure;
use crate::structure::*;
use algebraeon_structures::*;
use std::sync::Arc;

/// Q -> K
/// ↑    ↑
/// Z -> R
///
/// Where Q is the rationals, Z is the integers, K is an algebraic number field, R is its ring of integers
///
#[derive(Debug, Clone)]
pub struct RingOfIntegersIntegralExtension<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
> {
    r_to_k: Arc<RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>>,
}

impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
    RingOfIntegersIntegralExtension<K, R>
{
    pub fn new_integer_extension(
        r_to_k: Arc<RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>>,
    ) -> Arc<Self> {
        Self { r_to_k }.into()
    }

    pub fn with_ideals(
        self: &Arc<Self>,
    ) -> Arc<
        RingOfIntegersIntegralExtensionWithIdeals<
            K,
            R,
            IntegerIdealsStructure,
            <R as RingToIdealsSignature>::Ideals,
        >,
    >
    where
        R: RingToIdealsSignature,
    {
        let ideals_z = Integer::structure().ideals();
        let ideals_r = self.r_to_k.domain().ideals();
        RingOfIntegersIntegralExtensionWithIdeals::new_integer_extension(
            self.r_to_k.clone(),
            ideals_z,
            ideals_r,
        )
    }
}

impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>> IntegralClosureExtension
    for RingOfIntegersIntegralExtension<K, R>
{
    type QKBasis = K::Basis;
    type Z = IntegerCanonicalStructure;
    type Q = RationalCanonicalStructure;
    type R = R;
    type K = K;
    type ZQ = PrincipalIntegerMap<RationalCanonicalStructure>;
    type ZR = PrincipalIntegerMap<R>;
    type QK = K::RationalInclusion;
    type RK = RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>;

    fn z_ring(self: &Arc<Self>) -> Arc<Self::Z> {
        Integer::structure()
    }
    fn r_ring(self: &Arc<Self>) -> Arc<Self::R> {
        self.r_to_k.domain()
    }
    fn q_field(self: &Arc<Self>) -> Arc<Self::Q> {
        Rational::structure()
    }
    fn k_field(self: &Arc<Self>) -> Arc<Self::K> {
        self.r_to_k.range()
    }

    fn z_to_q(self: &Arc<Self>) -> Arc<Self::ZQ> {
        Rational::structure().inbound_principal_integer_map()
    }
    fn z_to_r(self: &Arc<Self>) -> Arc<Self::ZR> {
        self.r_ring().inbound_principal_integer_map()
    }
    fn q_to_k(self: &Arc<Self>) -> Arc<Self::QK> {
        self.k_field()
            .inbound_finite_dimensional_rational_extension()
    }
    fn r_to_k(self: &Arc<Self>) -> Arc<Self::RK> {
        self.r_to_k.clone()
    }

    fn integralize_multiplier(
        self: &Arc<Self>,
        alpha: &<Self::K as SetSignature>::Elem,
    ) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            self.k_field().min_poly_denominator_lcm(alpha)
        }
    }
}

/// Q -> K
/// ↑    ↑
/// Z -> R
///
/// Where Q is the rationals, Z is the integers, K is an algebraic number field, R is its ring of integers
///
#[derive(Debug, Clone)]
pub struct RingOfIntegersIntegralExtensionWithIdeals<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    IdealsZ: IdealsSignature<IntegerCanonicalStructure>,
    IdealsR: IdealsSignature<R>,
> {
    r_to_k: Arc<RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>>,
    ideals_z: Arc<IdealsZ>,
    ideals_r: Arc<IdealsR>,
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    IdealsZ: IdealsSignature<IntegerCanonicalStructure>,
    IdealsR: IdealsSignature<R>,
> RingOfIntegersIntegralExtensionWithIdeals<K, R, IdealsZ, IdealsR>
{
    pub fn new_integer_extension(
        r_to_k: Arc<RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>>,
        ideals_z: Arc<IdealsZ>,
        ideals_r: Arc<IdealsR>,
    ) -> Arc<Self> {
        Self {
            r_to_k,
            ideals_z,
            ideals_r,
        }
        .into()
    }

    pub fn z_ideals(self: &Arc<Self>) -> &Arc<IdealsZ> {
        &self.ideals_z
    }

    pub fn r_ideals(self: &Arc<Self>) -> &Arc<IdealsR> {
        &self.ideals_r
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<R>,
> IntegralClosureExtension for RingOfIntegersIntegralExtensionWithIdeals<K, R, IdealsZ, IdealsR>
{
    type QKBasis = K::Basis;
    type Z = IntegerCanonicalStructure;
    type Q = RationalCanonicalStructure;
    type R = R;
    type K = K;
    type ZQ = PrincipalIntegerMap<RationalCanonicalStructure>;
    type ZR = PrincipalIntegerMap<R>;
    type QK = K::RationalInclusion;
    type RK = RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>;

    fn z_ring(self: &Arc<Self>) -> Arc<Self::Z> {
        Integer::structure()
    }
    fn r_ring(self: &Arc<Self>) -> Arc<Self::R> {
        self.r_to_k.domain()
    }
    fn q_field(self: &Arc<Self>) -> Arc<Self::Q> {
        Rational::structure()
    }
    fn k_field(self: &Arc<Self>) -> Arc<Self::K> {
        self.r_to_k.range()
    }

    fn z_to_q(self: &Arc<Self>) -> Arc<Self::ZQ> {
        Rational::structure().inbound_principal_integer_map()
    }
    fn z_to_r(self: &Arc<Self>) -> Arc<Self::ZR> {
        self.r_ring().inbound_principal_integer_map()
    }
    fn q_to_k(self: &Arc<Self>) -> Arc<Self::QK> {
        self.k_field()
            .inbound_finite_dimensional_rational_extension()
    }
    fn r_to_k(self: &Arc<Self>) -> Arc<Self::RK> {
        self.r_to_k.clone()
    }

    fn integralize_multiplier(
        self: &Arc<Self>,
        alpha: &<Self::K as SetSignature>::Elem,
    ) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            self.k_field().min_poly_denominator_lcm(alpha)
        }
    }
}
