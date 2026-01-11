use crate::algebraic_number_field::{
    AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature,
    RingOfIntegersToAlgebraicNumberFieldInclusion,
};
use crate::integer::ideal::IntegerIdealsStructure;
use crate::structure::*;
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
pub struct RingOfIntegersIntegralExtension<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
> {
    _r: PhantomData<R>,
    _rb: PhantomData<RB>,
    _k: PhantomData<K>,
    r_to_k: RtoK,
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
> RingOfIntegersIntegralExtension<K, R, RB, RtoK>
{
    pub fn new_integer_extension(r_to_k: RtoK) -> Self {
        Self {
            _r: PhantomData,
            _rb: PhantomData,
            _k: PhantomData,
            r_to_k,
        }
    }

    pub fn with_ideals(
        &self,
    ) -> RingOfIntegersIntegralExtensionWithIdeals<
        K,
        R,
        RB,
        &RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        &R,
        <R as RingToIdealsSignature>::Ideals<&R>,
    >
    where
        R: RingToIdealsSignature,
    {
        let ideals_z = Integer::structure().into_ideals();
        let ideals_r = self.r_to_k.borrow().domain().ideals();
        RingOfIntegersIntegralExtensionWithIdeals::new_integer_extension(
            self.r_to_k.borrow(),
            ideals_z,
            ideals_r,
        )
    }

    pub fn into_with_ideals(
        self,
    ) -> RingOfIntegersIntegralExtensionWithIdeals<
        K,
        R,
        RB,
        RtoK,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        R,
        <R as RingToIdealsSignature>::Ideals<R>,
    >
    where
        R: RingToIdealsSignature,
    {
        let ideals_z = Integer::structure().into_ideals();
        let ideals_r = self.r_to_k().domain().clone().into_ideals();
        RingOfIntegersIntegralExtensionWithIdeals::new_integer_extension(
            self.r_to_k,
            ideals_z,
            ideals_r,
        )
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
> IntegralClosureExtension for RingOfIntegersIntegralExtension<K, R, RB, RtoK>
{
    type QKBasis = K::Basis;
    type Z = IntegerCanonicalStructure;
    type Q = RationalCanonicalStructure;
    type R = R;
    type K = K;
    type ZQ<BZ: BorrowedStructure<Self::Z>, BQ: BorrowedStructure<Self::Q>> =
        PrincipalIntegerMap<RationalCanonicalStructure, RationalCanonicalStructure>;
    type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>> =
        PrincipalIntegerMap<R, BR>;
    type QK<BQ: BorrowedStructure<Self::Q>, BK: BorrowedStructure<Self::K>> =
        K::RationalInclusion<BK>;
    type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>> =
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>;

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
        Cow::Owned(Rational::structure().into_inbound_principal_integer_map())
    }
    fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>> {
        Cow::Owned(self.r_ring().inbound_principal_integer_map())
    }
    fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>> {
        Cow::Owned(
            self.k_field()
                .inbound_finite_dimensional_rational_extension(),
        )
    }
    fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>> {
        Cow::Borrowed(self.r_to_k.borrow())
    }

    fn integralize_multiplier(&self, alpha: &<Self::K as SetSignature>::Set) -> Integer {
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
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
    IdealsZ: IdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RIB: BorrowedStructure<R>,
    IdealsR: IdealsSignature<R, RIB>,
> {
    _r: PhantomData<R>,
    _rb: PhantomData<RB>,
    _k: PhantomData<K>,
    r_to_k: RtoK,
    ideals_z: IdealsZ,
    _rib: PhantomData<RIB>,
    ideals_r: IdealsR,
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
    IdealsZ: IdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RIB: BorrowedStructure<R>,
    IdealsR: IdealsSignature<R, RIB>,
> RingOfIntegersIntegralExtensionWithIdeals<K, R, RB, RtoK, IdealsZ, RIB, IdealsR>
{
    pub fn new_integer_extension(r_to_k: RtoK, ideals_z: IdealsZ, ideals_r: IdealsR) -> Self {
        Self {
            _r: PhantomData,
            _rb: PhantomData,
            _k: PhantomData,
            r_to_k,
            ideals_z,
            _rib: PhantomData,
            ideals_r,
        }
    }

    pub fn z_ideals(&self) -> &IdealsZ {
        &self.ideals_z
    }

    pub fn r_ideals(&self) -> &IdealsR {
        &self.ideals_r
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    RIB: BorrowedStructure<R>,
    IdealsR: DedekindDomainIdealsSignature<R, RIB>,
> IntegralClosureExtension
    for RingOfIntegersIntegralExtensionWithIdeals<K, R, RB, RtoK, IdealsZ, RIB, IdealsR>
{
    type QKBasis = K::Basis;
    type Z = IntegerCanonicalStructure;
    type Q = RationalCanonicalStructure;
    type R = R;
    type K = K;
    type ZQ<BZ: BorrowedStructure<Self::Z>, BQ: BorrowedStructure<Self::Q>> =
        PrincipalIntegerMap<RationalCanonicalStructure, RationalCanonicalStructure>;
    type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>> =
        PrincipalIntegerMap<R, BR>;
    type QK<BQ: BorrowedStructure<Self::Q>, BK: BorrowedStructure<Self::K>> =
        K::RationalInclusion<BK>;
    type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>> =
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>;

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
        Cow::Owned(Rational::structure().into_inbound_principal_integer_map())
    }
    fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>> {
        Cow::Owned(self.r_ring().inbound_principal_integer_map())
    }
    fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>> {
        Cow::Owned(
            self.k_field()
                .inbound_finite_dimensional_rational_extension(),
        )
    }
    fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>> {
        Cow::Borrowed(self.r_to_k.borrow())
    }

    fn integralize_multiplier(&self, alpha: &<Self::K as SetSignature>::Set) -> Integer {
        if self.k_field().is_algebraic_integer(alpha) {
            Integer::ONE
        } else {
            self.k_field().min_poly_denominator_lcm(alpha)
        }
    }
}
