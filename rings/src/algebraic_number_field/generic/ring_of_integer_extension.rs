use crate::algebraic_number_field::{
    AlgebraicIntegerRingInAlgebraicNumberFieldSignature, AlgebraicIntegerRingSignature,
    AlgebraicNumberFieldSignature, RingOfIntegersToAlgebraicNumberFieldInclusion,
};
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
pub struct RingOfIntegersExtension<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<R, RB>,
> {
    _r: PhantomData<R>,
    _rb: PhantomData<RB>,
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    r_to_k: RtoK,
    ideals_z: IdealsZ,
    ideals_r: IdealsR,
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<R, RB>,
> RingOfIntegersExtension<R, RB, K, KB, RtoK, IdealsZ, IdealsR>
{
    pub fn new_integer_extension(r_to_k: RtoK, ideals_z: IdealsZ, ideals_r: IdealsR) -> Self {
        Self {
            _r: PhantomData,
            _rb: PhantomData,
            _k: PhantomData,
            _kb: PhantomData,
            r_to_k,
            ideals_z,
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
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    RtoK: BorrowedMorphism<R, K, RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>>,
    IdealsZ: DedekindDomainIdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
    IdealsR: DedekindDomainIdealsSignature<R, RB>,
> IntegralClosureExtension for RingOfIntegersExtension<R, RB, K, KB, RtoK, IdealsZ, IdealsR>
where
    RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>:
        AlgebraicIntegerRingInAlgebraicNumberFieldSignature<
                AlgebraicNumberField = K,
                RingOfIntegers = R,
            >,
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
        RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>;

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
