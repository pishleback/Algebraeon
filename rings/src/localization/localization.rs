use super::local_ring::{LocalRingIdealsSignature, LocalRingSignature};
use crate::structure::*;
use algebraeon_sets::structure::*;
use std::marker::PhantomData;

/// The ring obtained by localizing at the
/// multiplicative set (R - p)
#[derive(Debug, Clone)]
pub struct LocalizedRingAtPrime<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> {
    _ring: PhantomData<Ring>,
    ring: RingB,
    ideals: Ideals,
    prime_ideal: Ideals::Set,
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    pub fn ideals(&self) -> &Ideals {
        &self.ideals
    }

    /// Assuming that `prime_ideal` was actually a prime ideal
    /// for `ring`, then create an instance of this struct
    /// in order to implement `RingSignature` for the localized ring.
    /// This does not check that `prime_ideal` was actually prime,
    /// but we do the sanity checks that it is a proper ideal.
    pub fn new_unchecked(ring: RingB, ideals: Ideals, prime_ideal: Ideals::Set) -> Self {
        assert!(
            !ideals.ideal_equal(&ideals.unit_ideal(), &prime_ideal),
            "The provided prime ideal was the full ring"
        );
        debug_assert!(
            !ideals.ideal_is_zero(&prime_ideal),
            "The provided prime ideal was just 0"
        );
        Self {
            _ring: PhantomData,
            ring,
            ideals,
            prime_ideal,
        }
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> PartialEq for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn eq(&self, other: &Self) -> bool {
        self.ring == other.ring
            && self.ideals == other.ideals
            && self
                .ideals
                .ideal_equal(&self.prime_ideal, &other.prime_ideal)
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> Eq for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> Signature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> SetSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    // (s, r) represents the fraction r/s
    type Set = (Ring::Set, Ring::Set);

    // r can be anything in the ring while s must not belong to the prime ideal
    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let (s, _r) = x;
        if self.ring().is_zero(s) {
            return Err("s is zero. It is not in the multiplicative set R\\p".to_string());
        }
        if self.ideals.ideal_contains_element(&self.prime_ideal, s) {
            Err("s is not in the multiplicative set R\\p".to_string())
        } else {
            Ok(())
        }
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> ToStringSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: ToStringSignature,
{
    fn to_string(&self, elem: &Self::Set) -> String {
        let (s, r) = elem;
        format!(
            "{}^(-1) {}",
            self.ring().to_string(s),
            self.ring().to_string(r)
        )
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> EqSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: EqSignature,
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa = rb/sb
        let ra_sb = self.ring().mul(ra, sb);
        let rb_sa = self.ring().mul(rb, sa);
        self.ring().equal(&ra_sb, &rb_sa)
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> AdditiveMonoidSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn zero(&self) -> Self::Set {
        (self.ring().one(), self.ring().zero())
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa + rb/sb
        let mut num = self.ring().mul(ra, sb);
        let rb_sa = self.ring().mul(rb, sa);
        self.ring().add_mut(&mut num, &rb_sa);
        let den = self.ring().mul(sb, sa);
        (den, num)
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> AdditiveGroupSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        let (s, r) = a;
        (s.clone(), self.ring().neg(r))
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> SemiRingSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn one(&self) -> Self::Set {
        (self.ring().one(), self.ring().one())
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa * rb/sb
        let num = self.ring().mul(ra, rb);
        let den = self.ring().mul(sb, sa);
        (den, num)
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> RingSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> SemiRingUnitsSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        if self.is_zero(a) {
            Err(RingDivisionError::DivideByZero)
        } else {
            let rs = (a.1.clone(), a.0.clone());
            if self.is_element(&rs).is_ok() {
                Ok(rs)
            } else {
                Err(RingDivisionError::NotDivisible)
            }
        }
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
> IntegralDomainSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
{
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut to_return = self.inv(b)?;
        self.mul_mut(&mut to_return, a);
        Ok(to_return)
    }
}

/// R includes into S^-1 R
#[derive(Clone, Debug)]
pub struct LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    source: Ring,
    target: LocalizedRingAtPrime<Ring, RingB, Ideals>,
}

impl<Ring, RingB, Ideals> LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    pub fn new(source: Ring, target: LocalizedRingAtPrime<Ring, RingB, Ideals>) -> Self {
        Self { source, target }
    }
}

impl<Ring, RingB, Ideals> Morphism<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals>>
    for LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    fn domain(&self) -> &Ring {
        &self.source
    }

    fn range(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals> {
        &self.target
    }
}

impl<Ring, RingB, Ideals> Function<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals>>
    for LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    fn image(
        &self,
        x: &<Ring as SetSignature>::Set,
    ) -> <LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set {
        (self.source.one(), x.clone())
    }
}

impl<Ring, RingB, Ideals> InjectiveFunction<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals>>
    for LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    fn try_preimage(
        &self,
        y: &<LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set,
    ) -> Option<<Ring as SetSignature>::Set> {
        let (s, r) = y;
        self.source.div(r, s).ok()
    }
}

impl<Ring, RingB, Ideals> RingHomomorphism<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals>>
    for LocalizationInclusion<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
}

/// S^-1 R has an inclusion map to the field of fractions of R
#[derive(Clone, Debug)]
pub struct LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
    source: LocalizedRingAtPrime<Ring, RingB, Ideals>,
    target: Fof,
    // this is just R into its field of fractions
    ring_to_fof: FofInclusion,
}

impl<Ring, RingB, Ideals, Fof, FofInclusion>
    LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
    pub fn new(
        source: LocalizedRingAtPrime<Ring, RingB, Ideals>,
        ring_to_fof: FofInclusion,
    ) -> Self {
        let target = ring_to_fof.range().clone();
        Self {
            source,
            target,
            ring_to_fof,
        }
    }
}

impl<Ring, RingB, Ideals, Fof, FofInclusion>
    Morphism<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
    fn domain(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals> {
        &self.source
    }

    fn range(&self) -> &Fof {
        &self.target
    }
}

impl<Ring, RingB, Ideals, Fof, FofInclusion>
    Function<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
    fn image(
        &self,
        x: &<LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set,
    ) -> <Fof as SetSignature>::Set {
        let (s, r) = x;
        let num = self.ring_to_fof.image(r);
        let den = self.ring_to_fof.image(s);
        self.target
            .div(&num, &den)
            .expect("s was nonzero in the target field")
    }
}

impl<Ring, RingB, Ideals, Fof, FofInclusion>
    InjectiveFunction<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
    fn try_preimage(
        &self,
        y: &<Fof as SetSignature>::Set,
    ) -> Option<<LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set> {
        let try_this = self.ring_to_fof.numerator_and_denominator(y);
        if self.source.is_element(&try_this).is_ok() {
            Some(try_this)
        } else {
            None
        }
    }
}

impl<Ring, RingB, Ideals, Fof, FofInclusion>
    RingHomomorphism<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, Fof, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, Fof>,
{
}

/// (R - p)^-1 R has a quotient map to R/p
#[derive(Clone, Debug)]
pub struct LocalizationResidueField<Ring, RingB, Ideals, Fof, FofQ>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofQ: RingHomomorphism<Ring, Fof>,
{
    source: LocalizedRingAtPrime<Ring, RingB, Ideals>,
    target: Fof,
    ring_to_fof: FofQ,
}

impl<Ring, RingB, Ideals, Fof, FofQ> LocalizationResidueField<Ring, RingB, Ideals, Fof, FofQ>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofQ: RingHomomorphism<Ring, Fof>,
{
    pub fn new(source: LocalizedRingAtPrime<Ring, RingB, Ideals>, ring_to_fof: FofQ) -> Self {
        let target = ring_to_fof.range().clone();
        Self {
            source,
            target,
            ring_to_fof,
        }
    }
}

impl<Ring, RingB, Ideals, Fof, FofQ> Morphism<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationResidueField<Ring, RingB, Ideals, Fof, FofQ>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofQ: RingHomomorphism<Ring, Fof>,
{
    fn domain(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals> {
        &self.source
    }

    fn range(&self) -> &Fof {
        &self.target
    }
}

impl<Ring, RingB, Ideals, Fof, FofQ> Function<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationResidueField<Ring, RingB, Ideals, Fof, FofQ>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofQ: RingHomomorphism<Ring, Fof>,
{
    fn image(
        &self,
        x: &<LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set,
    ) -> <Fof as SetSignature>::Set {
        let (_s, r) = x;
        self.ring_to_fof.image(r)
    }
}

impl<Ring, RingB, Ideals, Fof, FofQ>
    RingHomomorphism<LocalizedRingAtPrime<Ring, RingB, Ideals>, Fof>
    for LocalizationResidueField<Ring, RingB, Ideals, Fof, FofQ>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    Fof: FieldSignature,
    FofQ: RingHomomorphism<Ring, Fof>,
{
}

/// This allows storage and manipulations of the ideals of (R-p)^-1 R
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    ambient_ring: LocalizedRingAtPrime<Ring, RingB, Ideals>,
}

impl<Ring, RingB, Ideals> Signature for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
}

impl<Ring, RingB, Ideals> SetSignature for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    type Set = Ideals::Set;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let localized_prime = &self.ambient_ring.prime_ideal;
        let _intersect = self.ambient_ring.ideals.ideal_intersect(x, localized_prime);
        todo!(
            "Do we just consider the extension of x and do not care about intersection with localized_prime?
            not considering if two different I_1,2 of R give the same extension S^-1 I_1,2
            "
        )
    }
}

impl<Ring, RingB, Ideals, LRB> IdealsSignature<LocalizedRingAtPrime<Ring, RingB, Ideals>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals>>,
{
    fn ring(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals> {
        &self.ambient_ring
    }
}

impl<Ring, RingB, Ideals, LRB>
    IdealsArithmeticSignature<LocalizedRingAtPrime<Ring, RingB, Ideals>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: DedekindDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: FactorableIdealsSignature<Ring, RingB>,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals>>,
{
    fn principal_ideal(
        &self,
        a: &<LocalizedRingAtPrime<Ring, RingB, Ideals> as SetSignature>::Set,
    ) -> Self::Set {
        let (_s, r) = a;
        self.ambient_ring.ideals.principal_ideal(r)
    }

    fn ideal_contains(&self, a: &Self::Set, b: &Self::Set) -> bool {
        // S^{-1} a >= S^{-1} b
        if self.ambient_ring.ideals.ideal_contains(a, b) {
            true
        } else if self.ambient_ring.ideals.ideal_is_zero(a) {
            false
        } else {
            todo!()
        }
    }

    fn ideal_intersect(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ambient_ring.ideals.ideal_intersect(a, b)
    }

    fn ideal_add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ambient_ring.ideals.ideal_add(a, b)
    }

    fn ideal_mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.ambient_ring.ideals.ideal_mul(a, b)
    }
}

impl<Ring, RingB, Ideals> LocalRingSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
}

impl<Ring, RingB, Ideals, LRB>
    LocalRingIdealsSignature<LocalizedRingAtPrime<Ring, RingB, Ideals>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals>>,
{
    fn unique_maximal(&self) -> <Self as SetSignature>::Set {
        self.ambient_ring.prime_ideal.clone()
    }
}

impl<Ring, RingB, Ideals> RingToIdealsSignature for LocalizedRingAtPrime<Ring, RingB, Ideals>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
{
    type Ideals<SelfB: BorrowedStructure<Self>> = IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        IdealsOfLocalizedRingAtPrime {
            ambient_ring: self.clone(),
        }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        IdealsOfLocalizedRingAtPrime { ambient_ring: self }
    }
}
