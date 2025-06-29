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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> {
    _ring: PhantomData<Ring>,
    _fof: PhantomData<FoF>,
    ring: RingB,
    fof_inclusion: FoFInclusion,
    ideals: Ideals,
    prime_ideal: Ideals::Set,
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    pub fn new_unchecked(
        ring: RingB,
        ideals: Ideals,
        prime_ideal: Ideals::Set,
        fof_inclusion: FoFInclusion,
    ) -> Self {
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
            _fof: PhantomData,
            ring,
            ideals,
            fof_inclusion,
            prime_ideal,
        }
    }

    pub fn simplify(&self, elt: &<Self as SetSignature>::Set) -> <Self as SetSignature>::Set {
        let fof_inclusion = LocalizationInclusionFoF::new(self.clone());
        let in_field = fof_inclusion.image(elt);
        let back_in_here = fof_inclusion.try_preimage(&in_field);
        if let Some(simplified) = back_in_here {
            simplified
        } else {
            elt.clone()
        }
    }
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> PartialEq for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> Eq for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> Signature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> SetSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> ToStringSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> EqSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> AdditiveMonoidSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> AdditiveGroupSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> SemiRingSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> RingSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
{
}

impl<
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> SemiRingUnitsSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
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
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
> IntegralDomainSignature for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
{
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut to_return = self.inv(b)?;
        self.mul_mut(&mut to_return, a);
        Ok(to_return)
    }
}

/// R includes into S^-1 R
#[derive(Clone, Debug)]
pub struct LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    source: Ring,
    target: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>,
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    pub fn new(
        source: Ring,
        target: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>,
    ) -> Self {
        Self { source, target }
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    Morphism<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>
    for LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn domain(&self) -> &Ring {
        &self.source
    }

    fn range(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> {
        &self.target
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    Function<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>
    for LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn image(
        &self,
        x: &<Ring as SetSignature>::Set,
    ) -> <LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set {
        (self.source.one(), x.clone())
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    InjectiveFunction<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>
    for LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn try_preimage(
        &self,
        y: &<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set,
    ) -> Option<<Ring as SetSignature>::Set> {
        let (s, r) = y;
        self.source.div(r, s).ok()
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    RingHomomorphism<Ring, LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>
    for LocalizationInclusion<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
}

/// S^-1 R has an inclusion map to the field of fractions of R
#[derive(Clone, Debug)]
pub struct LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    source: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FofInclusion>,
    target: FoF,
}

impl<Ring, RingB, Ideals, FoF, FofInclusion>
    LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FofInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FofInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    pub fn new(source: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FofInclusion>) -> Self {
        let target = source.fof_inclusion.range().clone();
        Self { source, target }
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    Morphism<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, FoF>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn domain(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> {
        &self.source
    }

    fn range(&self) -> &FoF {
        &self.target
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    Function<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, FoF>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn image(
        &self,
        x: &<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set,
    ) -> <FoF as SetSignature>::Set {
        let (s, r) = x;
        let num = self.source.fof_inclusion.image(r);
        let den = self.source.fof_inclusion.image(s);
        self.target
            .div(&num, &den)
            .expect("s was nonzero in the target field")
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    InjectiveFunction<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, FoF>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    fn try_preimage(
        &self,
        y: &<FoF as SetSignature>::Set,
    ) -> Option<<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set>
    {
        let mut try_this = self.source.fof_inclusion.numerator_and_denominator(y);
        core::mem::swap(&mut try_this.0, &mut try_this.1);
        if self.source.is_element(&try_this).is_ok() {
            Some(try_this)
        } else {
            None
        }
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion>
    RingHomomorphism<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, FoF>
    for LocalizationInclusionFoF<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
}

/// (R - p)^-1 R has a quotient map to R/p
#[derive(Clone, Debug)]
pub struct LocalizationResidueField<
    Ring,
    RingB,
    Ideals,
    FoF,
    FoFInclusion,
    ResidueField,
    ResidueFieldMap,
> where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
    ResidueField: FieldSignature,
    ResidueFieldMap: RingHomomorphism<Ring, ResidueField>,
{
    source: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>,
    target: ResidueField,
    ring_to_fof: ResidueFieldMap,
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, ResidueField, ResidueFieldMap>
    LocalizationResidueField<Ring, RingB, Ideals, FoF, FoFInclusion, ResidueField, ResidueFieldMap>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
    ResidueField: FieldSignature,
    ResidueFieldMap: RingHomomorphism<Ring, ResidueField>,
{
    pub fn new(
        source: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>,
        ring_to_fof: ResidueFieldMap,
    ) -> Self {
        let target = ring_to_fof.range().clone();
        Self {
            source,
            target,
            ring_to_fof,
        }
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, ResidueField, ResidueFieldMap>
    Morphism<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, ResidueField>
    for LocalizationResidueField<
        Ring,
        RingB,
        Ideals,
        FoF,
        FoFInclusion,
        ResidueField,
        ResidueFieldMap,
    >
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
    ResidueField: FieldSignature,
    ResidueFieldMap: RingHomomorphism<Ring, ResidueField>,
{
    fn domain(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> {
        &self.source
    }

    fn range(&self) -> &ResidueField {
        &self.target
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, ResidueField, ResidueFieldMap>
    Function<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, ResidueField>
    for LocalizationResidueField<
        Ring,
        RingB,
        Ideals,
        FoF,
        FoFInclusion,
        ResidueField,
        ResidueFieldMap,
    >
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
    ResidueField: FieldSignature,
    ResidueFieldMap: RingHomomorphism<Ring, ResidueField>,
{
    fn image(
        &self,
        x: &<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set,
    ) -> <ResidueField as SetSignature>::Set {
        let (_s, r) = x;
        self.ring_to_fof.image(r)
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, ResidueField, ResidueFieldMap>
    RingHomomorphism<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, ResidueField>
    for LocalizationResidueField<
        Ring,
        RingB,
        Ideals,
        FoF,
        FoFInclusion,
        ResidueField,
        ResidueFieldMap,
    >
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
    ResidueField: FieldSignature,
    ResidueFieldMap: RingHomomorphism<Ring, ResidueField>,
{
}

/// This allows storage and manipulations of the ideals of (R-p)^-1 R
#[derive(PartialEq, Eq, Clone, Debug)]
pub struct IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
    ambient_ring: LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>,
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion> Signature
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
{
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion> SetSignature
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
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

impl<Ring, RingB, Ideals, FoF, FoFInclusion, LRB>
    IdealsSignature<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>,
{
    fn ring(&self) -> &LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> {
        &self.ambient_ring
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, LRB>
    IdealsArithmeticSignature<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>,
{
    fn principal_ideal(
        &self,
        a: &<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion> as SetSignature>::Set,
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

impl<Ring, RingB, Ideals, FoF, FoFInclusion> LocalRingSignature
    for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF>,
{
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion, LRB>
    LocalRingIdealsSignature<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>, LRB>
    for IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
    LRB: BorrowedStructure<LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>>,
{
    fn unique_maximal(&self) -> <Self as SetSignature>::Set {
        self.ambient_ring.prime_ideal.clone()
    }
}

impl<Ring, RingB, Ideals, FoF, FoFInclusion> RingToIdealsSignature
    for LocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>
where
    Ring: IntegralDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: IdealsArithmeticSignature<Ring, RingB>,
    FoF: FieldSignature,
    FoFInclusion: FieldOfFractionsInclusion<Ring, FoF> + Eq,
{
    type Ideals<SelfB: BorrowedStructure<Self>> =
        IdealsOfLocalizedRingAtPrime<Ring, RingB, Ideals, FoF, FoFInclusion>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        IdealsOfLocalizedRingAtPrime {
            ambient_ring: self.clone(),
        }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        IdealsOfLocalizedRingAtPrime { ambient_ring: self }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Borrow;

    use algebraeon_nzq::{IntegerCanonicalStructure, Rational, RationalCanonicalStructure};

    use crate::integer::ideal::IntegerIdealsStructure;

    use super::*;

    type ZToQ = PrincipalSubringInclusion<RationalCanonicalStructure, RationalCanonicalStructure>;

    #[allow(type_alias_bounds)]
    type PLocalizeInts<IB: Borrow<IntegerCanonicalStructure>> = LocalizedRingAtPrime<
        IntegerCanonicalStructure,
        IB,
        IntegerIdealsStructure<IB>,
        RationalCanonicalStructure,
        ZToQ,
    >;

    #[test]
    fn two_localize() {
        let integers = IntegerCanonicalStructure {};
        let integer_ideals = integers.ideals();
        let two_ideal = integer_ideals.principal_ideal(&2.into());
        let z_to_q = PrincipalSubringInclusion::<
            RationalCanonicalStructure,
            RationalCanonicalStructure,
        >::new(RationalCanonicalStructure {});
        let two_localize = PLocalizeInts::<&IntegerCanonicalStructure>::new_unchecked(
            &integers,
            integer_ideals,
            two_ideal,
            z_to_q,
        );
        let z_includes =
            LocalizationInclusion::new(IntegerCanonicalStructure {}, two_localize.clone());
        let includes_q = LocalizationInclusionFoF::new(two_localize.clone());

        for no_denominator in -30..30 {
            let in_localization = (1.into(), no_denominator.into());
            let in_localization_string = two_localize.to_string(&in_localization);
            assert_eq!(in_localization_string, format!("1^(-1) {no_denominator}"));
            let in_localization_is_element = two_localize.is_element(&in_localization);
            assert!(in_localization_is_element.is_ok());

            let from_integer = no_denominator.into();
            let in_localization_via_map = z_includes.image(&from_integer);
            assert_eq!(in_localization_via_map, in_localization);
            assert_eq!(
                z_includes.try_preimage(&in_localization),
                Some(from_integer)
            );

            let as_q = includes_q.image(&in_localization);
            assert_eq!(as_q, Rational::from_integers(no_denominator, 1));
            assert_eq!(
                includes_q.try_preimage(&as_q),
                Some((1.into(), no_denominator.into()))
            );
        }

        for two_denominator in -30..30 {
            let in_localization = (2.into(), two_denominator.into());
            let in_localization_string = two_localize.to_string(&in_localization);
            assert_eq!(in_localization_string, format!("2^(-1) {two_denominator}"));
            let in_localization_is_element = two_localize.is_element(&in_localization);
            assert!(in_localization_is_element.is_err());
        }

        for three_denominator in -30..30 {
            let in_localization = (3.into(), three_denominator.into());
            let in_localization_string = two_localize.to_string(&in_localization);
            assert_eq!(
                in_localization_string,
                format!("3^(-1) {three_denominator}")
            );
            let in_localization_is_element: Result<(), String> =
                two_localize.is_element(&in_localization);
            assert!(in_localization_is_element.is_ok());
            let as_simplified = two_localize.simplify(&in_localization);
            assert!(two_localize.is_element(&as_simplified).is_ok());
            assert!(two_localize.equal(&as_simplified, &in_localization));

            let as_q = includes_q.image(&in_localization);
            assert_eq!(as_q, Rational::from_integers(three_denominator, 3));
            let non_reduced_form = (3.into(), three_denominator.into());
            let pre_image = includes_q.try_preimage(&as_q);
            if three_denominator % 3 != 0 {
                assert_eq!(pre_image, Some(non_reduced_form.clone()));
                assert!(two_localize.equal(
                    &pre_image.expect("Already know that it is Some"),
                    &non_reduced_form
                ));
            } else {
                assert_eq!(pre_image, Some((1.into(), (three_denominator / 3).into())));
                assert!(two_localize.equal(
                    &pre_image.expect("Already know that it is Some"),
                    &non_reduced_form
                ));
            }
        }

        let one_half: (_, _) = (2.into(), 1.into());
        let one_half_string = two_localize.to_string(&one_half);
        assert_eq!(one_half_string, "2^(-1) 1");
        let one_half_is_element = two_localize.is_element(&one_half);
        assert!(one_half_is_element.is_err());

        let one_third = (3.into(), 1.into());
        let one_third_string = two_localize.to_string(&one_third);
        assert_eq!(one_third_string, "3^(-1) 1");
        let one_third_is_element = two_localize.is_element(&one_third);
        assert!(one_third_is_element.is_ok());

        let one_third_squared = two_localize.mul(&one_third, &one_third);
        assert_eq!(one_third_squared, (9.into(), 1.into()));

        let mut one_third_squared_three_cubed = (1.into(), 3.into());
        two_localize.mul_mut(&mut one_third_squared_three_cubed, &(1.into(), 3.into()));
        two_localize.mul_mut(&mut one_third_squared_three_cubed, &(1.into(), 3.into()));
        two_localize.mul_mut(&mut one_third_squared_three_cubed, &one_third_squared);
        // This is not in reduced form. So this is where potential problems can occur.
        assert_eq!(one_third_squared_three_cubed, (9.into(), 27.into()));
    }
}
