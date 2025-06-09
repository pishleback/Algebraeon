use super::local_ring::{LocalRingIdealsSignature, LocalRingSignature};
use crate::structure::*;
use algebraeon_sets::structure::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct LocalizationPrime<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> {
    pub(crate) ring: RS,
    pub(crate) ring_borrow: PhantomData<RSB>,
    pub(crate) ideal_signature: PS,
    pub(crate) prime_ideal: PS::Set,
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> LocalizationPrime<RS, RSB, PS>
{
    pub fn new_unchecked(ring: RS, ideal_signature: PS, prime_ideal: PS::Set) -> Self {
        assert!(
            !ideal_signature.ideal_equal(&ideal_signature.unit_ideal(), &prime_ideal),
            "The provided prime ideal was the full ring"
        );
        debug_assert!(
            !ideal_signature.ideal_is_zero(&prime_ideal),
            "The provided prime ideal was just 0"
        );
        Self {
            ring,
            ring_borrow: PhantomData,
            ideal_signature,
            prime_ideal,
        }
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> PartialEq for LocalizationPrime<RS, RSB, PS>
{
    fn eq(&self, other: &Self) -> bool {
        self.ring == other.ring
            && self.ideal_signature == other.ideal_signature
            && self
                .ideal_signature
                .ideal_equal(&self.prime_ideal, &other.prime_ideal)
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> Eq for LocalizationPrime<RS, RSB, PS>
{
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> Signature for LocalizationPrime<RS, RSB, PS>
{
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> SetSignature for LocalizationPrime<RS, RSB, PS>
{
    type Set = (RS::Set, RS::Set);

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let (s, _r) = x;
        if self.ring.is_zero(s) {
            return Err(
                "s is zero. It is very much not in the multiplicative set R minus p".to_string(),
            );
        }
        if self
            .ideal_signature
            .ideal_contains_element(&self.prime_ideal, s)
        {
            Err("s is not in the multiplicative set R minus p".to_string())
        } else {
            Ok(())
        }
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> ToStringSignature for LocalizationPrime<RS, RSB, PS>
where
    RS: ToStringSignature,
{
    fn to_string(&self, elem: &Self::Set) -> String {
        let (s, r) = elem;
        format!("{}^(-1) {}", self.ring.to_string(s), self.ring.to_string(r))
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> EqSignature for LocalizationPrime<RS, RSB, PS>
where
    RS: EqSignature,
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa = rb/sb
        let ra_sb = self.ring.mul(ra, sb);
        let rb_sa = self.ring.mul(rb, sa);
        self.ring.equal(&ra_sb, &rb_sa)
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> AdditiveMonoidSignature for LocalizationPrime<RS, RSB, PS>
{
    fn zero(&self) -> Self::Set {
        (self.ring.one(), self.ring.zero())
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa + rb/sb
        let mut num = self.ring.mul(ra, sb);
        let rb_sa = self.ring.mul(rb, sa);
        self.ring.add_mut(&mut num, &rb_sa);
        let den = self.ring.mul(sb, sa);
        (den, num)
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> AdditiveGroupSignature for LocalizationPrime<RS, RSB, PS>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        let (s, r) = a;
        (s.clone(), self.ring.neg(r))
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> SemiRingSignature for LocalizationPrime<RS, RSB, PS>
{
    fn one(&self) -> Self::Set {
        (self.ring.one(), self.ring.one())
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let (sa, ra) = a;
        let (sb, rb) = b;
        // ra/sa * rb/sb
        let num = self.ring.mul(ra, rb);
        let den = self.ring.mul(sb, sa);
        (den, num)
    }
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> RingSignature for LocalizationPrime<RS, RSB, PS>
{
}

impl<
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> SemiRingUnitsSignature for LocalizationPrime<RS, RSB, PS>
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
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
> IntegralDomainSignature for LocalizationPrime<RS, RSB, PS>
{
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut to_return = self.inv(b)?;
        self.mul_mut(&mut to_return, a);
        Ok(to_return)
    }
}

#[derive(Clone, Debug)]
pub struct LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    source: RS,
    target: LocalizationPrime<RS, RSB, PS>,
}

impl<RS, RSB, PS> LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    pub fn new(source: RS, target: LocalizationPrime<RS, RSB, PS>) -> Self {
        Self { source, target }
    }
}

impl<RS, RSB, PS> Morphism<RS, LocalizationPrime<RS, RSB, PS>>
    for LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    fn domain(&self) -> &RS {
        &self.source
    }

    fn range(&self) -> &LocalizationPrime<RS, RSB, PS> {
        &self.target
    }
}

impl<RS, RSB, PS> Function<RS, LocalizationPrime<RS, RSB, PS>>
    for LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    fn image(
        &self,
        x: &<RS as SetSignature>::Set,
    ) -> <LocalizationPrime<RS, RSB, PS> as SetSignature>::Set {
        (self.source.one(), x.clone())
    }
}

impl<RS, RSB, PS> InjectiveFunction<RS, LocalizationPrime<RS, RSB, PS>>
    for LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    fn try_preimage(
        &self,
        y: &<LocalizationPrime<RS, RSB, PS> as SetSignature>::Set,
    ) -> Option<<RS as SetSignature>::Set> {
        let (s, r) = y;
        self.source.div(r, s).ok()
    }
}

impl<RS, RSB, PS> RingHomomorphism<RS, LocalizationPrime<RS, RSB, PS>>
    for LocalizationInclusion<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
}

#[derive(Clone, Debug)]
pub struct LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
    source: LocalizationPrime<RS, RSB, PS>,
    target: FS,
    map_restricted: FSI,
}

impl<RS, RSB, PS, FS, FSI> LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
    pub fn new(source: LocalizationPrime<RS, RSB, PS>, map_restricted: FSI) -> Self {
        let target = map_restricted.range().clone();
        Self {
            source,
            target,
            map_restricted,
        }
    }
}

impl<RS, RSB, PS, FS, FSI> Morphism<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
    fn domain(&self) -> &LocalizationPrime<RS, RSB, PS> {
        &self.source
    }

    fn range(&self) -> &FS {
        &self.target
    }
}

impl<RS, RSB, PS, FS, FSI> Function<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
    fn image(
        &self,
        x: &<LocalizationPrime<RS, RSB, PS> as SetSignature>::Set,
    ) -> <FS as SetSignature>::Set {
        let (s, r) = x;
        let num = self.map_restricted.image(r);
        let den = self.map_restricted.image(s);
        self.target
            .div(&num, &den)
            .expect("s was nonzero in the target field")
    }
}

impl<RS, RSB, PS, FS, FSI> InjectiveFunction<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
    fn try_preimage(
        &self,
        y: &<FS as SetSignature>::Set,
    ) -> Option<<LocalizationPrime<RS, RSB, PS> as SetSignature>::Set> {
        let try_this = self.map_restricted.numerator_and_denominator(y);
        if self.source.is_element(&try_this).is_ok() {
            Some(try_this)
        } else {
            None
        }
    }
}

impl<RS, RSB, PS, FS, FSI> RingHomomorphism<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationInclusionFoF<RS, RSB, PS, FS, FSI>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSI: FieldOfFractionsInclusion<RS, FS>,
{
}

#[derive(Clone, Debug)]
pub struct LocalizationResidueField<RS, RSB, PS, FS, FSQ>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSQ: RingHomomorphism<RS, FS>,
{
    source: LocalizationPrime<RS, RSB, PS>,
    target: FS,
    map_restricted: FSQ,
}

impl<RS, RSB, PS, FS, FSQ> LocalizationResidueField<RS, RSB, PS, FS, FSQ>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSQ: RingHomomorphism<RS, FS>,
{
    pub fn new(source: LocalizationPrime<RS, RSB, PS>, map_restricted: FSQ) -> Self {
        let target = map_restricted.range().clone();
        Self {
            source,
            target,
            map_restricted,
        }
    }
}

impl<RS, RSB, PS, FS, FSQ> Morphism<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationResidueField<RS, RSB, PS, FS, FSQ>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSQ: RingHomomorphism<RS, FS>,
{
    fn domain(&self) -> &LocalizationPrime<RS, RSB, PS> {
        &self.source
    }

    fn range(&self) -> &FS {
        &self.target
    }
}

impl<RS, RSB, PS, FS, FSQ> Function<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationResidueField<RS, RSB, PS, FS, FSQ>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSQ: RingHomomorphism<RS, FS>,
{
    fn image(
        &self,
        x: &<LocalizationPrime<RS, RSB, PS> as SetSignature>::Set,
    ) -> <FS as SetSignature>::Set {
        let (_s, r) = x;
        self.map_restricted.image(r)
    }
}

impl<RS, RSB, PS, FS, FSQ> RingHomomorphism<LocalizationPrime<RS, RSB, PS>, FS>
    for LocalizationResidueField<RS, RSB, PS, FS, FSQ>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    FS: FieldSignature,
    FSQ: RingHomomorphism<RS, FS>,
{
}

#[derive(PartialEq, Eq, Clone, Debug)]
pub struct LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    ambient_ring: LocalizationPrime<RS, RSB, PS>,
}

impl<RS, RSB, PS> Signature for LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
}

impl<RS, RSB, PS> SetSignature for LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    type Set = PS::Set;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        let localized_prime = &self.ambient_ring.prime_ideal;
        let _intersect = self
            .ambient_ring
            .ideal_signature
            .ideal_intersect(x, localized_prime);
        todo!(
            "Do we just consider the extension of x and do not care about intersection with localized_prime?
            not considering if two different I_1,2 of R give the same extension S^-1 I_1,2
            "
        )
    }
}

impl<RS, RSB, PS, LRB> IdealsSignature<LocalizationPrime<RS, RSB, PS>, LRB>
    for LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    LRB: BorrowedStructure<LocalizationPrime<RS, RSB, PS>>,
{
    fn ring(&self) -> &LocalizationPrime<RS, RSB, PS> {
        &self.ambient_ring
    }
}

impl<RS, RSB, PS, LRB> IdealsArithmeticSignature<LocalizationPrime<RS, RSB, PS>, LRB>
    for LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    LRB: BorrowedStructure<LocalizationPrime<RS, RSB, PS>>,
{
    fn principal_ideal(
        &self,
        a: &<LocalizationPrime<RS, RSB, PS> as SetSignature>::Set,
    ) -> Self::Set {
        let (_s, r) = a;
        self.ambient_ring.ideal_signature.principal_ideal(r)
    }

    fn ideal_contains(&self, _a: &Self::Set, _b: &Self::Set) -> bool {
        todo!()
    }

    fn ideal_intersect(&self, _a: &Self::Set, _b: &Self::Set) -> Self::Set {
        todo!()
    }

    fn ideal_add(&self, _a: &Self::Set, _b: &Self::Set) -> Self::Set {
        todo!()
    }

    fn ideal_mul(&self, _a: &Self::Set, _b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl<RS, RSB, PS> LocalRingSignature for LocalizationPrime<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
}

impl<RS, RSB, PS, LRB> LocalRingIdealsSignature<LocalizationPrime<RS, RSB, PS>, LRB>
    for LocalizedIdeals<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
    LRB: BorrowedStructure<LocalizationPrime<RS, RSB, PS>>,
{
    fn unique_maximal(&self) -> <Self as SetSignature>::Set {
        self.ambient_ring.prime_ideal.clone()
    }
}

impl<RS, RSB, PS> RingToIdealsSignature for LocalizationPrime<RS, RSB, PS>
where
    RS: IntegralDomainSignature,
    RSB: BorrowedStructure<RS>,
    PS: IdealsArithmeticSignature<RS, RSB>,
{
    type Ideals<SelfB: BorrowedStructure<Self>> = LocalizedIdeals<RS, RSB, PS>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        LocalizedIdeals {
            ambient_ring: self.clone(),
        }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        LocalizedIdeals { ambient_ring: self }
    }
}
