use crate::structure::{IdealsSignature, RingSignature};
use algebraeon_sets::structure::QuotientSetSignature;
use algebraeon_structures::Signature;
use std::borrow::Cow;
use algebraeon_structures::*;use algebraeon_structures::*;

/// A quotient of a ring by an ideal with its ring structure
pub trait QuotientRingSignature<PreQuoRing: RingSignature>:
    RingSignature + QuotientSetSignature<PreQuoRing>
{
}

/// A quotient ring where we can get a generator of the quotient principal ideal
pub trait QuotientRingGetPrincipalIdealSignature<PreQuoRing: RingSignature>:
    QuotientRingSignature<PreQuoRing>
{
    fn modulus<'a>(&'a self) -> Cow<'a, PreQuoRing::Elem>;
}

/// A quotient ring where we can get the quotient ideal
pub trait QuotientRingGetIdealSignature<
    PreQuoRing: RingSignature,
    PreQuoRingB: BorrowedStructure<PreQuoRing>,
    PreQuoRingIdeal: IdealsSignature<PreQuoRing, PreQuoRingB>,
>: QuotientRingSignature<PreQuoRing>
{
    fn quotient_ideal<'a>(&'a self) -> Cow<'a, PreQuoRingIdeal::Elem>;
}
