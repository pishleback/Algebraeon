use itertools::Itertools;

use super::*;
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait Morphism<Domain: Signature, Range: Signature>: Debug + Clone {
    fn domain(&self) -> &Domain;
    fn range(&self) -> &Range;
}

pub trait Function<Domain: SetSignature, Range: SetSignature>: Morphism<Domain, Range> {
    fn image(&self, x: &Domain::Set) -> Range::Set;
}

pub trait InjectiveFunction<Domain: SetSignature, Range: SetSignature>:
    Function<Domain, Range>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<Domain::Set>;
}

pub trait BijectiveFunction<Domain: SetSignature, Range: SetSignature>:
    InjectiveFunction<Domain, Range>
{
    fn preimage(&self, x: &Range::Set) -> Domain::Set {
        self.try_preimage(x).unwrap()
    }
}

/// The identity morphism X -> X
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityMorphism<X: Signature> {
    x: X,
}

impl<X: Signature> IdentityMorphism<X> {
    pub fn new(x: X) -> Self {
        Self { x }
    }
}

impl<X: Signature> Morphism<X, X> for IdentityMorphism<X> {
    fn domain(&self) -> &X {
        &self.x
    }

    fn range(&self) -> &X {
        &self.x
    }
}

impl<X: SetSignature> Function<X, X> for IdentityMorphism<X> {
    fn image(&self, x: &X::Set) -> X::Set {
        x.clone()
    }
}

impl<X: SetSignature> InjectiveFunction<X, X> for IdentityMorphism<X> {
    fn try_preimage(&self, x: &X::Set) -> Option<X::Set> {
        Some(x.clone())
    }
}

impl<X: SetSignature> BijectiveFunction<X, X> for IdentityMorphism<X> {}

/// The composition A -> B -> C of two morphisms A -> B and B -> C
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionMorphism<
    A: Signature,
    B: Signature,
    C: Signature,
    AB: Morphism<A, B>,
    BC: Morphism<B, C>,
> {
    a: PhantomData<A>,
    b: PhantomData<B>,
    c: PhantomData<C>,
    // required to satisfy ab.range() == bc.domain()
    a_to_b: AB,
    b_to_c: BC,
}

impl<A: Signature, B: Signature, C: Signature, AB: Morphism<A, B>, BC: Morphism<B, C>>
    CompositionMorphism<A, B, C, AB, BC>
{
    pub fn new(a_to_b: AB, b_to_c: BC) -> Self {
        assert_eq!(a_to_b.range(), b_to_c.domain());
        Self {
            a: PhantomData,
            b: PhantomData,
            c: PhantomData,
            a_to_b,
            b_to_c,
        }
    }

    pub fn a(&self) -> &A {
        self.a_to_b.domain()
    }

    pub fn b(&self) -> &B {
        let b = self.a_to_b.range();
        debug_assert_eq!(b, self.b_to_c.domain());
        b
    }

    pub fn c(&self) -> &C {
        self.b_to_c.range()
    }

    pub fn a_to_b(&self) -> &AB {
        &self.a_to_b
    }

    pub fn b_to_c(&self) -> &BC {
        &self.b_to_c
    }
}

impl<A: Signature, B: Signature, C: Signature, AB: Morphism<A, B>, BC: Morphism<B, C>>
    Morphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn domain(&self) -> &A {
        self.a()
    }

    fn range(&self) -> &C {
        self.c()
    }
}

impl<A: SetSignature, B: SetSignature, C: SetSignature, AB: Function<A, B>, BC: Function<B, C>>
    Function<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn image(&self, x: &A::Set) -> C::Set {
        self.b_to_c.image(&self.a_to_b.image(x))
    }
}

impl<
    A: SetSignature,
    B: SetSignature,
    C: SetSignature,
    AB: InjectiveFunction<A, B>,
    BC: InjectiveFunction<B, C>,
> InjectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn try_preimage(&self, x: &C::Set) -> Option<A::Set> {
        self.a_to_b.try_preimage(&self.b_to_c.try_preimage(x)?)
    }
}

impl<
    A: SetSignature,
    B: SetSignature,
    C: SetSignature,
    AB: BijectiveFunction<A, B>,
    BC: BijectiveFunction<B, C>,
> BijectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn preimage(&self, x: &C::Set) -> A::Set {
        self.a_to_b.preimage(&self.b_to_c.preimage(x))
    }
}

pub trait MorphismsSignature<Domain: Signature, Range: Signature>: SetSignature
where
    Self::Set: Morphism<Domain, Range>,
{
}

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Functions<Domain: SetSignature, Range: SetSignature> {
    domain: Domain,
    range: Range,
}

impl<Domain: SetSignature, Range: SetSignature> Functions<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: SetSignature, Range: SetSignature> Signature for Functions<Domain, Range> {}

impl<Domain: FiniteSetSignature, Range: EqSignature> SetSignature for Functions<Domain, Range> {
    type Set = Vec<Range::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        x.len() == self.domain.size() && x.iter().all(|y| self.range.is_element(y))
    }
}

impl<Domain: FiniteSetSignature, Range: EqSignature + FiniteSetSignature> CountableSetSignature
    for Functions<Domain, Range>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
        (0..self.domain.size())
            .map(|_| self.range.list_all_elements())
            .multi_cartesian_product()
    }
}

impl<Domain: FiniteSetSignature, Range: EqSignature + FiniteSetSignature> FiniteSetSignature
    for Functions<Domain, Range>
{
}
