use super::*;
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait Morphism<Domain: Structure, Range: Structure>: Structure {
    fn domain(&self) -> &Domain;
    fn range(&self) -> &Range;
}

pub trait Function<Domain: SetStructure, Range: SetStructure>: Morphism<Domain, Range> {
    fn image(&self, x: &Domain::Set) -> Range::Set;
}

pub trait InjectiveFunction<Domain: SetStructure, Range: SetStructure>:
    Function<Domain, Range>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<Domain::Set>;
}

pub trait BijectiveFunction<Domain: SetStructure, Range: SetStructure>:
    InjectiveFunction<Domain, Range>
{
    fn preimage(&self, x: &Range::Set) -> Domain::Set {
        self.try_preimage(x).unwrap()
    }
}

/// The identity morphism X -> X
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityMorphism<X: Structure> {
    x: X,
}

impl<X: Structure> IdentityMorphism<X> {
    pub fn new(x: X) -> Self {
        Self { x }
    }
}

impl<X: Structure> Structure for IdentityMorphism<X> {}

impl<X: Structure> Morphism<X, X> for IdentityMorphism<X> {
    fn domain(&self) -> &X {
        &self.x
    }

    fn range(&self) -> &X {
        &self.x
    }
}

impl<X: SetStructure> Function<X, X> for IdentityMorphism<X> {
    fn image(&self, x: &X::Set) -> X::Set {
        x.clone()
    }
}

impl<X: SetStructure> InjectiveFunction<X, X> for IdentityMorphism<X> {
    fn try_preimage(&self, x: &X::Set) -> Option<X::Set> {
        Some(x.clone())
    }
}

impl<X: SetStructure> BijectiveFunction<X, X> for IdentityMorphism<X> {}

/// The composition A -> B -> C of two morphisms A -> B and B -> C
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionMorphism<
    A: Structure,
    B: Structure,
    C: Structure,
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

impl<A: Structure, B: Structure, C: Structure, AB: Morphism<A, B>, BC: Morphism<B, C>>
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

impl<A: Structure, B: Structure, C: Structure, AB: Morphism<A, B>, BC: Morphism<B, C>> Structure
    for CompositionMorphism<A, B, C, AB, BC>
{
}

impl<A: Structure, B: Structure, C: Structure, AB: Morphism<A, B>, BC: Morphism<B, C>>
    Morphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn domain(&self) -> &A {
        self.a()
    }

    fn range(&self) -> &C {
        self.c()
    }
}

impl<A: SetStructure, B: SetStructure, C: SetStructure, AB: Function<A, B>, BC: Function<B, C>>
    Function<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn image(&self, x: &A::Set) -> C::Set {
        self.b_to_c.image(&self.a_to_b.image(x))
    }
}

impl<
    A: SetStructure,
    B: SetStructure,
    C: SetStructure,
    AB: InjectiveFunction<A, B>,
    BC: InjectiveFunction<B, C>,
> InjectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn try_preimage(&self, x: &C::Set) -> Option<A::Set> {
        self.a_to_b.try_preimage(&self.b_to_c.try_preimage(x)?)
    }
}

impl<
    A: SetStructure,
    B: SetStructure,
    C: SetStructure,
    AB: BijectiveFunction<A, B>,
    BC: BijectiveFunction<B, C>,
> BijectiveFunction<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn preimage(&self, x: &C::Set) -> A::Set {
        self.a_to_b.preimage(&self.b_to_c.preimage(x))
    }
}
