use super::*;
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait MorphismStructure<Domain: Structure, Range: Structure>: Structure {
    fn domain(&self) -> &Domain;
    fn range(&self) -> &Range;
}

pub trait FunctionStructure<Domain: SetStructure, Range: SetStructure>:
    MorphismStructure<Domain, Range>
{
    fn image(&self, x: &Domain::Set) -> Range::Set;
}

pub trait InjectiveFunctionStructure<Domain: SetStructure, Range: SetStructure>:
    FunctionStructure<Domain, Range>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<Domain::Set>;
}

pub trait BijectiveFunctionStructure<Domain: SetStructure, Range: SetStructure>:
    InjectiveFunctionStructure<Domain, Range>
{
    fn preimage(&self, x: &Range::Set) -> Domain::Set {
        self.try_preimage(x).unwrap()
    }
}

/// A general morphism between two objects Domain -> Range
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Morphism<Domain: Structure, Range: Structure> {
    domain: Domain,
    range: Range,
}

impl<Domain: Structure, Range: Structure> Morphism<Domain, Range> {
    pub fn new(domain: Domain, range: Range) -> Self {
        Self { domain, range }
    }
}

impl<Domain: Structure, Range: Structure> Structure for Morphism<Domain, Range> {}

impl<Domain: Structure, Range: Structure> MorphismStructure<Domain, Range>
    for Morphism<Domain, Range>
{
    fn domain(&self) -> &Domain {
        &self.domain
    }

    fn range(&self) -> &Range {
        &self.range
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

impl<X: Structure> MorphismStructure<X, X> for IdentityMorphism<X> {
    fn domain(&self) -> &X {
        &self.x
    }

    fn range(&self) -> &X {
        &self.x
    }
}

impl<X: SetStructure> FunctionStructure<X, X> for IdentityMorphism<X> {
    fn image(&self, x: &X::Set) -> X::Set {
        x.clone()
    }
}

impl<X: SetStructure> InjectiveFunctionStructure<X, X> for IdentityMorphism<X> {
    fn try_preimage(&self, x: &X::Set) -> Option<X::Set> {
        Some(x.clone())
    }
}

impl<X: SetStructure> BijectiveFunctionStructure<X, X> for IdentityMorphism<X> {}

/// The composition A -> B -> C of two morphisms A -> B and B -> C
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositionMorphism<
    A: Structure,
    B: Structure,
    C: Structure,
    AB: MorphismStructure<A, B>,
    BC: MorphismStructure<B, C>,
> {
    a: PhantomData<A>,
    b: PhantomData<B>,
    c: PhantomData<C>,
    // required to satisfy ab.range() == bc.domain()
    a_to_b: AB,
    b_to_c: BC,
}

impl<
    A: Structure,
    B: Structure,
    C: Structure,
    AB: MorphismStructure<A, B>,
    BC: MorphismStructure<B, C>,
> CompositionMorphism<A, B, C, AB, BC>
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

impl<
    A: Structure,
    B: Structure,
    C: Structure,
    AB: MorphismStructure<A, B>,
    BC: MorphismStructure<B, C>,
> Structure for CompositionMorphism<A, B, C, AB, BC>
{
}

impl<
    A: Structure,
    B: Structure,
    C: Structure,
    AB: MorphismStructure<A, B>,
    BC: MorphismStructure<B, C>,
> MorphismStructure<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn domain(&self) -> &A {
        self.a()
    }

    fn range(&self) -> &C {
        self.c()
    }
}

impl<
    A: SetStructure,
    B: SetStructure,
    C: SetStructure,
    AB: FunctionStructure<A, B>,
    BC: FunctionStructure<B, C>,
> FunctionStructure<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn image(&self, x: &A::Set) -> C::Set {
        self.b_to_c.image(&self.a_to_b.image(x))
    }
}

impl<
    A: SetStructure,
    B: SetStructure,
    C: SetStructure,
    AB: InjectiveFunctionStructure<A, B>,
    BC: InjectiveFunctionStructure<B, C>,
> InjectiveFunctionStructure<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn try_preimage(&self, x: &C::Set) -> Option<A::Set> {
        self.a_to_b.try_preimage(&self.b_to_c.try_preimage(x)?)
    }
}

impl<
    A: SetStructure,
    B: SetStructure,
    C: SetStructure,
    AB: BijectiveFunctionStructure<A, B>,
    BC: BijectiveFunctionStructure<B, C>,
> BijectiveFunctionStructure<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn preimage(&self, x: &C::Set) -> A::Set {
        self.a_to_b.preimage(&self.b_to_c.preimage(x))
    }
}

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct FunctionFiniteDomain<Domain: FiniteSetStructure, Range: SetStructure> {
//     domain: Domain,
//     range: Range,
//     // func.len() == domain.size()
//     func: Vec<Range>,
// }

// impl<Domain: FiniteSetStructure, Range: SetStructure> Structure
//     for FunctionFiniteDomain<Domain, Range>
// {
// }

// impl<Domain: FiniteSetStructure, Range: SetStructure> MorphismStructure<Domain, Range>
//     for FunctionFiniteDomain<Domain, Range>
// {
//     fn domain(&self) -> &Domain {
//         &self.domain
//     }

//     fn range(&self) -> &Range {
//         &self.range
//     }
// }

// impl<Domain: FiniteSetStructure, Range: SetStructure> SetStructure for Morphism<Domain, Range> {
//     type Set = Vec<Domain>;
// }

// impl<Domain: FiniteSetStructure, Range: SetStructure>  FiniteSetStructure for Morphism<Domain, Range>
