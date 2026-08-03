use crate::*;
use std::{marker::PhantomData, rc::Rc};

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
    a_to_b: Rc<AB>,
    b_to_c: Rc<BC>,
}

impl<A: Signature, B: Signature, C: Signature, AB: Morphism<A, B>, BC: Morphism<B, C>>
    CompositionMorphism<A, B, C, AB, BC>
{
    pub fn new(a_to_b: Rc<AB>, b_to_c: Rc<BC>) -> Rc<Self> {
        assert_eq!(a_to_b.range(), b_to_c.domain());
        Self {
            a: PhantomData,
            b: PhantomData,
            c: PhantomData,
            a_to_b,
            b_to_c,
        }
        .into()
    }

    pub fn a(&self) -> Rc<A> {
        self.a_to_b.domain()
    }

    pub fn b(&self) -> Rc<B> {
        let b = self.a_to_b.range();
        debug_assert_eq!(b, self.b_to_c.domain());
        b
    }

    pub fn c(&self) -> Rc<C> {
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
    fn domain(self: &Rc<Self>) -> Rc<A> {
        self.a()
    }

    fn range(self: &Rc<Self>) -> Rc<C> {
        self.c()
    }
}

impl<
    A: SetSignature,
    B: SetSignature,
    C: SetSignature,
    AB: FunctionMorphism<A, B>,
    BC: FunctionMorphism<B, C>,
> FunctionMorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn image(self: &Rc<Self>, x: &A::Elem) -> C::Elem {
        self.b_to_c.image(&self.a_to_b.image(x))
    }
}

impl<
    A: SetSignature,
    B: SetSignature,
    C: SetSignature,
    AB: InjectiveFunctionMorphism<A, B>,
    BC: InjectiveFunctionMorphism<B, C>,
> InjectiveFunctionMorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn try_preimage(self: &Rc<Self>, x: &C::Elem) -> Option<A::Elem> {
        self.a_to_b.try_preimage(&self.b_to_c.try_preimage(x)?)
    }
}

impl<
    A: SetSignature,
    B: SetSignature,
    C: SetSignature,
    AB: BijectiveFunctionMorphism<A, B>,
    BC: BijectiveFunctionMorphism<B, C>,
> BijectiveFunctionMorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
{
    fn preimage(self: &Rc<Self>, x: &C::Elem) -> A::Elem {
        self.a_to_b.preimage(&self.b_to_c.preimage(x))
    }
}
