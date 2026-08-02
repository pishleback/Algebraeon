use crate::*;
use std::{marker::PhantomData, sync::Arc};

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

    pub fn a(&self) -> Arc<A> {
        self.a_to_b.domain()
    }

    pub fn b(&self) -> Arc<B> {
        let b = self.a_to_b.range();
        debug_assert_eq!(b, self.b_to_c.domain());
        b
    }

    pub fn c(&self) -> Arc<C> {
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
    fn domain(&self) -> Arc<A> {
        self.a()
    }

    fn range(&self) -> Arc<C> {
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
    fn image(&self, x: &A::Elem) -> C::Elem {
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
    fn try_preimage(&self, x: &C::Elem) -> Option<A::Elem> {
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
    fn preimage(&self, x: &C::Elem) -> A::Elem {
        self.a_to_b.preimage(&self.b_to_c.preimage(x))
    }
}
