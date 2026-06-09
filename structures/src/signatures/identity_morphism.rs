use crate::*;

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
    fn image(&self, x: &X::Elem) -> X::Elem {
        x.clone()
    }
}

impl<X: SetSignature> InjectiveFunction<X, X> for IdentityMorphism<X> {
    fn try_preimage(&self, x: &X::Elem) -> Option<X::Elem> {
        Some(x.clone())
    }
}

impl<X: SetSignature> BijectiveFunction<X, X> for IdentityMorphism<X> {}
