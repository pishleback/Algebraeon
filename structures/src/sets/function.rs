use crate::*;

pub trait Function<Domain: SetSignature, Range: SetSignature>: Morphism<Domain, Range> {
    fn image(&self, x: &Domain::Elem) -> Range::Elem;
}

/// A function from a set into itself
pub trait Endofunction<X: SetSignature + EqSignature>: Function<X, X> {
    // TODO: remove EqSignature requirement and use specialization once it is stable.
    /// check if an element is fixed
    fn is_fixed_point(&self, x: X::Elem) -> bool {
        self.domain().equal(&self.image(&x), &x)
    }
}

pub trait InjectiveFunction<Domain: SetSignature, Range: SetSignature>:
    Function<Domain, Range>
{
    fn try_preimage(&self, y: &Range::Elem) -> Option<Domain::Elem>;
}

pub trait BijectiveFunction<Domain: SetSignature, Range: SetSignature>:
    InjectiveFunction<Domain, Range>
{
    fn preimage(&self, y: &Range::Elem) -> Domain::Elem {
        self.try_preimage(y).unwrap()
    }
}
