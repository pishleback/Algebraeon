use crate::*;
use std::sync::Arc;

pub trait Morphism<Domain: Signature, Range: Signature>: Debug + Clone + Send + Sync {
    fn domain(self: &Arc<Self>) -> Arc<Domain>;
    fn range(self: &Arc<Self>) -> Arc<Range>;
}

/// A morphism from an object to itself
pub trait Endomorphism<X: Signature>: Morphism<X, X> {}
impl<X: Signature, T: Morphism<X, X>> Endomorphism<X> for T {}

// pub trait BorrowedMorphism<Domain: Signature, Range: Signature, M: Morphism<Domain, Range>>:
//     Borrow<M> + Clone + Debug + Send + Sync
// {
// }
// impl<
//     Domain: Signature,
//     Range: Signature,
//     M: Morphism<Domain, Range>,
//     BM: Borrow<M> + Clone + Debug + Send + Sync,
// > BorrowedMorphism<Domain, Range, M> for BM
// {
// }
