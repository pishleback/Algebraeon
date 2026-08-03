use crate::*;

/// For anything which is a structure type
pub trait Signature: Clone + Debug + Eq {}

pub trait BorrowedElem<S>: Borrow<S> + Clone + Debug {}
impl<S, BS: Borrow<S> + Clone + Debug> BorrowedElem<S> for BS {}
