use crate::*;
use ambassador::delegatable_trait;

/// For anything which is a structure type
#[delegatable_trait]
pub trait Signature: Clone + Debug + Eq + Send + Sync {}

pub trait BorrowedElem<S>: Borrow<S> + Clone + Debug + Send + Sync {}
impl<S, BS: Borrow<S> + Clone + Debug + Send + Sync> BorrowedElem<S> for BS {}

/// Helper for borrowing a structure type
pub trait BorrowedStructure<S: Signature>: Borrow<S> + Clone + Debug + Eq + Send + Sync {}
impl<S: Signature, BS: Borrow<S> + Clone + Debug + Eq + Send + Sync> BorrowedStructure<S> for BS {}
