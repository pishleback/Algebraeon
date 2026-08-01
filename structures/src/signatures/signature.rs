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

// Future thing: replace BorrowedStructure with something like this:
//  - Need either this https://doc.rust-lang.org/stable/std/ops/trait.Receiver.html or use Arc directly

// #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
// pub struct WrappedStructure<S: Signature>(std::sync::Arc<S>);

// impl<S: Signature> From<S> for WrappedStructure<S> {
//     fn from(s: S) -> Self {
//         Self(std::sync::Arc::new(s))
//     }
// }

// impl<S: Signature> AsRef<S> for WrappedStructure<S> {
//     fn as_ref(&self) -> &S {
//         &self.0
//     }
// }

// impl<S: Signature> Borrow<S> for WrappedStructure<S> {
//     fn borrow(&self) -> &S {
//         &self.0
//     }
// }
