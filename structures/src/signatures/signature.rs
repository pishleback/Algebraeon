use std::ops::Deref;

use crate::*;

/// For anything which is a structure type
pub trait Signature: Clone + Debug + Eq + Send + Sync {}

pub trait BorrowedElem<S>: Borrow<S> + Clone + Debug + Send + Sync {}
impl<S, BS: Borrow<S> + Clone + Debug + Send + Sync> BorrowedElem<S> for BS {}

// /// Helper for borrowing a structure type
// pub trait BorrowedStructure<S: Signature>: Borrow<S> + Clone + Debug + Eq + Send + Sync {}
// impl<S: Signature, BS: Borrow<S> + Clone + Debug + Eq + Send + Sync> BorrowedStructure<S> for BS {}

pub enum Arg<'a, T> {
    Borrowed(&'a T),
    Owned(T),
}

impl<'a, T> From<T> for Arg<'a, T> {
    fn from(value: T) -> Self {
        Self::Owned(value)
    }
}

impl<'a, T> From<&'a T> for Arg<'a, T> {
    fn from(value: &'a T) -> Self {
        Self::Borrowed(value)
    }
}

impl<'a, T> From<&'a mut T> for Arg<'a, T> {
    fn from(value: &'a mut T) -> Self {
        Self::Borrowed(value)
    }
}

impl<'a, T> AsRef<T> for Arg<'a, T> {
    fn as_ref(&self) -> &T {
        match self {
            Arg::Borrowed(value) => value,
            Arg::Owned(value) => value,
        }
    }
}

impl<'a, T> Deref for Arg<'a, T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

impl<'a, T> Arg<'a, T> {
    pub fn into_owned(self) -> T
    where
        T: Clone,
    {
        match self {
            Arg::Owned(x) => x,
            Arg::Borrowed(x) => x.clone(),
        }
    }
}
