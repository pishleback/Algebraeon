use std::{borrow::Borrow, fmt::Debug};

pub trait Signature: Clone + Debug + PartialEq + Eq {}

/// Instances of a type implementing this trait represent
/// a set of elements of type `Self::Set` with some
/// structure, for example, the structure of a ring.
pub trait SetSignature: Signature {
    type Set: Clone + Debug;

    /// Some instances of Self::Set may not be valid to represent elements of this set.
    /// Return `true` if `x` is a valid element and `false` if not.
    fn is_element(&self, x: &Self::Set) -> bool;
}

pub trait MetaType: Clone + Debug {
    type Signature: SetSignature<Set = Self>;
    fn structure() -> Self::Signature;
}

pub fn common_structure<S: Signature>(structure1: impl Borrow<S>, structure2: impl Borrow<S>) -> S {
    if structure1.borrow() == structure2.borrow() {
        structure1.borrow().clone()
    } else {
        panic!("Unequal ring structures")
    }
}

pub trait ToStringSignature: SetSignature {
    fn to_string(&self, elem: &Self::Set) -> String;
}

pub trait EqSignature: SetSignature {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
}

pub trait CountableSetSignature: SetSignature {
    /// Yield distinct elements of the set such that every element eventually appears.
    /// Always yields elements in the same order.
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set>;
}

pub trait FiniteSetSignature: CountableSetSignature {
    /// A list of all elements in the set.
    /// Always returns elements in the same order.
    fn list_all_elements(&self) -> Vec<Self::Set> {
        self.generate_all_elements().collect()
    }
    fn size(&self) -> usize {
        self.list_all_elements().len()
    }
}

pub trait BorrowedStructure<S: Signature>: Borrow<S> + Clone + std::fmt::Debug + PartialEq + Eq {}
impl<S: Signature, BS: Borrow<S> + Clone + std::fmt::Debug + PartialEq + Eq> BorrowedStructure<S> for BS {}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_macros::CanonicalStructure;

    #[test]
    fn canonical_structure() {
        #[derive(Debug, Clone, PartialEq, Eq, CanonicalStructure)]
        pub struct A {
            x: i32,
        }

        impl ToString for A {
            fn to_string(&self) -> String {
                self.x.to_string()
            }
        }

        impl ToStringSignature for ACanonicalStructure {
            fn to_string(&self, elem: &Self::Set) -> String {
                elem.to_string()
            }
        }

        let a = A { x: 3 };
        let b = A { x: 4 };
        let v = A::structure().equal(&a, &b);
        assert_eq!(v, false);
        println!("{}", A::structure().to_string(&a));
    }

    #[test]
    fn foo() {
        #[derive(Debug, Clone, PartialEq, Eq)]
        struct A {
            t: usize,
        }

        impl Signature for A {}

        impl SetSignature for A {
            type Set = usize;

            fn is_element(&self, _x: &Self::Set) -> bool {
                true
            }
        }

        impl ToStringSignature for A {
            fn to_string(&self, elem: &Self::Set) -> String {
                elem.to_string()
            }
        }
    }
}
