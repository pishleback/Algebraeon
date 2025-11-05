use paste::paste;
use rand::{Rng, SeedableRng, rngs::StdRng};
use std::{borrow::Borrow, fmt::Debug};

macro_rules! make_maybe_trait {
    ($name:ident) => {
        paste! {
            pub trait [<Maybe $name Signature>]: SetSignature {
                type [<$name Structure>]: [<$name Signature>] <Set = Self::Set>;

                #[allow(clippy::result_unit_err)]
                fn [<$name:snake _structure>](
                    &self
                ) -> Result<Self::[<$name Structure>], ()>;
            }
        }
    };
}

pub trait Signature: Clone + Debug + Eq + Send + Sync {}

/// Instances of a type implementing this trait represent
/// a set of elements of type `Self::Set` with some
/// structure, for example, the structure of a ring.
pub trait SetSignature: Signature {
    type Set: Clone + Debug + Send + Sync;

    /// Some instances of `Self::Set` may not be valid to represent elements of this set.
    /// Return `true` if `x` is a valid element and `false` if not.
    fn is_element(&self, x: &Self::Set) -> Result<(), String>;
}

pub trait MetaType: Clone + Debug {
    type Signature: SetSignature<Set = Self>;
    fn structure() -> Self::Signature;
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
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone;
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
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Set> + Clone {
        let rng = StdRng::seed_from_u64(seed);
        FiniteSetRandomElementGenerator::<Self, StdRng> {
            all_elements: self.list_all_elements(),
            rng,
        }
    }
}
make_maybe_trait!(FiniteSet);

#[derive(Debug, Clone)]
pub struct FiniteSetRandomElementGenerator<S: FiniteSetSignature, R: Rng> {
    all_elements: Vec<S::Set>,
    rng: R,
}

impl<S: FiniteSetSignature, R: Rng> Iterator for FiniteSetRandomElementGenerator<S, R> {
    type Item = S::Set;

    fn next(&mut self) -> Option<Self::Item> {
        if self.all_elements.is_empty() {
            None
        } else {
            let idx = self.rng.random_range(0..self.all_elements.len());
            Some(self.all_elements[idx].clone())
        }
    }
}

pub trait BorrowedStructure<S: Signature>: Borrow<S> + Clone + Debug + Eq + Send + Sync {}
impl<S: Signature, BS: Borrow<S> + Clone + Debug + Eq + Send + Sync> BorrowedStructure<S> for BS {}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_macros::CanonicalStructure;

    #[test]
    fn canonical_structure() {
        #[derive(Debug, Clone, PartialEq, Eq, CanonicalStructure)]
        #[canonical_structure(eq)]
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
        assert!(!v);
        println!("{}", A::structure().to_string(&a));
    }

    #[test]
    fn to_string_structure_impl() {
        #[allow(dead_code)]
        #[derive(Debug, Clone, PartialEq, Eq)]
        struct A {
            t: usize,
        }

        impl Signature for A {}

        impl SetSignature for A {
            type Set = usize;

            fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
                Ok(())
            }
        }

        impl ToStringSignature for A {
            fn to_string(&self, elem: &Self::Set) -> String {
                elem.to_string()
            }
        }
    }
}
