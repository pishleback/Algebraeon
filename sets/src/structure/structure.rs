use paste::paste;
use rand::{Rng, RngExt, SeedableRng, rngs::StdRng};
use std::{borrow::Borrow, fmt::Debug};

macro_rules! make_maybe_trait {
    ($name:ident) => {
        paste! {
            pub trait [<Maybe $name Signature>]: SetSignature {
                type [<$name Structure>]: [<$name Signature>] <Elem = Self::Elem>;

                #[allow(clippy::result_unit_err)]
                fn [<$name:snake _structure>](
                    &self
                ) -> Result<Self::[<$name Structure>], ()>;
            }
        }
    };
}

pub trait Signature: Clone + Debug + Eq + Send + Sync {}

pub trait MetaType: Clone + Debug {
    type Signature: SetSignature<Elem = Self>;
    fn structure() -> Self::Signature;
}

pub trait BorrowedSet<S>: Borrow<S> + Clone + Debug + Send + Sync {}
impl<S, BS: Borrow<S> + Clone + Debug + Send + Sync> BorrowedSet<S> for BS {}

pub trait BorrowedStructure<S: Signature>: Borrow<S> + Clone + Debug + Eq + Send + Sync {}
impl<S: Signature, BS: Borrow<S> + Clone + Debug + Eq + Send + Sync> BorrowedStructure<S> for BS {}

/// Instances of a type implementing this trait represent
/// a set of elements of type `Self::Elem` with some
/// structure, for example, the structure of a ring.
pub trait SetSignature: Signature {
    type Elem: Clone + Debug + Send + Sync;

    /// Some instances of `Self::Elem` may not be valid to represent elements of this set.
    /// Return `Ok(())` if `x` is a valid element and an `Err` explaining why if not.
    fn validate_element(&self, x: &Self::Elem) -> Result<(), String>;

    /// Some instances of `Self::Elem` may not be valid to represent elements of this set.
    /// Return `true` if `x` is a valid element and an `false` if not.
    fn is_element(&self, x: &Self::Elem) -> bool {
        self.validate_element(x).is_ok()
    }
}

pub trait ToStringSignature: SetSignature {
    fn to_string(&self, elem: &Self::Elem) -> String;
}

pub trait EqSignature: SetSignature {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool;
}

pub trait CountableSetSignature: SetSignature {
    /// Yield distinct elements of the set such that every element eventually appears.
    /// Always yields elements in the same order.
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone;
}

pub trait FiniteSetSignature: CountableSetSignature {
    /// A list of all elements in the set.
    /// Always returns elements in the same order.
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.generate_all_elements().collect()
    }
    fn size(&self) -> usize {
        self.list_all_elements().len()
    }
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
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
    all_elements: Vec<S::Elem>,
    rng: R,
}

impl<S: FiniteSetSignature, R: Rng> Iterator for FiniteSetRandomElementGenerator<S, R> {
    type Item = S::Elem;

    fn next(&mut self) -> Option<Self::Item> {
        if self.all_elements.is_empty() {
            None
        } else {
            let idx = self.rng.random_range(0..self.all_elements.len());
            Some(self.all_elements[idx].clone())
        }
    }
}

/// Instances of a type implementing this trait represent
/// a set formed by a quotient of another set.
pub trait QuotientSetSignature<PreQuoSet: SetSignature>: SetSignature {
    fn pre_quotient_set(&self) -> &PreQuoSet;

    fn project(&self, x: PreQuoSet::Elem) -> Self::Elem;
    fn project_ref(&self, x: &PreQuoSet::Elem) -> Self::Elem;

    /// Return an element of the pre-quotient set which projects to the given element.
    fn unproject(&self, x: Self::Elem) -> PreQuoSet::Elem;
    fn unproject_ref(&self, x: &Self::Elem) -> PreQuoSet::Elem;
}

/// A quotient set where elements are represented using representative elements of the pre-quotient set
pub trait QuotientSetRepresentativesSignature<PreQuoSet: SetSignature<Elem = Self::Elem>>:
    QuotientSetSignature<PreQuoSet>
{
    /// Must satisfy x = y in the quotient set iff reduced_representative(x) = reduced_representative(y) in the pre quotient set.
    fn reduced_representative(&self, x: &Self::Elem) -> Self::Elem;
}

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
            fn to_string(&self, elem: &Self::Elem) -> String {
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
            type Elem = usize;

            fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
                Ok(())
            }
        }

        impl ToStringSignature for A {
            fn to_string(&self, elem: &Self::Elem) -> String {
                elem.to_string()
            }
        }
    }
}
