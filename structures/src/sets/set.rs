use crate::*;

use algebraeon_macros::{signature_meta_trait, skip_meta};
use ambassador::delegatable_trait;
use paste::paste;
use rand::{Rng, RngExt, SeedableRng, rngs::StdRng};
use std::fmt::Debug;

/// Instances of a type implementing this trait represent
/// a set of elements of type `Self::Elem` with some
/// structure, for example, the structure of a ring.
#[signature_meta_trait]
#[delegatable_trait]
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

pub trait MetaType: Clone + Debug {
    type Signature: SetSignature<Elem = Self> + 'static;

    fn structure() -> Self::Signature;
}

pub trait MetaTypeRef: MetaType {
    fn structure_ref() -> &'static Self::Signature;
}

#[signature_meta_trait]
#[delegatable_trait]
pub trait ToStringSignature: SetSignature {
    fn to_string(&self, elem: &Self::Elem) -> String;
}

#[signature_meta_trait]
#[delegatable_trait]
pub trait EqSignature: SetSignature {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool;
}

#[delegatable_trait]
pub trait CountableSetSignature: SetSignature {
    /// Yield distinct elements of the set such that every element eventually appears.
    /// Always yields elements in the same order.
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem>;
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem>;
}

/// A set with finitely many elements
#[signature_meta_trait]
#[delegatable_trait]
pub trait FiniteSetSignature: CountableSetSignature {
    /// A list of all elements in the set.
    /// Must always return elements in the same order.
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.generate_all_elements().collect()
    }

    fn size(&self) -> Natural {
        Natural::from(self.list_all_elements().len())
    }

    #[skip_meta]
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        let rng = StdRng::seed_from_u64(seed);
        FiniteSetRandomElementGenerator::<Self, StdRng> {
            all_elements: self.list_all_elements(),
            rng,
        }
    }
}
make_maybe_trait!(FiniteSet);

/// A set with N elements
#[signature_meta_trait]
pub trait ConstSizeFiniteSetSignature<const N: usize>: FiniteSetSignature {
    fn list_all_elements_sized(&self) -> [Self::Elem; N] {
        self.list_all_elements().try_into().unwrap()
    }
}

/// A finite set where the elements are numbered 0, 1, ..., n-1
/// self.list_all_elements is required to return elements in the correct order
/// The ordering on the set must also agree with the ordering given by the enumeration
#[signature_meta_trait]
#[delegatable_trait]
pub trait OrderedFiniteSetSignature: FiniteSetSignature + OrdSignature {
    /// List all elements in the order in which they are numbered
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem>;

    /// Return the numbering of an element
    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural;

    /// Return the numbering of an element
    /// None iff number is too large
    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem>;
}

// for testing the invariants of OrderedFiniteSetSignature
#[macro_export]
macro_rules! assert_enumerated_ord_finite_set {
    (
        $set:expr,
        $size:expr
    ) => {{
        let size = $size as usize;
        let set = $set;
        let elements = set.list_all_elements_ordered();

        assert_eq!(elements.len(), size);
        assert_eq!(set.size(), Natural::from(size));

        // all elements are valid
        for elem in &elements {
            println!("{:?}", elem);
            assert!(set.validate_element(elem).is_ok());
        }

        // all elements are distinct and ordered
        for i in 0..size {
            for j in (i + 1)..size {
                let si = &elements[i];
                let sj = &elements[j];
                assert!(set.cmp(si, sj).is_lt());
            }
        }

        // enumeration is correct
        for (i, s) in elements.iter().enumerate() {
            assert_eq!(Natural::from(i), set.element_to_enumeration(s));
            assert!(set.equal(&set.enumeration_to_element(&Natural::from(i)).unwrap(), s));
        }

        // enumeration past the end is invalid
        assert!(set.enumeration_to_element(&Natural::from(size)).is_none());
    }};
}

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

        #[allow(clippy::to_string_trait_impl)]
        impl ToString for A {
            fn to_string(&self) -> String {
                self.x.to_string()
            }
        }

        impl ToStringSignature for ACanonicalStructure {
            fn to_string(&self, elem: &Self::Elem) -> String {
                ToString::to_string(elem)
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
