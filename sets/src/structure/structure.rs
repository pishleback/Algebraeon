use std::{borrow::Borrow, fmt::Debug};

pub trait Structure: Clone + Debug + PartialEq + Eq {}

/// Instances of a type implementing this trait represent
/// a set of elements of type `Self::Set` with some
/// structure, for example, the structure of a ring.
pub trait SetStructure: Structure {
    type Set: Clone + Debug;
}

pub trait MetaType: Clone + Debug {
    type Structure: SetStructure<Set = Self>;
    fn structure() -> Self::Structure;
}

pub fn common_structure<S: SetStructure>(
    structure1: impl Borrow<S>,
    structure2: impl Borrow<S>,
) -> S {
    if structure1.borrow() == structure2.borrow() {
        structure1.borrow().clone()
    } else {
        panic!("Unequal ring structures")
    }
}

pub trait ToStringStructure: SetStructure {
    fn to_string(&self, elem: &Self::Set) -> String;
}

pub trait EqStructure: SetStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
}

pub trait CountableSetStructure: SetStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set>;
}

pub trait FiniteSetStructure: CountableSetStructure {
    fn list_all_elements(&self) -> Vec<Self::Set> {
        self.generate_all_elements().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_canonical_structure_derive::CanonicalStructure;

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

        impl ToStringStructure for ACanonicalStructure {
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

        impl Structure for A {}

        impl SetStructure for A {
            type Set = usize;
        }

        impl ToStringStructure for A {
            fn to_string(&self, elem: &Self::Set) -> String {
                elem.to_string()
            }
        }
    }
}
