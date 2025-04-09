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

// #[derive(Debug, Clone)]
// pub struct CanonicalStructure<T: MetaType> {
//     _ghost: std::marker::PhantomData<T>,
// }

// impl<T: MetaType> CanonicalStructure<T> {
//     pub fn new() -> Self {
//         Self {
//             _ghost: PhantomData::default(),
//         }
//     }
// }

// impl<T: MetaType> Structure for CanonicalStructure<T> {}

// impl<T: MetaType> SetStructure for CanonicalStructure<T> {
//     type Set = T;
// }

// impl<T: MetaType> PartialEq for CanonicalStructure<T> {
//     fn eq(&self, _: &Self) -> bool {
//         true
//     }
// }

// impl<T: MetaType> Eq for CanonicalStructure<T> {}

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
// impl<T: MetaType + ToString> ToStringStructure for CanonicalStructure<T> {
//     fn to_string(&self, elem: &Self::Set) -> String {
//         elem.to_string()
//     }
// }

pub trait PartialEqStructure: SetStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
}
// impl<T: MetaType + PartialEq> PartialEqStructure for CanonicalStructure<T> {
//     fn equal(&self, a: &T, b: &T) -> bool {
//         a == b
//     }
// }

pub trait EqStructure: PartialEqStructure {}
// impl<T: MetaType + Eq> EqStructure for CanonicalStructure<T> {}

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
    use algebraeon_canonical_structure_derive::CanonicalStructure;

    use super::*;

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

        impl PartialEqStructure for ACanonicalStructure {
            fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
                a == b
            }
        }

        impl EqStructure for ACanonicalStructure {}

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
