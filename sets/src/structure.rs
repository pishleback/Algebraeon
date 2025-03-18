use std::{borrow::Borrow, fmt::Debug, marker::PhantomData, rc::Rc};

pub trait Structure: Clone + Debug + PartialEq + Eq {
    type Set: Clone + Debug;
}

#[derive(Debug, Clone)]
pub struct CannonicalStructure<T: MetaType> {
    _ghost: std::marker::PhantomData<T>,
}
impl<T: MetaType> PartialEq for CannonicalStructure<T> {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}
impl<T: MetaType> Eq for CannonicalStructure<T> {}
pub trait MetaType: Clone + Debug {
    type Structure: Structure<Set = Self>;
    fn structure() -> Rc<Self::Structure>;
}

impl<T: MetaType> CannonicalStructure<T> {
    pub fn new() -> Self {
        Self {
            _ghost: PhantomData::default(),
        }
    }
}
impl<T: MetaType> Structure for CannonicalStructure<T> {
    type Set = T;
}

pub fn common_structure<S: Structure>(
    structure1: impl Borrow<Rc<S>>,
    structure2: impl Borrow<Rc<S>>,
) -> Rc<S> {
    if structure1.borrow() == structure2.borrow() {
        structure1.borrow().clone()
    } else {
        panic!("Unequal ring structures")
    }
}

pub trait ToStringStructure: Structure {
    fn to_string(&self, elem: &Self::Set) -> String;
}
impl<T: MetaType + ToString> ToStringStructure for CannonicalStructure<T> {
    fn to_string(&self, elem: &Self::Set) -> String {
        elem.to_string()
    }
}

pub trait PartialEqStructure: Structure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
}
impl<T: MetaType + PartialEq> PartialEqStructure for CannonicalStructure<T> {
    fn equal(&self, a: &T, b: &T) -> bool {
        a == b
    }
}

pub trait EqStructure: PartialEqStructure {}
impl<T: MetaType + Eq> EqStructure for CannonicalStructure<T> {}

// impl MetaType for Natural {
//     type Structure = CannonicalStructure<Natural>;

//     fn structure() -> Rc<Self::Structure> {
//         CannonicalStructure::new().into()
//     }
// }

// impl MetaType for Integer {
//     type Structure = CannonicalStructure<Integer>;

//     fn structure() -> Rc<Self::Structure> {
//         CannonicalStructure::new().into()
//     }
// }

// impl MetaType for Rational {
//     type Structure = CannonicalStructure<Rational>;

//     fn structure() -> Rc<Self::Structure> {
//         CannonicalStructure::new().into()
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cannonical_structure() {
        #[derive(Debug, Clone, PartialEq, Eq)]
        struct A {
            x: i32,
        }

        impl MetaType for A {
            type Structure = CannonicalStructure<A>;

            fn structure() -> Rc<Self::Structure> {
                CannonicalStructure::new().into()
            }
        }

        impl ToString for A {
            fn to_string(&self) -> String {
                self.x.to_string()
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
        struct Scary {
            t: usize,
        }

        impl Structure for Scary {
            type Set = usize;
        }

        impl ToStringStructure for Scary {
            fn to_string(&self, elem: &Self::Set) -> String {
                elem.to_string()
            }
        }
    }
}
