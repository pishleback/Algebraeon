use core::fmt::Debug;
use std::{borrow::Borrow, fmt::Display, marker::PhantomData, rc::Rc};

//Instances of this represent structure on Set represented by a type
//For example, instances of this might represent a group structure on a Set
pub trait Structure: Debug + Clone + PartialEq + Eq + 'static {
    type Set: Clone + Debug;
}

pub trait DisplayableStructure: Structure {
    fn elem_to_string(&self, elem: &Self::Set) -> String;
}

pub trait EqualityStructure: Structure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool;
}

pub trait InfiniteStructure: Structure {
    //generate an infinite sequence of distinct elements
    fn generate_distinct_elements(&self) -> Box<dyn Iterator<Item = Self::Set>>;
}

pub fn common_structure<RS: Structure>(
    structure1: impl Borrow<Rc<RS>>,
    structure2: impl Borrow<Rc<RS>>,
) -> Rc<RS> {
    if structure1.borrow() == structure2.borrow() {
        structure1.borrow().clone()
    } else {
        panic!("Unequal ring structures")
    }
}

pub trait StructuredType: Clone + Debug + 'static {
    type Structure: Structure<Set = Self>;

    fn structure() -> Rc<Self::Structure>;

    fn into_ring(self) -> StructuredElement<Self::Structure> {
        StructuredElement::new(Self::structure(), self)
    }
}

#[derive(Debug, Clone)]
pub struct StructuredElement<S: Structure> {
    structure: Rc<S>,
    elem: S::Set,
}

impl<S: Structure> StructuredElement<S> {
    pub fn new(structure: Rc<S>, elem: S::Set) -> Self {
        Self { structure, elem }
    }

    pub fn structure(&self) -> Rc<S> {
        self.structure.clone()
    }

    pub fn ref_set(&self) -> &S::Set {
        &self.elem
    }

    pub fn into_set(self) -> S::Set {
        self.elem
    }
}

impl<RS: DisplayableStructure> Display for StructuredElement<RS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.structure().elem_to_string(&self.elem))
    }
}

#[derive(Debug, Clone)]
pub struct CannonicalStructure<T: StructuredType> {
    _ghost: PhantomData<T>,
}

impl<T: StructuredType> PartialEq for CannonicalStructure<T> {
    fn eq(&self, _other: &Self) -> bool {
        true
    }
}

impl<T: StructuredType> Eq for CannonicalStructure<T> {}

impl<T: StructuredType> CannonicalStructure<T> {
    pub fn new() -> Self {
        Self {
            _ghost: PhantomData::default(),
        }
    }
}

impl<T: StructuredType> Structure for CannonicalStructure<T> {
    type Set = T;
}

impl<T: StructuredType + Display> DisplayableStructure for CannonicalStructure<T> {
    fn elem_to_string(&self, elem: &Self::Set) -> String {
        elem.to_string()
    }
}
