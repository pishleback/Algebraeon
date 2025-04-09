# The Structure System

## The Structure Trait

In mathematics there are many instances of sets of sets and often these sets of sets parameterized by other sets. For example the set of integers modulo \\(n\\) is defined for every natural number \\(n \in \mathbb{N}\\). A mathematician would like to think that there exists infinitely many types, \\(\frac{\mathbb{Z}}{n \mathbb{Z}}\\), one for each natural number \\(n\\). However, in the world of Rust we now have a problem since it is not in general possible (except in simple cases, for example by using const generics) to define infinitely many types, one for each instance of another type.

Algebraeons workaround for this problem is to use instances of certain types (types which implementing the `Structure` trait), rather than types themselves, to represent mathematical sets of sets. In the example of the integers modulo \\(n\\), one might proceed as follows

```
use algebraeon::nzq::{Natural, Integer};
use algebraeon::sets::structure::{Structure, SetStructure};

// Define a type whose instances will represent the set of integers modulo `n`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct IntegersModuloN {
    n : Natural
}

// Implement `Structure` to indicate that instances of `IntegersModuloN` represent abstract mathematical objects.
impl Structure for IntegersModuloN {}

// Implement `SetStructure` to indicate that instances of `IntegersModuloN` represent sets whose elements are represented by instances of `Integer`.
impl SetStructure for IntegersModuloN {
    type Set = Integer;
}
```


<!-- the `Structure` trait. Whenever a type implements `Structure` it means 

Compromise

In mathematics

Algebraeon has a system of traits

For example
 - Instances of types implementing `Structure` represent abstract mathematical objects which 
 - Instances of types implementing `SetStructure` represent mathematical sets.
 - Instances of types implementing `GroupStructure` represent mathematical sets with the additional structure of a group.
 - Instances of types implementing `RingStructure` represent mathematical sets with the additional structure of a ring. -->