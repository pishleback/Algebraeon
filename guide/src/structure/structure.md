# Structure Types

## Motivation

In mathematics there are many instances of sets of sets with some additional structure (think the set of all groups) and often these sets of sets are parameterized by other sets. For example the set of integers modulo \\(n\\) is defined for every natural number \\(n \in \mathbb{N}\\) and has the structure of a commutative ring. A mathematician would like to think that there exists infinitely many types, \\(\frac{\mathbb{Z}}{n \mathbb{Z}}\\), one for each natural number \\(n\\). However, in the world of Rust we now have a problem since it is not in general possible (except in simple cases, for example by using generics) to define infinitely many types, one for each instance of some mathematical set. 

The solution Algebraeon takes for this problem is structure types - types whose objects represent sets with additional structure. Algebraeon provides many structure traits for types whose objects represent mathematical sets with additional structure. A non-exhaustive list of structure traits is as follows:
 - `Structure`: This is the base trait for structure types. All structure types implement `Structure`.
 - `SetStructure: Structure`: Objects are sets whose elements are objects of the associated type `Set`.
 - `EqStructure: SetStructure`: It is possible to compare elements for equality.
 - `FiniteSetStructure : SetStructure`: Objects are finite sets. There is a way to obtain a `Vec<Set>` containing all the elements.
 - `RingStructure: SetStructure`: Objects are commutative rings and there are functions for doing ring operations (addition, subtraction, multiplication) with elements.
 - `FieldStructure: RingStructure`: Objects are fields and there are functions for doing field operations (addition, subtraction, multiplication, division) with elements.

## Example: Ring Structure on the Integers Modulo \\(n\\)

An example of using Algebraeon to define the ring structure on the set of integers modulo \\(n \in \mathbb{N}\\) is below.

```rust
use algebraeon::nzq::traits::AbsDiff;
use algebraeon::nzq::{Integer, Natural};

// Define a type whose instances will represent the set of integers modulo `n`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct IntegersModuloN {
    n: Natural,
}

use algebraeon::sets::structure::{EqStructure, SetStructure, Structure};

// Implement `Structure` to indicate that instances of `IntegersModuloN` represent abstract mathematical objects.
impl Structure for IntegersModuloN {}

// Implement `SetStructure` to indicate that instances of `IntegersModuloN` represent sets whose elements are represented by instances of `Integer`.
impl SetStructure for IntegersModuloN {
    type Set = Integer;

    fn is_element(&self, x : &Integer) -> bool {
        x < &self.n
    }
}

// Implement `EqStructure` so that integers can be compared for equality modulo `n`.
impl EqStructure for IntegersModuloN {
    fn equal(&self, a: &Integer, b: &Integer) -> bool {
        a.abs_diff(b) % &self.n == Natural::ZERO
    }
}

let mod_6 = IntegersModuloN {n : 6u32.into()};
assert!(mod_6.equal(&2.into(), &8.into()));
assert!(!mod_6.equal(&3.into(), &5.into()));

use algebraeon::rings::structure::{RingStructure, SemiRingStructure};

// Implement `SemiRingStructure` and `RingStructure` to equip the integers modulo `n` with the quotient ring structure.
impl SemiRingStructure for IntegersModuloN {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }

    fn one(&self) -> Self::Set {
        (Integer::ONE % &self.n).into()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        ((a + b) % &self.n).into()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        ((a * b) % &self.n).into()
    }
}
impl RingStructure for IntegersModuloN {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        (-a % &self.n).into()
    }
}

// Now `mod_6` now has the structure of a ring so Algebraeon implements the repeated squaring algorithm for taking very large powers modulo `n`.
assert!(mod_6.equal(
    &mod_6.nat_pow(&2.into(), &1000000000000u64.into()),
    &4.into()
));
```