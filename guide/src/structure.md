# The Structure System

## The Structure Trait

In mathematics there are many instances of sets of sets with some additional structure (think the set of all groups) and often these sets of sets are parameterized by other sets. For example the set of integers modulo \\(n\\) is defined for every natural number \\(n \in \mathbb{N}\\) and has the structure of a commutative ring. A mathematician would like to think that there exists infinitely many types, \\(\frac{\mathbb{Z}}{n \mathbb{Z}}\\), one for each natural number \\(n\\). However, in the world of Rust we now have a problem since it is not in general possible (except in simple cases, for example by using generics) to define infinitely many types, one for each instance of some mathematical set.

The workaround Algebraeon takes is to use instances of certain types (types which implement the `Structure` trait), rather than types themselves, to represent mathematical sets of sets with structure. In the case of integers modulo \\(n\\) the approach looks like this:

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

// Now `mod_6` has the structure of a ring, Algebraeon implements the repeated squaring algorithm for taking very large powers modulo `n`.
assert!(mod_6.equal(
    &mod_6.nat_pow(&2.into(), &1000000000000u64.into()),
    &4.into()
));
```

## Canonical Structures

Sometimes we don't want to implement a whole family of sets with structure, such as the ring of integers modulo \\(n\\) for all naturals \\(n\\), rather we only want to implement a single type. An example of this scenario is the field of rational numbers. Since sets with structure are supposed to be represented by instances of a type implementing `Structure`, we'd need to define a singleton type, say `RationalCanonicalStructure`, which implements `SetStructure<Set = Rational>`. Algebraeon provides a derive macro called `CanonicalStructure` to streamline this.

```rust
use algebraeon::{nzq::Rational, rings::structure::*, sets::structure::*};

#[derive(Debug, Clone, PartialEq, Eq, CanonicalStructure)]
pub struct MyRational {
    value: Rational,
}

impl SemiRingStructure for MyRationalCanonicalStructure {
    fn zero(&self) -> Self::Set {
        MyRational {
            value: Rational::ZERO,
        }
    }

    fn one(&self) -> Self::Set {
        MyRational {
            value: Rational::ONE,
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        MyRational {
            value: &a.value + &b.value,
        }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        MyRational {
            value: &a.value * &b.value,
        }
    }
}

impl RingStructure for MyRationalCanonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        MyRational { value: -&a.value }
    }
}

// Algebraeon redefines functions defined on `MyRationalCanonicalStructure` to functions defined on `MyRational`, for example `nat_pow`, so that we can use `.nat_pow(..)` on an instance of `MyRational` without going through `MyRationalCanonicalStructure`
use std::str::FromStr;
let a = MyRational {
    value: Rational::from_str("2/3").unwrap(),
};
assert_eq!(
    a.nat_pow(&3u32.into()),
    MyRational {
        value: Rational::from_str("8/27").unwrap()
    }
);
```