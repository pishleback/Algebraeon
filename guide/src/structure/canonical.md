# Canonical Structures

Sometimes the situation is simple and we only want to define one set with structure rather than a family of sets, for example, the set of all rational numbers. Since sets with structure are represented in Algebraeon objects of structure tyes we will need a structure type with exactly once instance. This can be done explicitly like so

```rust
use algebraeon::{nzq::Rational, rings::structure::*, sets::structure::*};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MyRational {
    value: Rational,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MyRationalCanonicalStructure {}

impl Signature for MyRationalCanonicalStructure {}

impl SetSignature for MyRationalCanonicalStructure {
    type Set = MyRational;

    fn is_element(&self, _x: &Self::Set) -> bool {
        true
    }
}

impl EqSignature for MyRationalCanonicalStructure {
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        x == y
    }
}
```

However, Algebraeon provides a derive macro `CanonicalStructure` which reduces the boilerplate above to

```rust
use algebraeon::{nzq::Rational, rings::structure::*, sets::structure::*};

#[derive(Debug, Clone, PartialEq, Eq, CanonicalStructure)]
pub struct MyRational {
    value: Rational,
}
```

In any case, once we have the structure type `MyRationalCanonicalStructure` implementing `Signature + SetSignature<Set = MyRational>` we can go on to implement more structure traits like  `RingSignature` and `FieldSignature` for `MyRationalCanonicalStructure` to give the set of instances of `MyRational` the structure of a ring or a field.

<!-- ```
// The CanonicalStructure derive macro defines a new type MyRationalCanonicalStructure with one value and implements `Structure`, `SetStructure` and `EqStructure` for it.
// We can proceed to implement more interesting structures.

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
``` -->

<!-- # Sets

# Rings and Fields

# Groups -->