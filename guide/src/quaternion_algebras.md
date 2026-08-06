# Quaternion Algebras

Quaternion algebras can be created over any field: the rationals, number fields, finite fields, ...

The main constructor for the quaternion algebra over F, such that \\(i^2 = a\\) and \\(j^2 = b\\) is `QuaternionAlgebraStructure::new(F, a, b)`.

```rust
use algebraeon::structures::Rational;
use algebraeon::rings::quaternion_algebra::QuaternionAlgebraStructure;
use algebraeon::rings::structure::*;
use algebraeon::structures::{EqSignature, MetaType};

let h = QuaternionAlgebraStructure::new(
    Rational::structure(),
    -Rational::ONE,
    -Rational::TWO,
);

let i = h.i();
let j = h.j();
let k = h.mul(&i, &j);

// ij = k and ji = -k
assert!(h.equal(&k, &h.k()));
assert!(h.equal(&h.mul(&j, &i), &h.neg(&k)));
```

Elements can also be written as strings and read with `parse_quaternion`. The factors of a product are multiplied in the order they are written, so `"i*j"` and `"j*i"` give different elements.

```rust
use algebraeon::structures::Rational;
use algebraeon::rings::parsing::parse_quaternion;
use algebraeon::rings::quaternion_algebra::QuaternionAlgebraStructure;
use algebraeon::rings::structure::*;
use algebraeon::structures::{EqSignature, MetaType};

let h = QuaternionAlgebraStructure::new(
    Rational::structure(),
    -Rational::ONE,
    -Rational::TWO,
);

let q = parse_quaternion("2 + 3i + 5j - 2k", &h).unwrap();
assert!(h.equal(&q, &h.from_components(
    Rational::from(2),
    Rational::from(3),
    Rational::from(5),
    Rational::from(-2),
)));

assert!(h.equal(&parse_quaternion("i*j", &h).unwrap(), &h.k()));
assert!(h.equal(&parse_quaternion("j*i", &h).unwrap(), &h.neg(&h.k())));
```
