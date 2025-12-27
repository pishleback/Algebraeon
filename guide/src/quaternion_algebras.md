# Quaternion Algebras

Quaternion algebras can be created over any field: the rationals, number fields, finite fields, ...

The main constructor for the quaternion algebra over F, such that \\(i^2 = a\\) and \\(j^2 = b\\) is `QuaternionAlgebraStructure::new(F, a, b)`.

```rust
use algebraeon::nzq::Rational;
use algebraeon::rings::quaternion_algebra::QuaternionAlgebraStructure;
use algebraeon::rings::structure::*;
use algebraeon::sets::structure::{EqSignature, MetaType};

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
