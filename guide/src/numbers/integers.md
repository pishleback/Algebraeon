# Naturals

## Constructing integers

There are several ways to construct integers. Some of them are:

- `Integer::ZERO`, `Integer::ONE`, `Integer::TWO`
- from a string: `Integer::from_str()`
- from builtin unsigned or signed types `Integer::from(-4i32)`
- from a `Natural`: `Integer::from(Natural::TWO)`

```rust
use std::str::FromStr;
use algebraeon::nzq::{Integer, Natural};

let one = Integer::ONE;
let two = Integer::TWO;
let n = Integer::from(42usize);
let m = Integer::from(-42);
let big = Integer::from_str("706000565581575429997696139445280900").unwrap();
let from_natural = Integer::from(Natural::from(5u32));
```

## Basic operations

Integer supports the following operators:

- `+` (addition)

- `-` (substraction or negation)

- `*` (multiplication)

- `%` (modulo)

For exponentiation, use the methods `.nat_pow(&exp)` or `.int_pow(&exp)`.

## Available functions

- `abs`

```rust
use algebraeon::nzq::{Integer, Natural};
use algebraeon::rings::structure::*;

let a = Integer::from(-12);
let b = Integer::from(5);

// Basic operations
let sum = &a + &b;                              // -7
let neg = -&a;                                  // 12
let sub = &a - &b;                              // -17
let product = &a * &b;                          // -60
let power = a.nat_pow(&Natural::from(5u32));    // -248832
let modulo = &a % &b;                           // 3

assert_eq!(sum, Integer::from(-7));
assert_eq!(neg, Integer::from(12));
assert_eq!(sub, Integer::from(-17));
assert_eq!(product, Integer::from(-60));
assert_eq!(power, Integer::from(-248832));
assert_eq!(modulo, Integer::from(3));

// abs
use algebraeon::nzq::traits::Abs;
let abs_a = a.abs();
assert_eq!(abs_a, Natural::from(12u32));
```
