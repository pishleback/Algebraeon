# Rational Numbers

The `Rational` type represents a number of the form `a/b`, where `a` is an integer and `b` is a non-zero integer. It is a wrapper around `malachite_q::Rational`.

## Constructing rationals

There are several ways to construct rational numbers. Some of them are:

* `Rational::ZERO`, `Rational::ONE`, `Rational::TWO`, `Rational::ONE_HALF`
* from a string: `Rational::from_str()`
* from builtin signed/unsigned types: `Rational::from(4u32)`, `Rational::from(-3i64)`
* from `Integer` or `Natural`: `Rational::from(Integer::from(-5))`
* from two integers: `Rational::from_integers(n, d)`

## Available functions

The following methods are available on `Rational`:

- `abs()` – returns the absolute value
- `floor()` – returns the greatest integer less than or equal to the rational
- `ceil()` – returns the smallest integer greater than or equal to the rational
- `approximate(max_denominator)` – approximates the rational with another having a bounded denominator
- `simplest_rational_in_closed_interval(a, b)` – finds the simplest rational between two bounds
- `simplest_rational_in_open_interval(a, b)` – finds the simplest rational strictly between two bounds
- `decimal_string_approx()` – returns a decimal string approximation of the rational
- `exhaustive_rationals()` – returns an infinite iterator over all reduced rational numbers
- `into_abs_numerator_and_denominator()` – returns the absolute numerator and denominator as `Natural`s
- `try_from_float_simplest(x: f64)` – converts a float into the simplest rational representation

```rust
use std::str::FromStr;
use algebraeon::nzq::{Rational, Integer, Natural};

let zero = Rational::ZERO;
let one = Rational::ONE;
let half = Rational::ONE_HALF;

let r1 = Rational::from(5u32);
let r2 = Rational::from(-3i64);
let r3 = Rational::from(Integer::from(-7));
let r4 = Rational::from(Natural::from(10u32));
let r5 = Rational::from_str("42/7").unwrap();
let r6 = Rational::from_integers(3, 4);

let a = Rational::from_str("2/3").unwrap();
let b = Rational::from_str("1/6").unwrap();

// Basic operations
let sum = &a + &b;        // 5/6
let diff = &a - &b;       // 1/2
let product = &a * &b;    // 1/9
let quotient = &a / &b;   // 4
let negated = -&a;        // -2/3

// Available functions
let r = Rational::from_str("-2/5").unwrap();

// abs
use algebraeon::nzq::traits::Abs;
assert_eq!(r.clone().abs(), Rational::from_str("2/5").unwrap());

// ceil
use algebraeon::nzq::traits::Ceil;
assert_eq!(r.clone().ceil(), Integer::from(0));

// floor
use algebraeon::nzq::traits::Floor;
assert_eq!(r.clone().floor(), Integer::from(-1));

// into_abs_numerator_and_denominator
use algebraeon::nzq::traits::Fraction;
let (n, d) = r.into_abs_numerator_and_denominator();
assert_eq!(n, Natural::from(2u32));
assert_eq!(d, Natural::from(5u32));

// approximate
let approx = Rational::from_str("355/113").unwrap()
    .approximate(&Natural::from(10u32));
assert_eq!(approx, Rational::from_str("22/7").unwrap());

// simplest_rational_in_closed_interval
let a = Rational::from_str("2/5").unwrap();
let b = Rational::from_str("3/5").unwrap();
let simple = Rational::simplest_rational_in_closed_interval(&a, &b);
assert_eq!(simple, Rational::from_str("1/2").unwrap());

// simplest_rational_in_open_interval
let open_a = Rational::from_str("1/3").unwrap();
let open_b = Rational::from_str("2/3").unwrap();
let open_simple = Rational::simplest_rational_in_open_interval(&open_a, &open_b);
assert!(open_simple > open_a && open_simple < open_b);

// try_from_float_simplest
let from_float = Rational::try_from_float_simplest(0.5).unwrap();
assert_eq!(from_float, Rational::ONE_HALF);

// decimal_string_approx
let dec = simple.decimal_string_approx();
assert_eq!(dec, "0.5");

// Iteration
let mut iter = Rational::exhaustive_rationals();
assert_eq!(iter.next(), Some(Rational::ZERO));
assert_eq!(iter.next(), Some(Rational::ONE));
```

