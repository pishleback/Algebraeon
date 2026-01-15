# Rational Bounds on Real Values

Compute the following rational bounds on \\(\pi + e\\), where the difference between the upper and lower bound is at most \\(\frac{1}{100}\\).

\\[\frac{163994429}{28005120} < \pi + e < \frac{492403663}{84015360}\\]

```rust
use algebraeon::nzq::Rational;
use algebraeon::rings::approximation::{RealApproximatePoint, e, pi};
use algebraeon::rings::structure::MetaAdditionSignature;
use std::str::FromStr;

let p = RealApproximatePoint::add(&pi(), &e());

p.lock()
    .refine_to_length(&Rational::from_str("1/100").unwrap());

println!("{:?}", p.lock().rational_interval_neighbourhood());

/*
Output:
    Interval(RationalInterval { a: Rational(163994429/28005120), b: Rational(492403663/84015360) })
*/
```