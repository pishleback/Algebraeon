# Continued Fractions

## Computing Continued Fraction Coefficients

This example computes the first \\(10\\) continued fraction coefficients of \\(\pi\\), \\(e\\), and \\(\sqrt[3]{2}\\).

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::real_intervals::{e, pi};
use algebraeon::rings::continued_fraction::{
    MetaToSimpleContinuedFraction, SimpleContinuedFraction,
};
use algebraeon::rings::isolated_algebraic::RealAlgebraic;
use algebraeon::rings::structure::{MetaPositiveRealNthRoot, MetaRing};

println!(
    "e: {:?}",
    e().simple_continued_fraction()
        .into_iter()
        .take(10)
        .collect::<Vec<_>>()
);

println!(
    "pi: {:?}",
    pi().simple_continued_fraction()
        .into_iter()
        .take(10)
        .collect::<Vec<_>>()
);

println!(
    "cube_root(2): {:?}",
    RealAlgebraic::from_int(Integer::from(2))
        .cube_root()
        .unwrap()
        .simple_continued_fraction()
        .into_iter()
        .take(10)
        .collect::<Vec<_>>()
);

/*
Output:
    e: [Integer(2), Integer(1), Integer(2), Integer(1), Integer(1), Integer(4), Integer(1), Integer(1), Integer(6), Integer(1)]
    pi: [Integer(3), Integer(7), Integer(15), Integer(1), Integer(292), Integer(1), Integer(1), Integer(1), Integer(2), Integer(1)]
    cube_root(2): [Integer(1), Integer(3), Integer(1), Integer(5), Integer(1), Integer(1), Integer(4), Integer(1), Integer(1), Integer(8)]
*/
```

## Defining a Real Number by Its Continued Fraction Coefficients

Defining the Golden Ratio in terms of its continued fraction. 

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::real_intervals::Point;
use algebraeon::rings::continued_fraction::SimpleContinuedFraction;
use algebraeon::rings::structure::MetaRealSubset;

#[derive(Debug, Clone)]
struct MyContinuedFraction {}

impl SimpleContinuedFraction for MyContinuedFraction {
    fn next(&mut self) -> Integer {
        Integer::ONE
    }
}

let phi = Point::from_continued_fraction(MyContinuedFraction {});

println!("phi = {}", phi.as_f64());
/*
Output:
    phi = 1.618033988749895
*/
```