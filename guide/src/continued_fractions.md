# Continued Fractions

## Computing Continued Fraction Coefficients

This example computes the first \\(10\\) continued fraction coefficients of \\(\pi\\), \\(e\\), and \\(\sqrt[3]{2}\\).

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::{e, pi};
use algebraeon::rings::continued_fraction::{
    MetaToSimpleContinuedFractionSignature, SimpleContinuedFraction,
};
use algebraeon::rings::isolated_algebraic::RealAlgebraic;
use algebraeon::rings::structure::{MetaPositiveRealNthRootSignature, MetaRingSignature};

println!(
    "e: {:?}",
    e().simple_continued_fraction()
        .iter()
        .take(10)
        .collect::<Vec<_>>()
);

println!(
    "pi: {:?}",
    pi().simple_continued_fraction()
        .iter()
        .take(10)
        .collect::<Vec<_>>()
);

println!(
    "cube_root(2): {:?}",
    RealAlgebraic::from_int(Integer::from(2))
        .cube_root()
        .unwrap()
        .simple_continued_fraction()
        .iter()
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

We can also compute continued fractions for rational numbers, for example

\\[-\frac{-5678}{1234} = -5 + \cfrac{1}{2 + \cfrac{1}{1 + \cfrac{1}{1 + \cfrac{1}{30 + \cfrac{1}{4}}}}}\\]

```rust
use algebraeon::{
    nzq::{Integer, Rational},
    rings::continued_fraction::{MetaToSimpleContinuedFractionSignature, SimpleContinuedFraction},
};
use std::str::FromStr;

assert_eq!(
    Rational::from_str("-5678/1234")
        .unwrap()
        .simple_continued_fraction()
        .iter()
        .map(|x| x.as_ref().clone())
        .collect::<Vec<_>>(),
    vec![
        Integer::from(-5),
        Integer::from(2),
        Integer::from(1),
        Integer::from(1),
        Integer::from(30),
        Integer::from(4)
    ]
);
```


## Defining a Real Number by Its Continued Fraction Coefficients

Defining the real number

\\[1 + \cfrac{1}{2 + \cfrac{1}{3 + \cfrac{1}{4 + \ddots}}} = 1.433127...\\]

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::RealApproximatePoint;
use algebraeon::rings::continued_fraction::IrrationalSimpleContinuedFractionGenerator;
use algebraeon::rings::structure::MetaRealSubsetSignature;

#[derive(Debug, Clone)]
struct MyContinuedFraction {
    n: usize,
}

impl IrrationalSimpleContinuedFractionGenerator for MyContinuedFraction {
    fn next(&mut self) -> Integer {
        self.n += 1;
        Integer::from(self.n)
    }
}

let value =
    RealApproximatePoint::from_continued_fraction(MyContinuedFraction { n: 0 }.into_continued_fraction());

println!("value = {}", value.as_f64());

/*
Output:
    value = 1.4331274267223117
*/
```

Defining the Golden Ratio

\\[\varphi = 1 + \cfrac{1}{1 + \cfrac{1}{1 + \cfrac{1}{1 + \ddots}}}\\]

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::RealApproximatePoint;
use algebraeon::rings::continued_fraction::PeriodicSimpleContinuedFraction;
use algebraeon::rings::structure::MetaRealSubsetSignature;

let phi = RealApproximatePoint::from_continued_fraction(
    PeriodicSimpleContinuedFraction::new(vec![], vec![Integer::from(1)]).unwrap(),
);

println!("phi = {}", phi.as_f64());

/*
Output:
    phi = 1.618033988749895
*/

```

Defining \\(\sqrt{2}\\)

\\[\sqrt{2} = 1 + \cfrac{1}{2 + \cfrac{1}{2 + \cfrac{1}{2 + \ddots}}}\\]

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::approximation::RealApproximatePoint;
use algebraeon::rings::continued_fraction::PeriodicSimpleContinuedFraction;
use algebraeon::rings::structure::MetaRealSubsetSignature;

let sqrt2 = RealApproximatePoint::from_continued_fraction(
    PeriodicSimpleContinuedFraction::new(vec![Integer::from(1)], vec![Integer::from(2)])
        .unwrap(),
);

println!("sqrt2 = {}", sqrt2.as_f64());

/*
Output:
    sqrt2 = 1.4142135623730951
*/
```