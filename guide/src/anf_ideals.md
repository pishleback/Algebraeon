# Ideals in Algebraic Rings of Integers

## Factoring example

```rust
use algebraeon::rings::polynomial::PolynomialFromStr;
use algebraeon::{
    nzq::*,
    rings::{polynomial::Polynomial, structure::*},
};

// Construct the ring of integers Z[i]
let anf = Polynomial::<Rational>::from_str("x^2+1", "x")
    .unwrap()
    .algebraic_number_field()
    .unwrap();
let roi = anf.compute_ring_of_integers();

// The ideal (27i - 9) in Z[i]
let ideal = roi.ideals().principal_ideal(
    &roi.try_anf_to_roi(&Polynomial::from_str("27*x-9", "x").unwrap())
        .unwrap(),
);

// Factor the ideal
for (prime_ideal, power) in roi
    .ideals()
    .factorizations()
    .into_powers(roi.ideals().factor_ideal(&ideal).unwrap())
{
    println!("power = {power} prime_ideal_factor = {:?}", prime_ideal);
}

// There's not yet a nice way to print ideals so the output is messy
// But it prints the following factorization into primes
// ideal = (1+i) * (1+2i) * (3)^2
```