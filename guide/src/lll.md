# LLL Lattice Basis Reduction

Algebraeon implements the [Lenstra–Lenstra–Lovász lattice basis reduction algorithm](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm).

## Example: Approximate Algebraic Numbers

This example shows how LLL can be used to find a polynomial with a real root whose value is only known approximately, in this case, the Golden ratio \\(\varphi \approx 1.6180\\).

```rust
use algebraeon::{
    nzq::{Integer, Rational},
    rings::{
        matrix::{Matrix, StandardInnerProduct},
        polynomial::Polynomial,
    },
    sets::structure::MetaType,
};
use std::str::FromStr;

// Find a quadratic polynomial with ~φ = ~1.6180 (the Golden ratio) as a root.
let m = Matrix::<Integer>::from_rows(vec![
    vec![1, 0, 0, 10000], //  1   * 10^4
    vec![0, 1, 0, 16180], // ~φ   * 10^4
    vec![0, 0, 1, 26180], // ~φ^2 * 10^4
]);

let (u, _) = m.lll_integral_row_reduction_algorithm(
    &StandardInnerProduct::new(Integer::structure()),
    &Rational::from_str("3/4").unwrap(),
);

println!(
    "poly = {}",
    Polynomial::<Integer>::from_coeffs(u.get_row(0))
);

/*
Output:
    poly = λ^2-λ-1
*/
```