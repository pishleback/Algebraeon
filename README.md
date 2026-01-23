<p align="center"><img width="256" height="256" alt="icon" src="https://github.com/user-attachments/assets/bb2e40e2-bd57-47ea-b7d8-2d16b2e6ca5b" /></p>

# What is it?
Algebraeon is a computational algebra system (CAS) written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations over approximate numerical solutions. Algebraeon is in early stages of development and the API subject to change. Algebraeon uses [Malachite](https://www.malachite.rs/) under the hood for arbitrary sized integer and rational numbers.

 - See the [User Guide](https://pishleback.github.io/Algebraeon/) (a work in progress) to get started.
 - There is a [![Discord](https://img.shields.io/badge/Discord-Join%20Chat-7289DA?logo=discord&logoColor=white)](https://discord.gg/DBqbPqPMKR) for informal discussions about Algebraeon.
 - [Published to crates.io](https://crates.io/crates/algebraeon).
 - [Formal documentation for the most recent release](https://docs.rs/algebraeon/latest/algebraeon/).

Contributions are welcome and I am happy to do some hand-holding to help out new contributors. See below for the best ways to contribute.

# As a Python Library

Algebraeon can be used as a Python library, though it is not fully featured in that case.
 - [GitHub](https://github.com/pishleback/Algebraeon-Python)
 - [User Guide](https://pishleback.github.io/Algebraeon-Python/)

# Examples

## Factoring Integers

To factor large integers using Algebraeon

```rust
use algebraeon::nzq::Natural;
use algebraeon::rings::structure::{MetaFactoringMonoid, UniqueFactorizationMonoidSignature};
use algebraeon::sets::structure::ToStringSignature;
use std::str::FromStr;

let n = Natural::from_str("706000565581575429997696139445280900").unwrap();
let f = n.clone().factor();
println!(
    "{} = {}",
    n,
    Natural::structure_ref().factorizations().to_string(&f)
);
/*
Output:
    706000565581575429997696139445280900 = 2^2 × 5^2 × 6988699669998001 × 1010203040506070809
*/
```

Algebraeon implements [Lenstra elliptic-curve factorization](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization) for quickly finding prime factors with around 20 digits.

## Factoring Polynomials

Factor the polynomials $x^2 - 5x + 6$ and $x^{15} - 1$.

```rust
use algebraeon::rings::{parsing::parse_integer_polynomial, structure::MetaFactoringMonoid};

let f = parse_integer_polynomial("x^2 - 5*x + 6", "x").unwrap();
println!("{} = {}", f, f.factor());
/*
Output:
    λ^2-5λ+6 = 1 * (λ-2) * (λ-3)
*/

let f = parse_integer_polynomial("x^15 - 1", "x").unwrap();
println!("{} = {}", f, f.factor());
/*
Output:
    λ^15-1 = 1 * (λ-1) * (λ^2+λ+1) * (λ^4+λ^3+λ^2+λ+1) * (λ^8-λ^7+λ^5-λ^4+λ^3-λ+1)
*/
```

so

```math
x^2 - 5x + 6 = (x-2)(x-3)
```

```math
x^{15}-1 = (x-1)(x^2+x+1)(x^4+x^3+x^2+x+1)(x^8-x^7+x^5-x^4+x^3-x+1)
```

## Linear Systems of Equations

Find the general solution to the linear system

```math
a \begin{pmatrix}3 \\ 4 \\ 1\end{pmatrix} + b \begin{pmatrix}2 \\ 1 \\ 2\end{pmatrix} + c \begin{pmatrix}1 \\ 3 \\ -1\end{pmatrix} = \begin{pmatrix}5 \\ 5 \\ 3\end{pmatrix}
```

for integers $a$, $b$ and $c$.

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::linear::finitely_free_module::RingToFinitelyFreeModuleSignature;
use algebraeon::rings::matrix::Matrix;
use algebraeon::sets::structure::MetaType;

let m = Matrix::<Integer>::from_rows(vec![vec![3, 4, 1], vec![2, 1, 2], vec![1, 3, -1]]);
let y = vec![5.into(), 5.into(), 3.into()];
for x in Integer::structure()
    .free_module(3)
    .affine_subsets()
    .affine_basis(&m.row_solution_set(&y))
{
    println!("{:?}", x);
}
/*
Output:
    [Integer(0), Integer(2), Integer(1)]
    [Integer(1), Integer(1), Integer(0)]
*/
```

so two solutions are given by $(a, b, c) = (0, 2, 1)$ and $(a, b, c) = (1, 1, 0)$ and _every_ solution is a linear combination of these two solutions; The general solution is given by all $(a, b, c)$ such that

```math
\begin{pmatrix}a \\ b \\ c\end{pmatrix} = s\begin{pmatrix}0 \\ 2 \\ 1\end{pmatrix} + t\begin{pmatrix}1 \\ 1 \\ 0\end{pmatrix}
```

where $s$ and $t$ are integers such that $s + t = 1$.

## Complex Root Isolation

Find all complex roots of the polynomial
$$f(x) = x^5 + x^2 - x + 1$$

```rust
use algebraeon::rings::parsing::parse_integer_polynomial;

let f = parse_integer_polynomial("x^5 + x^2 - x + 1", "x").unwrap();
// Find the complex roots of f(x)
for root in f.all_complex_roots() {
    println!("root {} of degree {}", root, root.degree());
}

/*
Output:
    root ≈-1.328 of degree 3
    root ≈0.662-0.559i of degree 3
    root ≈0.662+0.559i of degree 3
    root -i of degree 2
    root i of degree 2
*/

```

Despite the output, the roots found are _not_ numerical approximations. Rather, they are stored internally as exact algebraic numbers by using isolating boxes in the complex plane and isolating intervals on the real line.

<img width="2158" height="1308" alt="Capture" src="https://github.com/user-attachments/assets/c727c71a-2345-45e4-81a5-35c8066024ea" />