# Algebraeon

Algebraeon is a computational algebra system (CAS) written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations over approximate numerical solutions. Algebraeon is in early stages of development and the API subject to change. Algebraeon uses [Malachite](https://www.malachite.rs/) under the hood for arbitrary sized integer and rational numbers.

See the [user guide](https://pishleback.github.io/Algebraeon/) (a work in progress) to get started. 

There is [![Discord](https://img.shields.io/badge/Discord-Join%20Chat-7289DA?logo=discord&logoColor=white)](https://discord.gg/DBqbPqPMKR) for informal discussions on the project.

The latest published version of Algebraeon is hosted on crates.io [here](https://crates.io/crates/algebraeon) and the formal documentation for the most recent release is [here](https://docs.rs/algebraeon/latest/algebraeon/).

# Examples

## Factoring Integers

To factor large integers using Algebraeon

```
use std::str::FromStr;
use algebraeon::{nzq::Natural, rings::rings::natural::factorization::factor};

let n = Natural::from_str("706000565581575429997696139445280900").unwrap();
let f = factor(n.clone()).unwrap();
println!("{} = {}", n, f);
/*
Output:
    706000565581575429997696139445280900 = 2^2 × 5^2 × 6988699669998001 × 1010203040506070809
*/
```

Algebraeon implements [Lenstra elliptic-curve factorization](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization) for quickly finding prime factors with around 20 digits.

## Factoring Polynomials

Factor the polynomials $x^2 - 5x + 6$ and $x^{15} - 1$.

```
use algebraeon::rings::{polynomial::*, structure::*};
use algebraeon::nzq::Integer;

let x = &Polynomial::<Integer>::var().into_ergonomic();
let f = (x.pow(2) - 5*x + 6).into_verbose();
println!("f(λ) = {}", f.factor().unwrap());
/*
Output:
    f(λ) = 1 * ((-2)+λ) * ((-3)+λ)
*/

let f = (x.pow(15) - 1).into_verbose();
println!("f(λ) = {}", f.factor().unwrap());
/*
Output:
    f(λ) = 1 * ((-1)+λ) * (1+λ+λ^2) * (1+λ+λ^2+λ^3+λ^4) * (1+(-1)λ+λ^3+(-1)λ^4+λ^5+(-1)λ^7+λ^8)
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

```
use algebraeon::nzq::Integer;
use algebraeon::rings::linear::finitely_free_modules::RingToFinitelyFreeModuleStructure;
use algebraeon::rings::linear::matrix::Matrix;
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

```
use algebraeon::rings::{polynomial::*, structure::*};
use algebraeon::nzq::Integer;

let x = &Polynomial::<Integer>::var().into_ergonomic();
let f = (x.pow(5) + x.pow(2) - x + 1).into_verbose();
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

Despite the output, the roots found are _not_ numerical approximations. Rather, they are stored internally as exact algebraic numbers by using isolating boxes in the complex plane.

# Contributing

Contributions are welcome and I am happy to do some hand-holding to help out new contributors. See below for the best ways to get started with contributions.

## Using the issue tracker

Use the issue tracker if you have questions, feature requests, bugs reports, etc.

## What to work on?

If you're unsure what you could do to help:
 - Have a look through the issue tracker for unassigned issues.
 - Ask in the discord server.

## Changing the code-base

Fork this repository, make changes, and submit a pull request.

## Getting started

Algebraeon is organized as a [cargo workspace](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html). Run `cargo test` in the root directory to build and run all tests.

A suggested workflow for getting started:
- Checkout the repository.
- Create a new binary in `examples/src/bin`, for example `my_main.rs`.
- Copy an example into `my_main.rs`.
- Run it with `cargo run --bin my_main` in the root directory.

## CLA

If you wish to contribute code we ask that you sign a CLA. The CLA allows us to relicense future versions of Algebraeon, including your contribution, under any licences the Free Software Foundation classifies as a Free Software Licence and which are approved by the Open Source Initiative as Open Source licences. It does not allow us to relicense of your contribution under any other more restrictive licences. The CLA can be signed on GitHub once you submit a pull request.
