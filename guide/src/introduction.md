# Introduction
## What is Algebraeon?

Algebraeon is a computer algebra system written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations over approximate numerical solutions. A subset of the features are also available as a [Python module](https://pishleback.github.io/Algebraeon-Python/).

## Source Code and Contributing

The project is open source and is hosted [on GitHub here](https://github.com/pishleback/Algebraeon). See the [contributing](contributing.md) section for more information.

## Stability

The API is subject to large breaking changes at this time. I hope to stabilize things more in the not too distant future.

## Links

 - [GitHub](https://github.com/pishleback/Algebraeon)
 - [Crates.io](https://crates.io/crates/algebraeon)
 - [Algebraeon as a Python module](https://pishleback.github.io/Algebraeon-Python/)

## Crates

Algebraeon is published to crates.io under an umbrella crate [algebraeon](https://crates.io/crates/algebraeon) and re-exports the sub-crates:

- [algebraeon-sets](https://crates.io/crates/algebraeon-sets)
- [algebraeon-nzq](https://crates.io/crates/algebraeon-nzq)
- [algebraeon-groups](https://crates.io/crates/algebraeon-groups)
- [algebraeon-rings](https://crates.io/crates/algebraeon-rings)
- [algebraeon-geometry](https://crates.io/crates/algebraeon-geometry)

## Algorithms

To give a flavour of what Algebraeon can do some of the algorithms it implements are listed:

- Euclids algorithm for GCD and the extended version for obtaining Bezout coefficients.
- Polynomial GCD computations using subresultant pseudo-remainder sequences.
- AKS algorithm for natural number primality testing.
- Matrix algorithms including:
  - Putting a matrix into Hermite normal form. In particular putting it into echelon form.
  - Putting a matrix into Smith normal form.
  - Gram–Schmidt algorithm for orthogonalization and orthonormalization.
  - Putting a matrix into Jordan normal.
  - Finding the general solution to a linear or affine system of equations.
- Polynomial factoring algorithms including:
  - Kronecker's method for factoring polynomials over the integers (slow).
  - Berlekamp-Zassenhaus algorithm for factoring polynomials over the integers.
  - Berlekamp's algorithm for factoring polynomials over finite fields.
  - Cantor–Zassenhaus algorithm for factoring polynomials over finite fields.
  - Trager's algorithm for factoring polynomials over algebraic number fields.
- Expressing symmetric polynomials in terms of elementary symmetric polynomials.
- Computations with algebraic numbers:
  - Real root isolation and arithmetic.
  - Complex root isolation and arithmetic.
- Computations with multiplication tables for small finite groups.
- Todd-Coxeter algorithm for the enumeration of finite index cosets of a finitely generated groups.