# Algebraeon
Algebraeon is a computer algebra system written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations rather than approximate numerical solutions.

## What it can do
Algebraeon can currently do the following:
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
   - Trager's algorithm for factoring polynomials over algebraic number fields.
 - Expressing symmetric polynomials in terms of elementary symmetric polynomials.
 - Computations with algebraic numbers:
   - Real root isolation and arithmetic.
   - Complex root isolation and arithmetic.
 - Computations with multiplication tables for small finite groups.
 - Todd-Coxeter algorithm for the enumeration of finite index cosets of a finitely generated groups.
and much more.

## Planned algorithms
 - Fast integer factorization.
 - LLL basis reduction algorithm.
 - Universal cyclotomic field.
 - Ideals in algebraic number fields.
 - Algebraic closure and Galois theory of finite fields.
 - Splitting fields of algebraic number fields.
 - Galois groups of algebraic number fields.
 - P-adic root approximation and arithmetic.

# Usage
Add the required crates to your ``cargo.toml`` file
 - [algebraeon-sets](https://crates.io/crates/algebraeon-sets)
 - [algebraeon-groups](https://crates.io/crates/algebraeon-groups)
 - [algebraeon-rings](https://crates.io/crates/algebraeon-rings)

## Factoring a Polynomial
```
use algebraeon_rings::{
    polynomial::polynomial::*,
    ring_structure::{elements::*, structure::*},
};
use malachite_nz::integer::Integer;

let x = &Polynomial::<Integer>::var().into_ergonomic();
let f = (x.pow(15) - 1).into_verbose();
println!("f = {}", f.factor().unwrap());
/*
Output:
    f = 1 * ((-1)+λ) * (1+λ+λ^2) * (1+λ+λ^2+λ^3+λ^4) * (1+(-1)λ+λ^3+(-1)λ^4+λ^5+(-1)λ^7+λ^8)
*/
```

## Jordan Normal Form of a Matrix
```
use algebraeon_rings::{linear::matrix::*, number::algebraic::isolated_roots::complex::*};
use algebraeon_sets::structure::*;
use malachite_q::Rational;
// Construct a matrix
let a = Matrix::<Rational>::from_rows(vec![
    vec![5, 4, 2, 1],
    vec![0, 1, -1, -1],
    vec![-1, -1, 3, 0],
    vec![1, 1, -1, 2],
]);
// Put it into Jordan Normal Form
let j = MatrixStructure::new(ComplexAlgebraic::structure()).jordan_normal_form(&a);
j.pprint();
/*
Output:
    / 2    0    0    0 \
    | 0    1    0    0 |
    | 0    0    4    1 |
    \ 0    0    0    4 /
*/
```


