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

## Stability
Algebraeon is still in an early experimental stage of development and as such the API will be unstable.