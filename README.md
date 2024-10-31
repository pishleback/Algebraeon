# Rust-Computational-Algebra
Rust implementation of some computational algebra algorithms. 

This project began as a faster and more maintainable replacement for an older attempt of mine written in Python: https://github.com/pishleback/Python-Computational-Algebra

## Sets
 - k-element subsets of an n-element set
 - Partitions of natural numbers

## Groups
 - Permutations
 - Groups represented by multiplication tables
   - Conjugacy classes
   - Quotients
   - Tests for isomorphism
 - Todd-Coxeter algorithm coset enumeration

## Rings
### Rings
 - Euclids algorithm for gcd
 - Extended Euclids algorithm for obtaining Bezout coefficients

### Linear Algebra
 - Hermite normal form of a matrix over a PID
 - Smith normal form of a matrix over a PID
 - General solution to a linear or affine system of equations over a PID
 - Gram–Schmidt algorithm
 - Jordan normal form of a matrix

### Polynomials
 - Kronecker's method for factoring polynomials over the integers (very slow)
 - Zassenhaus algorithm for factoring polynomials over the integers
 - Berlekamp's algorithm for factoring polynomials over finite fields
 - Factor polynomials over algebraic number fields
 - Express symmetric polynomials in terms of elementary symmetric polynomial

### Algebraic Numbers
 - Real root isolation and arithmetic
 - Complex root isolation and arithmetic

## Combinatorics
 - Binary Golay Codes
   - Find isomorphisms
 - Ternary Golay Codes

## Geometry
 - Geometric simplicial complexes
   - Intersection
   - Union
   - Difference
   - Simplification

## Planned Features
 - Fast integer factorization
 - LLL basis reduction algorithm
 - Universal cyclotomic field
 - Ideals in algebraic number fields
 - Algebraic closure and Galois theory of finite fields
 - Splitting fields of algebraic number fields
 - Galois groups of algebraic number fields
