# Algebraeon
Algebraeon is a computer algebra system written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations over approximate numerical solutions. Algebraeon is in early stages of development and the API is currently highly unstable and subject to change. Algebraeon uses [Malachite](https://www.malachite.rs/) for arbitrary sized integer and rational numbers.

# Usage
Add
```ignore
[dependencies]
algebraeon = "0.1.3"
```
to your `cargo.toml` to make Algebraeon available. Copy an example below to get started.

## Factoring Integers
To factor large integers using Algebraeon

```
use std::str::FromStr;
use algebraeon::{nzq::natural::Natural, rings::number::natural::factorization::factor};

let n = Natural::from_str("706000565581575429997696139445280900").unwrap();
let f = factor(n.clone()).unwrap();
println!("{} = {}", n, f);
/*
Output:
    706000565581575429997696139445280900 = 2^2 × 5^2 × 6988699669998001 × 1010203040506070809
*/
```

Algebraeon implements [Lenstra elliptic-curve factorization](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization), the third-fastest known factoring algorithm known, for quickly finding prime factors with around 20 digits.


## Factoring Polynomials
Factor the polynomials $x^2 - 5x + 6$ and $x^{15} - 1$.
```
use algebraeon::rings::{
    polynomial::polynomial::*,
    structure::{elements::*, structure::*},
};
use algebraeon::nzq::integer::*;

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
use algebraeon::rings::linear::matrix::Matrix;
use algebraeon::nzq::integer::*;
let x = Matrix::<Integer>::from_rows(vec![vec![3, 4, 1], vec![2, 1, 2], vec![1, 3, -1]]);
let y = Matrix::<Integer>::from_rows(vec![vec![5, 5, 3]]);
let s = x.row_solution_lattice(y);
s.pprint();
/*
Output:
    Start Affine Lattice
    Offset
    ( 2    0    -1 )
    Start Linear Lattice
    ( 1    -1    -1 )
    End Linear Lattice
    End Affine Lattice
*/
```
so the general solution is all $a$, $b$, $c$ such that
```math
\begin{pmatrix}a \\ b \\ c\end{pmatrix} = \begin{pmatrix}2 \\ 0 \\ -1\end{pmatrix} + t\begin{pmatrix}1 \\ -1 \\ -1\end{pmatrix}
```
for some integer $t$.

## Complex Root Isolation
Find all complex roots of the polynomial
$$f(x) = x^5 + x^2 - x + 1$$
```
use algebraeon::rings::{polynomial::polynomial::*, structure::elements::*};
use algebraeon::nzq::integer::*;

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

## Factoring Multivariable Polynomials
Factor the following multivariable polynomial with integer coefficients
```math
f(x, y) = 6x^4 - 6x^3y^2 + 6xy - 6x - 6y^3 + 6y^2
```
```
use algebraeon::{nzq::integer::*, rings::{polynomial::multipoly::*, structure::{elements::*, structure::*}}};

let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();

let f = (6 * (x.pow(4) - x.pow(3) * y.pow(2) + x * y - x - y.pow(3) + y.pow(2))).into_verbose();
println!("f(x, y) = {}", f.factor().unwrap());

/*
Output:
    f(x, y) = 1 * ((3)1) * ((2)1) * (x+(-1)y^2) * (x^3+y+(-1)1)
*/
```
so the factorization of $f(x, y)$ is
```math
f(x, y) = 2 \times 3 \times (x^3 + y - 1) \times (y^2 - x)
```


## P-adic Root Finding
Find the $2$-adic square roots of $17$.
```
use algebraeon::nzq::natural::*;
use algebraeon::nzq::integer::*;
use algebraeon::rings::{polynomial::polynomial::*, structure::elements::*};
let x = Polynomial::<Integer>::var().into_ergonomic();
let f = (x.pow(2) - 17).into_verbose();
for mut root in f.all_padic_roots(&Natural::from(2u32)) {
    println!("{}", root.truncate(&20.into()).string_repr()); // Show 20 2-adic digits
}
/*
Output:
    ...00110010011011101001
    ...11001101100100010111
*/
```
Truncating to the last 16 bits it can be verified that, modulo $2^{16}$, the square of these values is $17$.
```
let a = 0b0010011011101001u16;
assert_eq!(a.wrapping_mul(a), 17u16);
let b = 0b1101100100010111u16;
assert_eq!(b.wrapping_mul(b), 17u16);
```

## Enumerating a Finitely Generated Group
Let $G$ be the finitely generated group generated by $3$ generators $a, b, c$ subject to the relations $a^2 = b^2 = c^2 = (ab)^3 = (bc)^5 = (ac)^2 = e$.
```math
G = \langle a, b, c : a^2 = b^2 = c^2 = (ab)^3 = (bc)^5 = (ac)^2 = e \rangle
```
Using Algebraeon, $G$ is found to be a finite group of order $120$.
```
use algebraeon::groups::free_group::todd_coxeter::*;
let mut g = FinitelyGeneratedGroupPresentation::new();
// Add the 3 generators
let a = g.add_generator();
let b = g.add_generator();
let c = g.add_generator();
// Add the relations
g.add_relation(a.pow(2));
g.add_relation(b.pow(2));
g.add_relation(c.pow(2));
g.add_relation((&a * &b).pow(3));
g.add_relation((&b * &c).pow(5));
g.add_relation((&a * &c).pow(2));
// Count elements
let (n, _) = g.enumerate_elements();
assert_eq!(n, 120);
```

## Jordan Normal Form of a Matrix
```
use algebraeon::nzq::rational::*;
use algebraeon::rings::{linear::matrix::*, number::algebraic::complex::*};
use algebraeon::sets::structure::*;
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

## Computing Discriminants
Algebraeon can find an expression for the discriminant of a polynomial in terms of the polynomials coefficients.
```
use algebraeon::rings::polynomial::{
    multipoly::{MultiPolynomial, Variable},
    polynomial::Polynomial,
};
use algebraeon::nzq::integer::*;

let a_var = Variable::new("a");
let b_var = Variable::new("b");
let c_var = Variable::new("c");
let d_var = Variable::new("d");
let e_var = Variable::new("e");

let a = MultiPolynomial::<Integer>::var(a_var);
let b = MultiPolynomial::<Integer>::var(b_var);
let c = MultiPolynomial::<Integer>::var(c_var);
let d = MultiPolynomial::<Integer>::var(d_var);
let e = MultiPolynomial::<Integer>::var(e_var);

let p =
    Polynomial::<MultiPolynomial<Integer>>::from_coeffs(vec![c.clone(), b.clone(), a.clone()]);
println!("p(λ) = {}", p);
println!("disc(p) = {}", p.discriminant().unwrap());

println!();

let p = Polynomial::<MultiPolynomial<Integer>>::from_coeffs(vec![
    d.clone(),
    c.clone(),
    b.clone(),
    a.clone(),
]);
println!("p(λ) = {}", p);
println!("disc(p) = {}", p.discriminant().unwrap());

println!();

let p = Polynomial::<MultiPolynomial<Integer>>::from_coeffs(vec![
    e.clone(),
    d.clone(),
    c.clone(),
    b.clone(),
    a.clone(),
]);
println!("p(λ) = {}", p);
println!("disc(p) = {}", p.discriminant().unwrap());

/*
Output:
    p(λ) = (c)+(b)λ+(a)λ^2
    disc(p) = (-4)ac+b^2

    p(λ) = (d)+(c)λ+(b)λ^2+(a)λ^3
    disc(p) = (-27)a^2d^2+(18)abcd+(-4)ac^3+(-4)b^3d+b^2c^2
*/
```
so
```math
\mathop{\text{disc}}(ax^2 + bx + c) = b^2 - 4ac
```

```math
\mathop{\text{disc}}(ax^3 + bx^2 + cx + d) = b^2c^2 - 4ac^3 - 4b^3d - 27a^2d^2 + 18abcd
```


# Crates
Algebraeon is published to crates.io under an umbrella crate [algebraeon](https://crates.io/crates/algebraeon) made up of:
 - [algebraeon-sets](https://crates.io/crates/algebraeon-sets)
 - [algebraeon-nzq](https://crates.io/crates/algebraeon-nzq)
 - [algebraeon-groups](https://crates.io/crates/algebraeon-groups)
 - [algebraeon-rings](https://crates.io/crates/algebraeon-rings)
 - [algebraeon-geometry](https://crates.io/crates/algebraeon-geometry)

# Algorithms
Algebraeon currently implements the following algorithms:
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

# Getting Help
If you have questions, concerns, bug reports, etc, please file an issue in this repository's Issue Tracker.

# Contributing
Contributions are welcome. There are two primary ways to contribute:

## Using the issue tracker
Use the issue tracker to suggest feature requests, report bugs, and ask questions.

## Changing the code-base
You should fork this repository, make changes in your own fork, and then submit a pull request. New code should have associated unit tests that validate implemented features and the presence or lack of defects.

Algebraeon is organized as a [cargo workspace](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html). Run `cargo test` in the root directory to build and run all tests.

A suggested workflow for testing new features:
 - Create a new binary in `examples/src/bin`, for example `my_main.rs`.
 - To run, use `cargo run --bin my_main` in the root directory.
 - Test any changes to the codebase with unit tests and/or using `my_main.rs`.

## CLA
Anyone who contributes code is required to sign the CLA. You can sign the CLA when you submit a pull request. The CLA allows for us to relicense future versions of Algebraeon, including your contribution, under any licences the Free Software Foundation classifies as Free Software Licence and which are approved by the Open Source Initiative as Open Source licences. It does not allow us to relicense of your contribution under any other more restrictive licences.
