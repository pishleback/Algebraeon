# Quickstart Examples

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

Factor the polynomials \\(x^2 - 5x + 6\\) and \\(x^{15} - 1\\).

```rust
use algebraeon::rings::{polynomial::*, structure::*};
use algebraeon::nzq::Integer;

let x = &Polynomial::<Integer>::var().into_ergonomic();
let f = (x.pow(2) - 5*x + 6).into_verbose();
println!("f(λ) = {}", f.factor());
/*
Output:
    f(λ) = 1 * ((-2)+λ) * ((-3)+λ)
*/

let f = (x.pow(15) - 1).into_verbose();
println!("f(λ) = {}", f.factor());
/*
Output:
    f(λ) = 1 * ((-1)+λ) * (1+λ+λ^2) * (1+λ+λ^2+λ^3+λ^4) * (1+(-1)λ+λ^3+(-1)λ^4+λ^5+(-1)λ^7+λ^8)
*/
```

so

\\[x^2 - 5x + 6 = (x-2)(x-3)\\]

\\[x^{15}-1 = (x-1)(x^2+x+1)(x^4+x^3+x^2+x+1)(x^8-x^7+x^5-x^4+x^3-x+1)\\]

## Linear Systems of Equations

Find the general solution to the linear system

\\[a \begin{pmatrix}3 \\\\ 4 \\\\ 1\end{pmatrix} + b \begin{pmatrix}2 \\\\ 1 \\\\ 2\end{pmatrix} + c \begin{pmatrix}1 \\\\ 3 \\\\ -1\end{pmatrix} = \begin{pmatrix}5 \\\\ 5 \\\\ 3\end{pmatrix}\\]

for integers \\(a\\), \\(b\\) and \\(c\\).

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

so two solutions are given by \\((a, b, c) = (0, 2, 1)\\) and \\((a, b, c) = (1, 1, 0)\\) and _every_ solution is a linear combination of these two solutions; The general solution is given by all \\((a, b, c)\\) such that

\\[\begin{pmatrix}a \\\\ b \\\\ c\end{pmatrix} = s\begin{pmatrix}0 \\\\ 2 \\\\ 1\end{pmatrix} + t\begin{pmatrix}1 \\\\ 1 \\\\ 0\end{pmatrix}\\]

where \\(s\\) and \\(t\\) are integers such that \\(s + t = 1\\).

## Complex Root Isolation

Find all complex roots of the polynomial
\\[f(x) = x^5 + x^2 - x + 1\\]

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

Despite the output, the roots found are _not_ numerical approximations. Rather, they are stored internally as exact algebraic numbers by using isolating boxes in the complex plane.

## Factoring Multivariable Polynomials

Factor the following multivariable polynomial with integer coefficients

\\[f(x, y) = 6x^4 - 6x^3y^2 + 6xy - 6x - 6y^3 + 6y^2\\]

```rust
use algebraeon::{nzq::Integer, rings::{polynomial::*, structure::*}};

let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();

let f = (6 * (x.pow(4) - x.pow(3) * y.pow(2) + x * y - x - y.pow(3) + y.pow(2))).into_verbose();
println!("f(x, y) = {}", f.factor());

/*
Output:
    f(x, y) = 1 * ((3)1) * ((2)1) * (x+(-1)y^2) * (x^3+y+(-1)1)
*/
```

so the factorization of \\(f(x, y)\\) is

\\[f(x, y) = 2 \times 3 \times (x^3 + y - 1) \times (y^2 - x)\\]

## P-adic Root Finding

Find the \\(2\\)-adic square roots of \\(17\\).

```rust
use algebraeon::nzq::{Natural, Integer};
use algebraeon::rings::{polynomial::*, structure::*};
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

Truncating to the last 16 bits it can be verified that, modulo \\(2^{16}\\), the square of these values is \\(17\\).

```rust
let a = 0b0010011011101001u16;
assert_eq!(a.wrapping_mul(a), 17u16);
let b = 0b1101100100010111u16;
assert_eq!(b.wrapping_mul(b), 17u16);
```

## Enumerating a Finitely Generated Group

Let \\(G\\) be the finitely generated group generated by \\(3\\) generators \\(a\\), \\(b\\), \\(c\\) subject to the relations \\(a^2 = b^2 = c^2 = (ab)^3 = (bc)^5 = (ac)^2 = e\\).

\\[G = \langle a, b, c : a^2 = b^2 = c^2 = (ab)^3 = (bc)^5 = (ac)^2 = e \rangle\\]

Using Algebraeon, \\(G\\) is found to be a finite group of order \\(120\\):

```rust
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

```rust
use algebraeon::nzq::{Rational};
use algebraeon::rings::{matrix::*, isolated_algebraic::*};
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

```rust
use algebraeon::rings::polynomial::*;
use algebraeon::nzq::Integer;

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

\\[\mathop{\text{disc}}(ax^2 + bx + c) = b^2 - 4ac\\]

\\[\mathop{\text{disc}}(ax^3 + bx^2 + cx + d) = b^2c^2 - 4ac^3 - 4b^3d - 27a^2d^2 + 18abcd\\]