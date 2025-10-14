# Matrices

## Creating Matrices

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::matrix::Matrix;

let a = Matrix::<Integer>::from_rows(vec![
    vec![Integer::from(1), Integer::from(2), Integer::from(3)],
    vec![Integer::from(0), Integer::from(-1), Integer::from(4)],
]);

assert_eq!(a.rows(), 2);
assert_eq!(a.cols(), 3);

let zero = Matrix::<Integer>::zero(2, 3);
let ident = Matrix::<Integer>::ident(3);
let diag = Matrix::<Integer>::diag(&vec![
    Integer::from(1),
    Integer::from(5),
    Integer::from(9),
]);

assert_eq!(zero.rows(), 2);
assert_eq!(ident.cols(), 3);
assert_eq!(
    diag.get_row(0),
    vec![Integer::from(1), Integer::from(0), Integer::from(0)]
);
```

## Basic Arithmetic

```rust
use algebraeon::nzq::Integer;
use algebraeon::rings::matrix::Matrix;

let a = Matrix::<Integer>::from_rows(vec![
    vec![Integer::from(1), Integer::from(2), Integer::from(3)],
    vec![Integer::from(0), Integer::from(1), Integer::from(4)],
]);
let b = Matrix::<Integer>::from_rows(vec![
    vec![Integer::from(1), Integer::from(0)],
    vec![Integer::from(0), Integer::from(1)],
    vec![Integer::from(2), Integer::from(3)],
]);

let sum = Matrix::add(&a, &a).unwrap();
let product = Matrix::mul(&a, &b).unwrap();

assert_eq!(
    sum.get_row(0),
    vec![Integer::from(2), Integer::from(4), Integer::from(6)]
);
assert_eq!(
    product.get_row(0),
    vec![Integer::from(7), Integer::from(11)]
);
assert_eq!(
    product.get_row(1),
    vec![Integer::from(8), Integer::from(13)]
);

let scaled = a.clone().mul_scalar(&Integer::from(3));
let transposed = b.clone().transpose();

assert_eq!(
    scaled.get_row(1),
    vec![Integer::from(0), Integer::from(3), Integer::from(12)]
);
assert_eq!(
    transposed.get_row(0),
    vec![Integer::from(1), Integer::from(0), Integer::from(2)]
);
```

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::matrix::Matrix;
let v = Matrix::<Integer>::from_rows(vec![vec![
    Integer::from(1),
    Integer::from(2),
    Integer::from(3),
]]);
let w = Matrix::<Integer>::from_rows(vec![vec![
    Integer::from(4),
    Integer::from(5),
    Integer::from(6),
]]);
assert_eq!(Matrix::dot(&v, &w), Integer::from(32));
```

## Determinants, Rank, and Inverses

```rust
use algebraeon::nzq::{Integer, Rational};
use algebraeon::rings::matrix::Matrix;

let integer_matrix = Matrix::<Integer>::from_rows(vec![
    vec![Integer::from(1), Integer::from(2)],
    vec![Integer::from(3), Integer::from(5)],
]);
assert_eq!(
    integer_matrix.clone().det().unwrap(),
    Integer::from(-1)
);
assert_eq!(integer_matrix.rank(), 2);

let rational_matrix = Matrix::<Rational>::from_rows(vec![
    vec![Rational::from(2), Rational::from(4), Rational::from(4)],
    vec![Rational::from(-6), Rational::from(6), Rational::from(12)],
    vec![Rational::from(10), Rational::from(7), Rational::from(17)],
]);
let inverse = Matrix::inv(&rational_matrix).unwrap();
assert_eq!(Matrix::mul(&rational_matrix, &inverse).unwrap(), Matrix::ident(3));
```
