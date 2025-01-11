#![allow(dead_code, warnings)]

fn main() {
    use algebraeon_rings::{linear::matrix::*, number::algebraic::isolated_roots::complex::*};
    use algebraeon_sets::structure::*;
    use malachite_q::Rational;
    // Construct a matrix
    let a = Matrix::<Rational>::from_rows(vec![
        vec![7, 5, -3, -2],
        vec![1, -1, -1, -1],
        vec![7, 4, -3, -6],
        vec![-1, 5, 1, 5],
    ]);
    // Put it into Jordan Normal Form
    #[cfg(not(debug_assertions))]
    let j = MatrixStructure::new(ComplexAlgebraic::structure()).jordan_normal_form(&a);
    #[cfg(not(debug_assertions))]
    j.pprint();
    /*
    Output:
        / -i√3    0      0    0 \
        | 0       i√3    0    0 |
        | 0       0      4    1 |
        \ 0       0      0    4 /
    */
}
