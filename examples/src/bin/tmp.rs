#![allow(dead_code, warnings)]


fn main() {
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
}
