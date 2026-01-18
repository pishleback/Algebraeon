use crate::polynomial::*;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;

mod gram_schmidt;
mod hermite_reduction;
mod inner_product;
mod jordan_normal_form;
mod lll_reduction;
#[allow(clippy::module_inception)]
mod matrix;
mod polynomial;
mod primitive;
mod row_operations;
mod smith_normal_form;
mod symmetric_matrix;

pub use hermite_reduction::*;
pub use inner_product::*;
pub use jordan_normal_form::*;
pub use matrix::*;
use row_operations::*;
pub use symmetric_matrix::*;

/*
//for LLL testing from wikipeida
let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
    vec![
        (-2 + 2 * i).into_verbose(),
        (7 + 3 * i).into_verbose(),
        (7 + 3 * i).into_verbose(),
        (-5 + 4 * i).into_verbose(),
    ],
    vec![
        (3 + 3 * i).into_verbose(),
        (-2 + 4 * i).into_verbose(),
        (6 + 2 * i).into_verbose(),
        (-1 + 4 * i).into_verbose(),
    ],
    vec![
        (2 + 2 * i).into_verbose(),
        (8 + 0 * i).into_verbose(),
        (-9 + 1 * i).into_verbose(),
        (-7 + 5 * i).into_verbose(),
    ],
    vec![
        (8 + 2 * i).into_verbose(),
        (-9 + 0 * i).into_verbose(),
        (6 + 3 * i).into_verbose(),
        (-4 + 4 * i).into_verbose(),
    ],
]);
*/
