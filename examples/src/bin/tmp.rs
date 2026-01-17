use std::str::FromStr;

use algebraeon::{
    nzq::Rational,
    rings::matrix::{Matrix, StandardInnerProduct},
    sets::structure::MetaType,
};

fn main() {
    let m = Matrix::<Rational>::from_rows(vec![
        vec![
            Rational::from_str("1").unwrap(),
            Rational::from_str("0").unwrap(),
            Rational::from_str("0").unwrap(),
            Rational::from_str("10000").unwrap(),
        ],
        vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("1").unwrap(),
            Rational::from_str("0").unwrap(),
            Rational::from_str("16180").unwrap(),
        ],
        vec![
            Rational::from_str("0").unwrap(),
            Rational::from_str("0").unwrap(),
            Rational::from_str("1").unwrap(),
            Rational::from_str("26180").unwrap(),
        ],
    ]);

    m.pprint();

    let h = m.clone().lll_row_reduction_algorithm(
        &StandardInnerProduct::new(Rational::structure()),
        &Rational::from_str("3/4").unwrap(),
    );

    h.pprint();

    Matrix::<Rational>::mul(&h, &m).unwrap().pprint();
}
