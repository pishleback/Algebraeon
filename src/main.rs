use std::str::FromStr;

use malachite_base::num::conversion::traits::IntegerMantissaAndExponent;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rings::algebraic::*;
use rings::ergonomic::*;
use rings::multipoly::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {
    let f = ZZ_POLY.from_coeffs(vec![
        Integer::from(-1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let n = ZZ_POLY.count_complex_roots(
        &f,
        &Rational::from(-2),
        &Rational::from(2),
        &Rational::from(-2),
        &Rational::from(2),
    );
    println!("{:?}", n);
}
