use std::str::FromStr;

use malachite_base::num::conversion::traits::IntegerMantissaAndExponent;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rings::root_tools::*;
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
    //f has roots 2+3i, 2-3i
    let f = ZZ_POLY.from_coeffs(vec![
        Integer::from(13),
        Integer::from(-4),
        Integer::from(1),
    ]);
    let n = ZZ_POLY.count_complex_roots(
        &f,
        &Rational::from(2),
        &Rational::from(3),
        &Rational::from(3),
        &Rational::from(4),
    );
    println!("{:?}", n);
}
