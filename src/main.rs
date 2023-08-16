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
        Integer::from(-2),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = ZZ_POLY.all_real_roots(&f);

    let a = &roots[0];
    let b = &roots[1];

    println!("a = {}", QQ_BAR_REAL.to_string(a));
    println!("b = {}", QQ_BAR_REAL.to_string(b));

    let c = QQ_BAR_REAL.mul(a.clone(), b.clone());

    println!("c = {}", QQ_BAR_REAL.to_string(&c));
    println!("{:?}", c);
}
