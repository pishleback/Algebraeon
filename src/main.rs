use std::str::FromStr;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rings::ergonomic::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {
    let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());

    let f = ((x.pow(2) + 3 * x + 1).pow(300)).elem();

    println!("{}", ZZ_POLY.to_string(&f));

    let a = ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap();

    println!("{}", a.to_string(&ZZ_POLY))
}
