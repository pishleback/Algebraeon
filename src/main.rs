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
    // let f = interpolate_by_lagrange_basis::<Rational>(&vec![
    //     (Rational::from(1), Rational::from(2)),
    //     (Rational::from(2), Rational::from(1)),
    //     (Rational::from(3), Rational::from(4)),
    // ])
    // .unwrap();
    // println!("{}", f.to_string());

    let x = &Ergonomic::new(Polynomial::<Integer>::var());
    let a = ((1 + 2 * x).pow(10)).elem();

    let fs = factor_by_kroneckers_method(&a);
    for d in fs.divisors().unwrap() {
        println!("{}", d.to_string());
    }
}
