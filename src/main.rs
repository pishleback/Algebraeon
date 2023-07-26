use std::str::FromStr;

use malachite_nz::integer::Integer;
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
    let a = Integer::from(120);
    println!("{:?}", a.factor());

    // let x = &Ergonomic::new(Polynomial::<Integer>::var());
    // let a = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();

    // let fs = a.factor().unwrap();
    // println!("{}", fs.unit().to_string());
    // for (p, k) in fs.factors() {
    //     println!("{} ^ {}", p.to_string(), k.to_string());
    // }
}
