#![allow(dead_code, warnings)]

use algebraeon_nzq::integer::*;
use algebraeon_nzq::natural::*;
use algebraeon_nzq::rational::*;
use algebraeon_rings::{
    polynomial::polynomial::*,
    structure::{elements::*, structure::*},
};

fn main() {
    let x = &Polynomial::<Integer>::var().into_ergonomic();
    let f = (7 * x.pow(1) + 7 * x.pow(3) + x.pow(4) + x.pow(6)).into_verbose();

    // let f = p1();

    println!("{}", f);
    println!("{}", f.factor().unwrap());

    // println!(
    //     "{:?}",
    //     aks_primality_test(&Natural::from_str("10498718947").unwrap())
    // );
}
