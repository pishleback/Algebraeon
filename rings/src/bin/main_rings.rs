#![allow(dead_code, warnings)]

use algebraeon_rings::{polynomial::polynomial::*, ring_structure::{elements::*, structure::*}};
use malachite_nz::integer::Integer;


fn main() {
    let x = &Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(30) - 1).into_verbose();

    // let f = p1();

    println!("{}", f);


    println!("{}", f.factor().unwrap());

    // println!(
    //     "{:?}",
    //     aks_primality_test(&Natural::from_str("10498718947").unwrap())
    // );
}
