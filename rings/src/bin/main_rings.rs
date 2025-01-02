#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon_rings::elements::*;
use algebraeon_rings::number::natural::*;
use algebraeon_rings::polynomial::polynomial::Polynomial;
use algebraeon_rings::ring_structure::structure::*;
use functions::*;
use malachite_base::num::basic::traits::One;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use primes::aks_primality_test;

fn main() {
    let x = &Polynomial::<Integer>::var().into_ring();
    // println!("{}", (x.pow(60) - 1).into_set().factor().unwrap());
    println!("{}", (x.pow(30) - 1).into_set().factor().unwrap());

    // println!(
    //     "{:?}",
    //     aks_primality_test(&Natural::from_str("10498718947").unwrap())
    // );
}
