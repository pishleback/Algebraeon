#![allow(dead_code, warnings)]

use algebraeon_rings::elements::*;
use algebraeon_rings::{polynomial::polynomial::Polynomial, ring_structure::cannonical::*};
use malachite_nz::integer::Integer;
use malachite_q::Rational;

fn main() {
    let x = &Polynomial::<Integer>::var().into_ring();
    println!("{}", (x.pow(12) - 1).into_set().factor().unwrap());
}
