#![allow(dead_code, warnings)]

use algebraeon_rings::elements::*;
use algebraeon_rings::number::natural::*;
use algebraeon_rings::{polynomial::polynomial::Polynomial, ring_structure::cannonical::*};
use malachite_base::num::basic::traits::One;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;

fn main() {
    // let x = &Polynomial::<Integer>::var().into_ring();
    // println!("{}", (x.pow(12) - 1).into_set().factor().unwrap());

    for x in 0..10usize {
        for n in 1..10usize {
            let x = Natural::from(x);
            let n = Natural::from(n);
            let r = nth_root(&x, &n);
            assert!(pow(&r, &n) <= n);
            assert!(x == 0 || n < pow(&(r + Natural::ONE), &n));
        }
    }
}
