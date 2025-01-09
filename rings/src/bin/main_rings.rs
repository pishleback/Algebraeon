#![allow(dead_code, warnings)]

use algebraeon_rings::{
    number::algebraic::isolated_roots::padic::*,
    polynomial::polynomial::*,
    ring_structure::{elements::*, structure::*},
};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

fn main() {
    // let x = &Polynomial::<Integer>::var().into_ergonomic();
    // let f = (x.pow(30) - 1).into_verbose();

    // // let f = p1();

    // println!("{}", f);

    // println!("{}", f.factor().unwrap());

    // println!(
    //     "{:?}",
    //     aks_primality_test(&Natural::from_str("10498718947").unwrap())
    // );

    let mut a = PAdicAlgebraic::Rational(PAdicRational {
        p: Natural::from(2u32),
        rat: Rational::from_integers(Integer::from(-2), Integer::from(3)),
    });

    println!("{:?}", a.reduce_modulo_valuation(10).digits());
}
