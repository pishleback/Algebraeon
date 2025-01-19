#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon_rings::{
    number::{
        algebraic::padic::*,
        finite_fields::modulo::Modulo,
        natural::{
            factor::factor_by_try_divisors,
            primes::{aks_primality_test, is_prime},
        },
    },
    polynomial::polynomial::Polynomial,
    structure::{elements::*, structure::*},
};
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;
use structure::PAdicAlgebraicStructure;

fn main() {
    let x = &Polynomial::<Modulo<3>>::var().into_ergonomic();
    let p = x.pow(101) - 1;
    let p = p.into_verbose();
    let f = p.factor().unwrap();
    println!("{} = {}", p, f);
}
