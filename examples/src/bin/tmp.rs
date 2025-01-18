#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon_rings::{
    number::{algebraic::padic::*, natural::{factor::factor_by_try_divisors, primes::{aks_primality_test, is_prime}}},
    polynomial::polynomial::Polynomial,
    structure::{elements::*, structure::*},
};
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;
use structure::PAdicAlgebraicStructure;

fn main() {
    println!(
        "{:?}",
        is_prime(&Natural::from_str("10947810912894721490981347547827").unwrap())
    );
}
