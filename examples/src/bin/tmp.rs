#![allow(dead_code, warnings)]

use algebraeon_rings::{
    number::algebraic::isolated_roots::padic::*,
    polynomial::polynomial::Polynomial,
    structure::{elements::*, structure::*},
};
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

fn main() {
    let x = Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(6) - 2).into_verbose();
    println!("f = {}", f);
    assert_eq!(f.all_padic_roots(&Natural::from(2u32)).len(), 0);
    assert_eq!(f.all_padic_roots(&Natural::from(7u32)).len(), 0);
    assert_eq!(f.all_padic_roots(&Natural::from(727u32)).len(), 6);
    for mut root in f.all_padic_roots(&Natural::from(727u32)) {
        println!("{}", root);
    }
}
