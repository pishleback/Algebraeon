#![allow(dead_code, warnings)]

use algebraeon_rings::{
    number::algebraic::padic::*,
    polynomial::polynomial::Polynomial,
    structure::{elements::*, structure::*},
};
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

fn main() {
    use algebraeon_rings::{polynomial::polynomial::*, structure::elements::*};
    use malachite_nz::{integer::Integer, natural::Natural};
    let x = Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(2) - 17).into_verbose();
    for mut root in f.all_padic_roots(&Natural::from(2u32)) {
        println!("{}", root.truncate(&20.into()).string_repr()); // Show 20 2-adic digits
    }
}
