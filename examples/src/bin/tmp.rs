#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::{integer::*, rational::*},
    rings::{
        number::anf::number_field::new_anf,
        polynomial::{multipoly::*, polynomial::Polynomial},
        structure::{elements::*, structure::*},
    },
};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x = &Polynomial::<Rational>::var().into_ergonomic();
    let f = (x.pow(5) - x + 1).into_verbose();
    let anf = new_anf(f);
    println!("(r, s) = {:?}", anf.signature());
}
