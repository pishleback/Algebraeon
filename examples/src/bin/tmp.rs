#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::{integer::*, rational::*},
    rings::{
        polynomial::{multipoly::*, polynomial::Polynomial},
        structure::{elements::*, structure::*},
    },
};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x = &Polynomial::<Rational>::var().into_ergonomic();
    let f = (x.pow(3) - x + 1).into_verbose();

    println!("f = {}", f.factor().unwrap());

    for root in f.all_complex_roots() {
        println!("{}", root);
    }
}
