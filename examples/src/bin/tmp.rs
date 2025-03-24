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
    let p = (x.pow(4) - x + 1).into_verbose();
    let (f, roots) = p.splitting_field();
    println!("{}", f.modulus());
    for r in roots {
        println!("{}", r);
    }
}
