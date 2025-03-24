#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::integer::*,
    rings::{
        polynomial::multipoly::*,
        structure::{elements::*, structure::*},
    },
};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
    let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();

    let f = (x.pow(120) - 1).into_verbose();
    println!("f = {}", f.factor().unwrap());
}
