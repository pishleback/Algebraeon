#![allow(dead_code, warnings)]

use algebraeon::{nzq::integer::*, rings::{polynomial::multipoly::*, structure::{elements::*, structure::*}}};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
    let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();

    let f = (6 * (x.pow(4) - x.pow(3) * y.pow(2) + x * y - x - y.pow(3) + y.pow(2))).into_verbose();
    println!("f(x, y) = {}", f.factor().unwrap());

    /*
    Output:
        f(x, y) = 1 * ((3)1) * ((2)1) * (x+(-1)y^2) * (x^3+y+(-1)1)
    */
}
