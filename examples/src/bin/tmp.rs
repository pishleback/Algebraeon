#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon::nzq::integer::*;
use algebraeon::nzq::natural::*;
use algebraeon::nzq::rational::*;
use algebraeon::rings::number::finite_fields::extension::f9;
use algebraeon::rings::polynomial::multipoly::MultiPolynomial;
use algebraeon::rings::polynomial::multipoly::Variable;
use algebraeon::rings::polynomial::polynomial::*;
use algebraeon::rings::structure::elements::IntoErgonomic;
use algebraeon::rings::structure::structure::*;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x_var = Variable::new("x");
    let y_var = Variable::new("y");

    let x = MultiPolynomial::<Integer>::var(x_var.clone()).into_ergonomic();
    let y = MultiPolynomial::<Integer>::var(y_var.clone()).into_ergonomic();

    let f = ((y.pow(2) - x.pow(2)) * 24).into_verbose();
    let (u, g) = f.factor_fav_assoc();

    println!("f = {}", f);
    println!("u g = {} {}", u, g);

    println!("{}", f.factor().unwrap());

    // let x = Integer::from_str("32198573291847394729843798732185472398457").unwrap();
    // let f = x.factor().unwrap();
    // println!("{}", f);
}
