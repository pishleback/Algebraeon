#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon::nzq::integer::*;
use algebraeon::nzq::natural::*;
use algebraeon::nzq::rational::*;
use algebraeon::rings::number::finite_fields::extension::f9;
use algebraeon::rings::polynomial::polynomial::*;
use algebraeon::rings::structure::structure::*;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    // let x = Integer::from_str("32198573291847394729843798732185472398457").unwrap();
    // let f = x.factor().unwrap();
    // println!("{}", f);
}
