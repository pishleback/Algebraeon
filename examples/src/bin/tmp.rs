#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon_nzq::integer::*;
use algebraeon_nzq::natural::*;
use algebraeon_nzq::rational::*;
use algebraeon_rings::structure::structure::*;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let x = Integer::from_str("32198573291847394729843798732185472398457").unwrap();
    let f = x.factor().unwrap();
    println!("{}", f);
}
