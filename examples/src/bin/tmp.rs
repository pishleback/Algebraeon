#![allow(dead_code, warnings)]

use std::str::FromStr;
use algebraeon::{nzq::natural::Natural, rings::number::natural::factorization::factor};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let n = Natural::from_str("11111111111111111111111111111111111111111111111111111111111111").unwrap();
    let f = factor(n.clone()).unwrap();
    println!("{} = {}", n, f);
}
