#![allow(dead_code, warnings)]

use std::str::FromStr;

use algebraeon_rings::structure::structure::*;
use malachite_nz::{integer::Integer, natural::Natural};

fn main() {
    println!(
        "{}",
        Integer::from_str("2")
            .unwrap()
            .nat_pow(&Natural::from_str("321985732985472398457").unwrap())
            .factor()
            .unwrap()
    );
}
