#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::{integer::*, natural::*, rational::*},
    rings::{
        number::natural::factorization::factor,
        polynomial::{multipoly::*, polynomial::Polynomial},
        structure::{elements::*, structure::*},
    },
};
use std::str::FromStr;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let n = Natural::from_str("89754902345089754074587943578934757489025340982435043289079823459780598711")
        .unwrap();
    println!("{}", n);
    println!("{:?}", factor(n));
}
