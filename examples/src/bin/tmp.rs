#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::{integer::*, natural::*, rational::*},
    rings::{
        number::natural::factored::factor,
        polynomial::{multipoly::*, polynomial::Polynomial},
        structure::{elements::*, structure::*},
    },
};
use std::str::FromStr;

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let n = Natural::from_str("435879897345023450789532407859897340783457803453420987534278053978411111").unwrap();
    println!("{}", n);
    println!("{:?}", factor(n));
}
