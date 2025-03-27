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

    let n =
        // Natural::from_str("3029750799235790328974398724798327893798253709351343459177777").unwrap();
        Natural::from_str("3146531246531241245132451321").unwrap();
    println!("{}", n);
    println!("{:?}", factor(n));
}
