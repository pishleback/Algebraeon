#![allow(dead_code, warnings)]

use std::str::FromStr;
use algebraeon::{
    nzq::{integer::*, natural::*, rational::*},
    rings::{
        number::natural::factored::factor, polynomial::{multipoly::*, polynomial::Polynomial}, structure::{elements::*, structure::*}
    },
};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }

    let n = Natural::from_str("3029750799235790328974398724798327893798253709351343459177777").unwrap();
    println!("{}", n);
    println!("{:?}", factor(n));
}
