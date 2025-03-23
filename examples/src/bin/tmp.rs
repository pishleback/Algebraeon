#![allow(dead_code, warnings)]

use algebraeon::{
    nzq::integer::*,
    rings::{
        polynomial::multipoly::*,
        structure::{elements::*, structure::*},
    },
};

fn main() {
    unsafe {
        std::env::set_var("RUST_BACKTRACE", "1");
    }
}
