use std::str::FromStr;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rings::algebraic::*;
use rings::ergonomic::*;
use rings::multipoly::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {

    // for (c, k, h) in ZZ_POLY.collin_akritas(f) {
    //     println!("c={} k={} h={}", c, k, h);
    //     let h = {
    //         if h {
    //             Natural::from(1u8)
    //         } else {
    //             Natural::from(0u8)
    //         }
    //     };
    //     println!(
    //         "{}/{} {}/{}",
    //         c.clone(),
    //         Natural::from(1u8) << k,
    //         c + h,
    //         Natural::from(1u8) << k
    //     );
    // }
}
