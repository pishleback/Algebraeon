use std::str::FromStr;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use rings::ergonomic::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {
    // let f = interpolate_by_lagrange_basis::<Rational>(&vec![
    //     (Rational::from(1), Rational::from(2)),
    //     (Rational::from(2), Rational::from(1)),
    //     (Rational::from(3), Rational::from(4)),
    // ])
    // .unwrap();
    // println!("{}", f.to_string());
    // let x = &Ergonomic::new(Polynomial::<Rational>::var());

    // let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
    // let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();

    // let r1 = Polynomial::rem_refs(&f, &g).unwrap();
    // println!("r1 = {}", r1.to_string());
    // let r2 = Polynomial::rem_refs(&g, &r1).unwrap();
    // println!("r2 = {}", r2.to_string());
    // let r3 = Polynomial::rem_refs(&r1, &r2).unwrap();
    // println!("r3 = {}", r3.to_string());
    // let r4 = Polynomial::rem_refs(&r2, &r3).unwrap();
    // println!("r4 = {}", r4.to_string());

    // println!("---------------");

    // let x = &Ergonomic::new(Polynomial::<Integer>::var());
    // let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
    // let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();

    // println!("{}", f.to_string());
    // println!("{}", g.to_string());
    // println!(
    //     "subres6 = {}",
    //     subresultant_naive(8, 6, &f, &g, 6).to_string()
    // );
    // println!(
    //     "subres5 = {}",
    //     subresultant_naive(8, 6, &f, &g, 5).to_string()
    // );
    // println!(
    //     "subres4 = {}",
    //     subresultant_naive(8, 6, &f, &g, 4).to_string()
    // );
    // println!(
    //     "subres3 = {}",
    //     subresultant_naive(8, 6, &f, &g, 3).to_string()
    // );
    // println!(
    //     "subres2 = {}",
    //     subresultant_matrix(8, 6, &f, &g, 2)
    //         .apply_map(|x| x.apply_map(|y| Rational::from(y)))
    //         .det()
    //         .unwrap()
    //         .to_string()
    // );
    // println!(
    //     "subres1 = {}",
    //     subresultant_matrix(8, 6, &f, &g, 1)
    //         .apply_map(|x| x.apply_map(|y| Rational::from(y)))
    //         .det()
    //         .unwrap()
    //         .to_string()
    // );
    // println!(
    //     "subres0 = {}",
    //     subresultant_matrix(8, 6, &f, &g, 0)
    //         .apply_map(|x| x.apply_map(|y| Rational::from(y)))
    //         .det()
    //         .unwrap()
    //         .to_string()
    // );

    // let f = (x+7).elem();
    // let g = ((x + 3).pow(3)).elem();

    // let res = Polynomial::pseudo_gcd(f, g);
    // println!("{}", res.to_string());

    println!("hello 71 tests");
}
