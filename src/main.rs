use std::str::FromStr;

use malachite_base::num::conversion::traits::IntegerMantissaAndExponent;
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

fn todo() {
    // let s = NaturalPrimeGenerator::new();
    // for p in s {
    //     println!("{}", p);
    // }

    // let x = Integer::from(10000);
    // let f = ZZ.factor(&x);
    // println!("{:?}", f);

    let f = ZZ_POLY.from_coeffs(vec![
        Integer::from(1),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(0),
        Integer::from(1),
    ]);
    let roots = ZZ_POLY.all_complex_roots(&f);

    let a = QQ_BAR.sum(roots);

    println!("{:?}", a);
}

fn main() {
    let modn = EuclideanQuotient::new_ring(ZZ, Integer::from(100));

    println!("{:?}", modn.div_refs(&Integer::from(75), &Integer::from(85)));

    let poly_modn = PolynomialRing::new(&modn);

    let f = poly_modn.from_coeffs(vec![Integer::from(-2), Integer::from(3), Integer::from(4)]);
    let g = poly_modn.from_coeffs(vec![Integer::from(1), Integer::from(4), Integer::from(6)]);

    println!("f = {}", poly_modn.to_string(&f));
    println!("g = {}", poly_modn.to_string(&g));

    let h = poly_modn.add_refs(&f, &g);
    println!("h = {}", poly_modn.to_string(&h));
}
