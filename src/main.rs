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

    let a = QQ_BAR.sum(roots.iter().collect());

    println!("{:?}", a);
}

fn main() {
    // (x+11)*(x-222)*(x-3333)*(x^2+x+1)*(x+67)

    let f = ZZ_POLY.product(vec![
        &ZZ_POLY.from_coeffs(vec![Integer::from(11), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(-222), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(-3333), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(1), Integer::from(1), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![
            Integer::from(-1),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
        &ZZ_POLY.from_coeffs(vec![Integer::from(67), Integer::from(1)]),
        &ZZ_POLY.from_coeffs(vec![
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ]),
    ]);
    println!("{}", ZZ_POLY.to_string(&f));

    let fs = ZZ_POLY.factorize_by_hensel_lifting(f);
    println!("{:?}", fs);

    // let modn = EuclideanQuotient::new_field_unchecked(ZZ, Integer::from(5));

    // println!("{:?}", modn.all_units());

    // let poly_modn = PolynomialRing::new(&modn);

    // // let f = poly_modn.from_coeffs(vec![
    // //     Integer::from(-1),
    // //     Integer::from(1),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(0),
    // //     Integer::from(1),
    // // ]);

    // // let fs = ZZ_POLY.factor(&f);
    // // println!("{:?}", fs);

    // let fs = poly_modn.factorize_by_trying_all_factors(f);
    // println!("{:?}", fs);
}
