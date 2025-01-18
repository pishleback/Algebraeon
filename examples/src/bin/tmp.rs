#![allow(dead_code, warnings)]

use algebraeon_rings::{
    number::algebraic::padic::*,
    polynomial::polynomial::Polynomial,
    structure::{elements::*, structure::*},
};
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;
use structure::PAdicAlgebraicStructure;

fn main() {
    let ring = PAdicAlgebraicStructure::new(Natural::from(5u32));
    let x = Polynomial::<Integer>::var().into_ergonomic();

    let a = {
        let f = (x.pow(3) - 3 * x.pow(2) - x.pow(1) + 1).into_verbose();
        let r = f.all_padic_roots(&Natural::from(5u32));
        assert_eq!(r.len(), 1);
        r.into_iter().next().unwrap()
    };

    let b = {
        let f = (x.pow(4) + x.pow(2) - 2 * x.pow(1) - 1).into_verbose();
        let r = f.all_padic_roots(&Natural::from(5u32));
        assert_eq!(r.len(), 1);
        r.into_iter().next().unwrap()
    };

    let c = {
        let f = (x.pow(5) + x.pow(2) + 2 * x.pow(1) + 1).into_verbose();
        let r = f.all_padic_roots(&Natural::from(5u32));
        assert_eq!(r.len(), 1);
        r.into_iter().next().unwrap()
    };

    let d = PAdicAlgebraic::from_rational(
        Natural::from(5u32),
        Rational::from_integers(Integer::from(2), Integer::from(7)),
    );

    let e = PAdicAlgebraic::from_rational(
        Natural::from(5u32),
        Rational::from_integers(Integer::from(-2), Integer::from(1)),
    );

    println!("a = {}", a);
    println!("b = {}", b);
    println!("c = {}", c);
    println!("d = {}", d);
    println!("e = {}", e);

    println!("-a = {}", ring.neg(&a));
    println!("-b = {}", ring.neg(&b));
    println!("-c = {}", ring.neg(&c));
    println!("-d = {}", ring.neg(&d));
    println!("-e = {}", ring.neg(&e));

    println!("a+b = {}", ring.add(&a, &b));
    println!("a+c = {}", ring.add(&a, &c));
    println!("d+b = {}", ring.add(&d, &b));
    println!("d+c = {}", ring.add(&d, &c));
    println!("c+c = {}", ring.add(&c, &c));

    println!("a*b = {}", ring.mul(&a, &b));
    println!("a*c = {}", ring.mul(&a, &c));
    println!("d*b = {}", ring.mul(&d, &b));
    println!("d*c = {}", ring.mul(&d, &c));
    println!("c*c = {}", ring.mul(&c, &c));
    println!("a*e = {}", ring.mul(&a, &e));

    println!("a^-1 = {}", ring.inv(&a).unwrap());
    println!("b^-1 = {}", ring.inv(&b).unwrap());
    println!("c^-1 = {}", ring.inv(&c).unwrap());
    println!("d^-1 = {}", ring.inv(&d).unwrap());
    println!("e^-1 = {}", ring.inv(&e).unwrap());
}
