use std::rc::Rc;

use malachite_base::num::conversion::traits::FromSciString;
use malachite_nz::integer::Integer;
use malachite_q::Rational;
use orthoclase_all::rings::{
    polynomial::polynomial::Polynomial, ring_structure::cannonical::*, structure::StructuredType,
};

#[derive(Debug)]
struct A {
    n: usize,
}

#[derive(Debug)]
struct B<ABorrow : std::borrow::Borrow<A>> {
    a: ABorrow,
}

fn thing() -> (Rc<A>, B<Rc<A>>, B<Rc<A>>) {
    let a = Rc::new(A { n: 4 });
    (a.clone(), B { a : a.clone() }, B { a : a.clone() })
}

fn main() {
    let (a, b1, b2) = thing();
    let b3 = B{a : a.clone()};

    println!("{:?}", b1);
    println!("{:?}", b2);
    println!("{:?}", b3);
    println!("{:?}", a);


    println!("hi");

    println!("{:?}", Rational::abs(&Rational::from(-10)));

    let x = &Polynomial::<Integer>::var().into_ring();
    println!("{}", (x.pow(12) - 1).into_set().factor().unwrap());
}
