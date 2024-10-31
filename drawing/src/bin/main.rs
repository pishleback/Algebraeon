use malachite_nz::integer::Integer;
use malachite_q::Rational;
use orthoclase_rings::{
    polynomial::polynomial::Polynomial, ring_structure::cannonical::*, structure::StructuredType,
};

fn main() {
    println!("{:?}", Rational::abs(&Rational::from(-10)));

    let x = &Polynomial::<Integer>::var().into_ring();
    println!("{}", (x.pow(12) - 1).into_set().factor().unwrap());
}
