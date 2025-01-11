#![allow(dead_code, warnings)]

use algebraeon_groups::free_group::todd_coxeter::*;

fn main() {
    let mut g = FinitelyGeneratedGroupPresentation::new();
    // Add the 3 generators
    let a = g.add_generator();
    let b = g.add_generator();
    let c = g.add_generator();
    // Add the relations
    g.add_relation(a.pow(2));
    g.add_relation(b.pow(2));
    g.add_relation(c.pow(2));
    g.add_relation((&a * &b).pow(3));
    g.add_relation((&b * &c).pow(5));
    g.add_relation((&a * &c).pow(2));
    // Count elements
    let (n, _) = g.enumerate_elements();
    assert_eq!(n, 120);
}
