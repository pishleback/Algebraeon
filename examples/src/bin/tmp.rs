#![allow(dead_code, warnings)]

use algebraeon_rings::{polynomial::polynomial::*, structure::elements::*};
use malachite_nz::integer::Integer;

fn main() {
    let x = &Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(5) + x.pow(2) - x + 1).into_verbose();
    // Find the complex roots of f(x)
    for root in f.all_complex_roots() {
        println!("root {} of degree {}", root, root.degree());
    }
    /*
    Output:
        root ≈-1.328 of degree 3
        root ≈0.662-0.559i of degree 3
        root ≈0.662+0.559i of degree 3
        root -i of degree 2
        root i of degree 2
    */
}
