#![allow(dead_code, warnings)]

use algebraeon_rings::{polynomial::polynomial::*, structure::elements::*};
use malachite_nz::integer::Integer;

fn main() {
    let x = &Polynomial::<Integer>::var().into_ergonomic();
    let f = ((x.pow(5) - x + 1) * (x.pow(2) + 1)).into_verbose();
    // Find the complex roots of f(x)
    for root in f.all_complex_roots() {
        println!("root {} of degree {}", root, root.degree());
    }
    /*
    Output:
        root -i of degree 2
        root i of degree 2
        root ≈-1.172 of degree 5
        root ≈-0.179-1.087i of degree 5
        root ≈-0.179+1.087i of degree 5
        root ≈0.767-0.352i of degree 5
        root ≈0.767+0.352i of degree 5
    */
}
