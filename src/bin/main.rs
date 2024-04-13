#![allow(dead_code, warnings, unused)]

#[macro_use]
extern crate glium;

use std::marker::PhantomData;
use std::str::FromStr;
use std::task::Poll;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use orthoclase_all::drawing::canvas2d::*;
use orthoclase_all::drawing::Canvas;
use orthoclase_all::geometry::convex_simplicial_complex::*;
use orthoclase_all::geometry::vector::*;
use orthoclase_all::geometry::*;
use orthoclase_all::groups::group::*;
use orthoclase_all::groups::permutation::*;
use orthoclase_all::rings::linear::matrix::*;
use orthoclase_all::rings::number::algebraic::isolated_roots::*;
use orthoclase_all::rings::number::algebraic::number_field::*;
use orthoclase_all::rings::number::modulo::*;
use orthoclase_all::rings::polynomial::multipoly::*;
use orthoclase_all::rings::polynomial::polynomial::*;
use orthoclase_all::rings::ring_structure::cannonical::*;
use orthoclase_all::rings::ring_structure::quotient::*;
use orthoclase_all::rings::ring_structure::structure::*;
use orthoclase_all::rings::structure::*;
use orthoclase_all::rings::polynomial::polynomial::*;

fn main() {
    let x = &Polynomial::<Rational>::var().into_ring();
    let f = (x.pow(2) + 1).into_set();
    let (k, roots) = splitting_field_anf(&f);
    println!("{}", k.modulus());
    for root in roots {
        println!("root = {}", root);
    }
}
