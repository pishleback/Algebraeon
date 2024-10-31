#![allow(dead_code, warnings, unused)]

#[macro_use]
extern crate glium;

use std::marker::PhantomData;
use std::str::FromStr;
use std::task::Poll;

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;
use orthoclase_drawing::canvas::canvas2d::*;
use orthoclase_drawing::canvas::Canvas;
use orthoclase_groups::examples::symmetric::*;
use orthoclase_groups::group::*;
use orthoclase_rings::linear::matrix::*;
use orthoclase_rings::number::algebraic::isolated_roots::*;
use orthoclase_rings::number::algebraic::number_field::*;
use orthoclase_rings::number::modulo::*;
use orthoclase_rings::polynomial::multipoly::*;
use orthoclase_rings::polynomial::polynomial::*;
use orthoclase_rings::polynomial::polynomial::*;
use orthoclase_rings::ring_structure::cannonical::*;
use orthoclase_rings::ring_structure::quotient::*;
use orthoclase_rings::ring_structure::structure::*;
use orthoclase_rings::structure::*;

fn main() {
    let a = orthoclase_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ring() + 1).into_set().into();
            x.into()
        }
    });

    let a = orthoclase_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ring() + 1).into_set().into();
            x.into()
        }
    })
    .unwrap();
    let b = orthoclase_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ring() * 2).into_set().into();
            x.into()
        }
    })
    .unwrap();
    let c = orthoclase_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            0.into()
        } else if x == 0 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = x.inv().unwrap().neg().into();
            x.into()
        }
    })
    .unwrap();
    let d = Permutation::<24>::new_from_cycles(vec![
        vec![14, 17, 11, 19, 22],
        vec![20, 10, 7, 5, 21],
        vec![18, 4, 2, 6, 1],
        vec![8, 16, 13, 9, 12],
    ])
    .unwrap();

    println!("a = {}", a);
    println!("b = {}", b);
    println!("c = {}", c);
    println!("d = {}", d);

    let sg = Permutation::generated_finite_subgroup(vec![a, b, c, d]);
    println!("{:?}", sg.size())
}
