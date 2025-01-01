#![allow(dead_code, warnings, unused)]

use std::marker::PhantomData;
use std::str::FromStr;
use std::task::Poll;

use algebraeon_groups::examples::symmetric::*;
use algebraeon_groups::group::*;
use algebraeon_rings::elements::*;
use algebraeon_rings::linear::matrix::*;
use algebraeon_rings::number::finite_fields::modulo::*;
use algebraeon_rings::polynomial::multipoly::*;
use algebraeon_rings::polynomial::polynomial::*;
use algebraeon_rings::polynomial::polynomial::*;
use algebraeon_rings::ring_structure::quotient::*;
use algebraeon_rings::ring_structure::structure::*;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;

fn main() {
    let a = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ring() + 1).into_set().into();
            x.into()
        }
    });

    let a = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
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
    let b = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
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
    let c = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
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
