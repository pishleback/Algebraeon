use std::str::FromStr;

use malachite_nz::integer::Integer;
use malachite_q::Rational;
use rings::ergonomic::*;
use rings::nzq::*;
use rings::poly::*;
use rings::ring::*;

mod groups;
mod numbers;
mod rings;
mod sets;

fn main() {
    // let (grp_g, _perms, _elems) = enumerated_groups::symmetric_group_structure(6);
    // let (grp_h, _perms, _elems) = enumerated_groups::symmetric_group_structure(6);
    // let grp_h = todd_coxeter::enumerate_group(
    //     3,
    //     vec![
    //         vec![0, 0],
    //         vec![2, 2],
    //         vec![4, 4],
    //         vec![0, 4, 0, 4],
    //         vec![0, 2, 0, 2, 0, 2],
    //         vec![2, 4, 2, 4, 2, 4, 2, 4, 2, 4],
    //     ],
    // ); //A5 x C2

    // let mut grp = groups::IsoRep::Symmetric(6).to_group().unwrap();
    // let grp = groups::quaternion_group_structure();
    // let grp_h = enumerated_groups::cyclic_group_structure(720);

    // let mut grp = todd_coxeter::enumerate_group(
    //     3,
    //     vec![
    //         vec![0, 2, 1, 3],
    //         vec![0, 4, 1, 5],
    //         vec![2, 4, 3, 5],
    //         vec![0, 0, 0, 0, 0],
    //         vec![2, 2, 2, 2, 2],
    //         vec![4, 4, 4, 4, 4],
    //     ],
    // );

    // let iso = groups::iso_rep::IsoRep::Symmetric(5);// * groups::IsoRep::Dihedral(3);
    // println!("{:?} ", iso);
    // let grp = iso.to_group().unwrap();

    // println!("");

    // let mut isom_classes = HashSet::new();
    // for (sg, _gens) in grp.subgroups() {
    //     let isom_class = groups::iso_rep::IsoRep::from_group(&sg.to_group());
    //     isom_classes.insert(isom_class);
    // }
    // for isom_class in isom_classes {
    //     print!("{} ", isom_class.to_string());
    // }

    // println!("");

    // let a = Matrix::from_rows(vec![
    //     vec![Integer::from(-2), Integer::from(0), Integer::from(0)],
    //     vec![Integer::from(0), Integer::from(-6), Integer::from(0)],
    //     vec![Integer::from(0), Integer::from(0), Integer::from(-120)],
    // ]);

    // a.pprint();

    // let (u, s, v, k) = a.smith_algorithm();

    // u.pprint();
    // s.pprint();
    // v.pprint();

    // let alat1 = AffineLattice::from_offset_and_linear_lattice(
    //     2,
    //     1,
    //     Matrix::from_rows(vec![vec![Integer::from(3)], vec![Integer::from(2)]]),
    //     LinearLattice::from_span(
    //         2,
    //         1,
    //         vec![
    //             Matrix::from_rows(vec![vec![Integer::from(5)], vec![Integer::from(0)]]),
    //             Matrix::from_rows(vec![vec![Integer::from(0)], vec![Integer::from(5)]]),
    //         ],
    //     ),
    // );

    // let alat2 = AffineLattice::from_offset_and_linear_lattice(
    //     2,
    //     1,
    //     Matrix::from_rows(vec![vec![Integer::from(1)], vec![Integer::from(1)]]),
    //     LinearLattice::from_span(
    //         2,
    //         1,
    //         vec![
    //             Matrix::from_rows(vec![vec![Integer::from(7)], vec![Integer::from(0)]]),
    //             Matrix::from_rows(vec![vec![Integer::from(0)], vec![Integer::from(7)]]),
    //         ],
    //     ),
    // );

    // alat1.check_invariants().unwrap();
    // alat2.check_invariants().unwrap();

    // alat1.pprint();
    // println!();
    // alat2.pprint();

    // let alat3 = AffineLattice::intersect_pair(2, 1, alat1, alat2);
    // println!();
    // alat3.pprint();

    // let f = Integer::from_str("168").unwrap().factor();
    // f.clone().unwrap().check_invariants().unwrap();
    // println!("{:?}", f);


    impl UniquelyFactorable for Integer {
        fn make_factorizer() -> Box<dyn UniqueFactorizer<Self>> {
            Box::new(NaiveIntegerFactorizer())
        }
    }

    // let a = Integer::from(120);
    // println!("{:?}", a.factor());

    impl UniquelyFactorable for Polynomial<Integer> {
        fn make_factorizer() -> Box<dyn UniqueFactorizer<Self>> {
            Box::new(KroneckerFactorizer())
        }
    }

    impl UniquelyFactorable for Polynomial<Rational> {
        fn make_factorizer() -> Box<dyn UniqueFactorizer<Self>> {
            Box::new(FractionFieldPolynomialFactorizer::new(KroneckerFactorizer()))
        }
    }

    let x = &Ergonomic::new(Polynomial::<Integer>::var());
    let a = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();

    let fs = a.factor().unwrap();
    println!("{}", fs.unit().to_string());
    for (p, k) in fs.factors() {
        println!("{} ^ {}", p.to_string(), k.to_string());
    }
}
