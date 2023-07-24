use std::collections::HashSet;

use malachite_nz::integer::Integer;
use malachite_q::Rational;
use rings::lattice::LinearLattice;

use crate::rings::matrix::Matrix;

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

    let mat = Matrix::from_rows(vec![
        vec![Integer::from(2), Integer::from(1), Integer::from(2), Integer::from(2)],
        vec![Integer::from(2), Integer::from(2), Integer::from(1), Integer::from(2)],
        vec![Integer::from(2), Integer::from(2), Integer::from(2), Integer::from(1)],
    ]);

    let y = Matrix::from_rows(vec![
        vec![Integer::from(1)],
        vec![Integer::from(1)],
        vec![Integer::from(1)],
    ]);

    mat.pprint();

    let x = mat.col_solve(&y).unwrap();
    x.pprint();

    Matrix::mul_refs(&mat, &x).unwrap().pprint();

    // mat.pprint();

    // let lattice = mat.col_span();

    // lattice.pprint();
}
