use std::{collections::{HashSet}};


mod sets;
mod groups;
mod numbers;

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

    let iso = groups::iso_rep::IsoRep::Symmetric(5);// * groups::IsoRep::Dihedral(3);
    println!("{:?} ", iso);
    let grp = iso.to_group().unwrap();

    println!("");

    let mut isom_classes = HashSet::new();
    for (sg, _gens) in grp.subgroups() {
        let isom_class = groups::iso_rep::IsoRep::from_group(&sg.to_group());
        isom_classes.insert(isom_class);
    }
    for isom_class in isom_classes {
        print!("{} ", isom_class.to_string());
    }

    println!("");
}
