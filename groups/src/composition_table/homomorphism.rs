use std::borrow::Borrow;
use std::collections::BTreeSet;
use std::collections::HashMap;

use super::generating_set::GeneratingSet;
use super::group::FiniteGroupMultiplicationTable;

#[derive(Clone)]
pub struct Homomorphism<
    DomainT: Borrow<FiniteGroupMultiplicationTable> + Clone,
    RangeT: Borrow<FiniteGroupMultiplicationTable> + Clone,
> {
    domain: DomainT,
    range: RangeT,
    func: Vec<usize>, //func : domain -> range    y = func[x]
}

impl<
    DomainT: Borrow<FiniteGroupMultiplicationTable> + Clone,
    RangeT: Borrow<FiniteGroupMultiplicationTable> + Clone,
> Homomorphism<DomainT, RangeT>
{
    pub fn check_state(&self) -> Result<(), &'static str> {
        //is function
        if self.func.len() != self.domain.borrow().size() {
            return Err("func size does not match domain size");
        }

        for x in self.domain.borrow().elems() {
            if self.func[x] >= self.range.borrow().size() {
                return Err("func image is too big for an element of the range");
            }
        }

        //is homomorphism
        for x in self.domain.borrow().elems() {
            for y in self.domain.borrow().elems() {
                if self.func[self.domain.borrow().mul(x, y)]
                    != self.range.borrow().mul(self.func[x], self.func[y])
                {
                    return Err("homomorphism does not respect composition");
                }
            }
        }

        Ok(())
    }

    pub fn new_unchecked(domain: DomainT, range: RangeT, func: Vec<usize>) -> Self {
        Self {
            domain,
            range,
            func,
        }
    }

    pub fn to_isomorphism(self) -> Option<Isomorphism<DomainT, RangeT>> {
        let n = self.domain.borrow().size();
        if n != self.range.borrow().size() {
            return None;
        }

        let mut inv = vec![None; n];
        #[allow(clippy::indexing_slicing)]
        for x in 0..n {
            match inv[self.func[x]] {
                Some(_y) => {
                    //self.func is not injective
                    return None;
                }
                None => inv[self.func[x]] = Some(x),
            }
        }

        let inv = inv
            .into_iter()
            .map(|x| if let Some(x_val) = x { x_val } else { panic!() })
            .collect::<Vec<usize>>();

        Some(Isomorphism {
            left_group: self.domain,
            right_group: self.range,
            left_func: self.func.clone(),
            right_func: inv,
        })
    }
}

#[derive(Clone)]
pub struct Isomorphism<
    LeftGrpT: Borrow<FiniteGroupMultiplicationTable> + Clone,
    RightGrpT: Borrow<FiniteGroupMultiplicationTable> + Clone,
> {
    left_group: LeftGrpT,
    right_group: RightGrpT,
    left_func: Vec<usize>,
    right_func: Vec<usize>,
}

impl<
    LeftGrpT: Borrow<FiniteGroupMultiplicationTable> + Clone,
    RightGrpT: Borrow<FiniteGroupMultiplicationTable> + Clone,
> Isomorphism<LeftGrpT, RightGrpT>
{
    pub fn check_state(&self) -> Result<(), &'static str> {
        let left_hom = Homomorphism {
            domain: self.left_group.borrow(),
            range: self.right_group.borrow(),
            func: self.left_func.clone(),
        };
        match left_hom.check_state() {
            Ok(()) => {}
            Err(msg) => {
                return Err(msg);
            }
        }

        let right_hom = Homomorphism {
            domain: self.right_group.borrow(),
            range: self.left_group.borrow(),
            func: self.right_func.clone(),
        };
        match right_hom.check_state() {
            Ok(()) => {}
            Err(msg) => {
                return Err(msg);
            }
        }

        if self.left_group.borrow().size() != self.right_group.borrow().size() {
            return Err("isomorphism group sizes dont match");
        }

        //are mutually inverse
        //only need to check one of left/right inverse because injective/surjective individually imply bijective once the sizes are the same
        for x in self.left_group.borrow().elems() {
            if x != self.right_func[self.left_func[x]] {
                return Err("isomorphism not inv");
            }
        }

        Ok(())
    }
}

//return an isomorphism from domain to range if one exists
//return None if no isomorphism exists
#[allow(clippy::too_many_lines)]
pub fn find_isomorphism<'a, 'b>(
    domain: &'a FiniteGroupMultiplicationTable,
    range: &'b FiniteGroupMultiplicationTable,
) -> Option<Isomorphism<&'a FiniteGroupMultiplicationTable, &'b FiniteGroupMultiplicationTable>> {
    let group_size = domain.size();
    if group_size != range.size() {
        return None;
    }

    //an element profile can be generated for each element in the range and domain
    //an isomorphism can only every map elements of the same profile to each other
    //use this is restrict the number of things to check
    #[allow(clippy::items_after_statements)]
    #[derive(Debug, PartialEq, Eq, Hash)]
    struct ElementProfile {
        order: usize,
        mul_order: BTreeSet<usize>, //helped a lot when comparing non-isomorphic groups of order 144
    }

    impl ElementProfile {
        fn new(x: usize, group: &FiniteGroupMultiplicationTable) -> Self {
            debug_assert!(x < group.size());
            ElementProfile {
                order: group.order(x).unwrap(),
                mul_order: group
                    .elems()
                    .map(|y| group.order(group.mul(x, y)).unwrap())
                    .collect(),
            }
        }
    }

    //want to quickly lookup [element -> profile] on the domain side
    let domain_elem_profiles: Vec<ElementProfile> = domain
        .elems()
        .map(|x| ElementProfile::new(x, domain))
        .collect();

    // println!("domain_elem_profiles {:?}", domain_elem_profiles);

    //want to quickly lookup [profile -> elements] on the range side
    let mut range_elem_profiles: HashMap<ElementProfile, Vec<usize>> = HashMap::new();
    for x in range.elems() {
        let p = ElementProfile::new(x, range);
        match range_elem_profiles.get_mut(&p) {
            Some(elems) => {
                elems.push(x);
            }
            None => {
                range_elem_profiles.insert(p, vec![x]);
            }
        }
    }

    // println!("range_elem_profiles {:?}", range_elem_profiles);

    struct GenInfo<'c> {
        gen_list: Vec<usize>,
        gen_set: GeneratingSet<'c>,
        image_options: Vec<Vec<usize>>,
        quant_to_check: usize,
    }

    fn find_new_gen_info<'c>(
        domain: &'c FiniteGroupMultiplicationTable,
        domain_elem_profiles: &[ElementProfile],
        range_elem_profiles: &HashMap<ElementProfile, Vec<usize>>,
    ) -> Result<GenInfo<'c>, ()> {
        let domain_gen_set = domain.generating_set();
        let domain_gens = domain_gen_set.gens();

        //TODO: compute a smaller subset of possible image points
        //e.g. each gen can only map to something of the same order
        let mut gen_image_options: Vec<Vec<usize>> = vec![];
        for g in domain_gens {
            match range_elem_profiles.get(&domain_elem_profiles[*g]) {
                Some(r_elems) => gen_image_options.push(r_elems.clone()),
                None => {
                    return Err(()); //there is no isomorphism
                }
            }
        }

        let mut num_to_check: usize = 1;
        for gen_image_option in &gen_image_options {
            num_to_check = match num_to_check.checked_mul(gen_image_option.len()) {
                Some(prod) => prod,
                None => usize::MAX,
            };
        }

        Ok(GenInfo {
            gen_list: domain_gens.clone(),
            gen_set: domain_gen_set,
            image_options: gen_image_options,
            quant_to_check: num_to_check,
        })
    }

    let mut current_gen_info;
    match find_new_gen_info(domain, &domain_elem_profiles, &range_elem_profiles) {
        Ok(gen_info) => current_gen_info = gen_info,
        Err(()) => {
            return None;
        }
    }

    // println!("group size {}", group_size);
    // println!("generating set {}", current_gen_info.quant_to_check);

    'outer_loop: loop {
        //loop over all possible images of the domain generators
        let mut image_option_counter = vec![0; current_gen_info.gen_list.len()];
        let mut already_checked: usize = 0;
        'loop_label: loop {
            //try to find a better generating set every few loops
            //better meaning there are fewer possible images to check than we currently have left to go
            if already_checked % 128 == 127 {
                // println!("maybe better generating set {}", already_checked);
                match find_new_gen_info(domain, &domain_elem_profiles, &range_elem_profiles) {
                    Ok(gen_info) => {
                        if gen_info.quant_to_check
                            < current_gen_info.quant_to_check - already_checked
                        {
                            // println!("better generating set {} {}", gen_info.quant_to_check, gen_info.gen_list.len());
                            current_gen_info = gen_info;
                            continue 'outer_loop;
                        }
                    }
                    Err(()) => {
                        return None;
                    }
                }
            }

            //compute homomorphism sending domain_gens -> image_gens
            if let Some(f) = current_gen_info
                .gen_set
                .generated_homomorphism(
                    &image_option_counter
                        .iter()
                        .enumerate()
                        .map(|(i, j)| current_gen_info.image_options[i][*j])
                        .collect(),
                    range,
                )
                .unwrap()
                && let Some(f_iso) = f.to_isomorphism()
            {
                return Some(f_iso);
            }
            already_checked += 1;

            //compute next gen_img
            'for_label: {
                for i in 0..current_gen_info.gen_list.len() {
                    image_option_counter[i] += 1;
                    if image_option_counter[i] == current_gen_info.image_options[i].len() {
                        image_option_counter[i] = 0;
                    } else {
                        break 'for_label;
                    }
                }
                break 'loop_label;
            }
        }
        return None;
    }
}

#[cfg(test)]
mod homomorphism_tests {
    use crate::{
        composition_table::group::examples,
        free_group::todd_coxeter::FinitelyGeneratedGroupPresentation,
    };

    use super::*;

    #[test]
    fn homomorphism_state() {
        {
            //identity map
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 1, 2, 3, 4, 5],
            };
            if f.check_state().is_err() {
                panic!()
            }
        }

        {
            //injective map
            let grp_g = examples::cyclic_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 2, 4],
            };
            if f.check_state().is_err() {
                panic!()
            }
        }

        {
            //surjective map
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(3);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 1, 2, 0, 1, 2],
            };
            if f.check_state().is_err() {
                panic!()
            }
        }

        {
            //bad func
            let grp_g = examples::cyclic_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 1, 2, 3],
            };
            if let Ok(()) = f.check_state() {
                panic!()
            }
        }

        {
            //bad func
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(3);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 1, 2, 3, 4, 5, 6],
            };
            if let Ok(()) = f.check_state() {
                panic!()
            }
        }

        {
            //bad hom
            let grp_g = examples::cyclic_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 1, 2],
            };
            if let Ok(()) = f.check_state() {
                panic!()
            }
        }
    }

    #[test]
    fn isomorphism_state() {
        {
            //happy
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Isomorphism {
                left_group: &grp_g,
                right_group: &grp_h,
                left_func: vec![0, 5, 4, 3, 2, 1],
                right_func: vec![0, 5, 4, 3, 2, 1],
            };
            if f.check_state().is_err() {
                panic!()
            }
        }

        {
            //size mismatch
            let grp_g = examples::cyclic_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Isomorphism {
                left_group: &grp_g,
                right_group: &grp_h,
                left_func: vec![0, 2, 4],
                right_func: vec![0, 1, 2, 0, 1, 2],
            };
            if let Ok(()) = f.check_state() {
                panic!()
            }
        }

        {
            //not inv
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Isomorphism {
                left_group: &grp_g,
                right_group: &grp_h,
                left_func: vec![0, 2, 4, 0, 2, 4],
                right_func: vec![0, 5, 4, 3, 2, 1],
            };
            if let Ok(()) = f.check_state() {
                panic!()
            }
        }
    }

    #[test]
    fn homomorphism_to_isomorphism() {
        {
            //happy
            let grp_g = examples::cyclic_group_structure(7);
            let grp_h = examples::cyclic_group_structure(7);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 3, 6, 2, 5, 1, 4],
            };
            f.check_state().unwrap();
            if let Some(_f_iso) = f.to_isomorphism() {
            } else {
                panic!()
            }
        }

        {
            //size mismatch
            let grp_g = examples::cyclic_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 2, 4],
            };
            f.check_state().unwrap();
            if let Some(_f_iso) = f.to_isomorphism() {
                panic!()
            }
        }

        {
            //not surjective
            let grp_g = examples::cyclic_group_structure(6);
            let grp_h = examples::cyclic_group_structure(6);
            let f = Homomorphism {
                domain: &grp_g,
                range: &grp_h,
                func: vec![0, 2, 4, 0, 2, 4],
            };
            f.check_state().unwrap();
            if let Some(__iso) = f.to_isomorphism() {
                panic!()
            }
        }
    }

    #[test]
    fn test_find_isomorphism() {
        {
            let grp_g = examples::symmetric_group_structure(3);
            let grp_h = examples::dihedral_group_structure(3);

            if let Some(f) = find_isomorphism(&grp_g, &grp_h) {
                assert!(std::ptr::eq(&grp_g, f.left_group));
                assert!(std::ptr::eq(&grp_h, f.right_group));
            } else {
                panic!()
            }
        }

        {
            let grp_g = examples::symmetric_group_structure(3);
            let grp_h = examples::cyclic_group_structure(6);

            if let Some(_f) = find_isomorphism(&grp_g, &grp_h) {
                panic!();
            }
        }

        {
            let grp_g = examples::symmetric_group_structure(5);
            let mut grp_h = FinitelyGeneratedGroupPresentation::new();
            let a = grp_h.add_generator();
            let b = grp_h.add_generator();
            let c = grp_h.add_generator();
            grp_h.add_relation(a.pow(2));
            grp_h.add_relation(b.pow(2));
            grp_h.add_relation(c.pow(2));
            grp_h.add_relation((&a * &c).pow(2));
            grp_h.add_relation((&a * &b).pow(3));
            grp_h.add_relation((&b * &c).pow(5));
            let grp_h = grp_h.into_finite_group(); //A5 x C2

            if let Some(_f) = find_isomorphism(&grp_g, &grp_h) {
                panic!()
            }
        }
    }
}
