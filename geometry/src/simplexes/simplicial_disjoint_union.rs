use std::collections::{HashMap, HashSet};

use super::*;

#[derive(Clone)]
pub struct SimplicialDisjointUnion<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> where
    FS::Set: Hash,
{
    ambient_space: SP,
    simplexes: HashSet<Simplex<FS, SP>>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    From<&SimplicialComplex<FS, SP>> for SimplicialDisjointUnion<FS, SP>
where
    FS::Set: Hash,
{
    fn from(sc: &SimplicialComplex<FS, SP>) -> Self {
        Self {
            ambient_space: sc.ambient_space(),
            simplexes: sc.simplexes().into_iter().map(|s| s.clone()).collect(),
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    From<&PartialSimplicialComplex<FS, SP>> for SimplicialDisjointUnion<FS, SP>
where
    FS::Set: Hash,
{
    fn from(sc: &PartialSimplicialComplex<FS, SP>) -> Self {
        Self {
            ambient_space: sc.ambient_space(),
            simplexes: sc.simplexes().into_iter().map(|s| s.clone()).collect(),
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialDisjointUnion<FS, SP>
where
    FS::Set: Hash,
{
    pub(super) fn check(&self) {
        for spx_a in &self.simplexes {
            for spx_b in &self.simplexes {
                let bdry_a = spx_a
                    .sub_simplices_not_null()
                    .into_iter()
                    .collect::<HashSet<_>>();
                let bdry_b = spx_b
                    .sub_simplices_not_null()
                    .into_iter()
                    .collect::<HashSet<_>>();
                if !bdry_a.contains(spx_b) && !bdry_b.contains(spx_a) {
                    let overlap = ConvexHull::intersect(
                        &ConvexHull::from_simplex(spx_a.clone()),
                        &ConvexHull::from_simplex(spx_b.clone()),
                    );

                    if !(overlap.affine_span_dimension() < spx_a.n()
                        && overlap.affine_span_dimension() < spx_b.n())
                    {
                        println!("spx_a = {:?}", spx_a);
                        println!("spx_b = {:?}", spx_b);
                        panic!("simplicial complex simplex overlap");
                    }
                }
            }
        }
    }

    pub fn new_unchecked(ambient_space: SP, simplexes: HashSet<Simplex<FS, SP>>) -> Self {
        Self {
            ambient_space,
            simplexes,
        }
    }

    pub fn simplexes(&self) -> &HashSet<Simplex<FS, SP>> {
        &self.simplexes
    }

    pub fn into_simplexes(self) -> HashSet<Simplex<FS, SP>> {
        self.simplexes
    }

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    pub fn refine_to_partial_simplicial_complex(mut self) -> PartialSimplicialComplex<FS, SP> {
        let ambient_space = self.ambient_space();

        //maintain a list of pairs of simplexes which may intersect on their boundary
        let mut pairs_todo: HashMap<Simplex<FS, SP>, HashSet<Simplex<FS, SP>>> = HashMap::new();
        let simplexes = self
            .simplexes()
            .iter()
            .map(|spx| spx.clone())
            .collect::<Vec<_>>();
        for i in 0..simplexes.len() {
            for j in 0..simplexes.len() {
                if i != j {
                    let spx_i = &simplexes[i];
                    let spx_j = &simplexes[j];
                    pairs_todo
                        .entry(spx_i.clone())
                        .or_insert(HashSet::new())
                        .insert(spx_j.clone());
                }
            }
        }

        while !pairs_todo.is_empty() {
            let spx1 = pairs_todo.keys().into_iter().next().unwrap().clone();
            match pairs_todo.get(&spx1).unwrap().iter().next().cloned() {
                None => {
                    pairs_todo.remove(&spx1);
                }
                Some(spx2) => {
                    //The pair (spx1, spx2) no longer needs to be checked
                    pairs_todo.get_mut(&spx1).unwrap().remove(&spx2);
                    pairs_todo.get_mut(&spx2).unwrap().remove(&spx1);

                    debug_assert_ne!(spx1, spx2);
                    debug_assert!(self.simplexes.contains(&spx1));
                    debug_assert!(self.simplexes.contains(&spx2));

                    let overlap = ConvexHull::intersect(
                        &ConvexHull::from_simplex(spx1.clone()),
                        &ConvexHull::from_simplex(spx2.clone()),
                    );

                    if !overlap.is_empty() {
                        if match Simplex::new(
                            ambient_space.clone(),
                            overlap.defining_points().into_iter().collect(),
                        ) {
                            Ok(overlap_spx) => {
                                let spx1_points = spx1.points().into_iter().collect::<HashSet<_>>();
                                let spx2_points = spx2.points().into_iter().collect::<HashSet<_>>();
                                !overlap_spx
                                    .points()
                                    .into_iter()
                                    .all(|pt| spx1_points.contains(pt) && spx2_points.contains(pt))
                            }
                            Err(_) => true,
                        } {
                            //there is a bad overlap between spx1 and spx2
                            let mut spx1_replacement = overlap.clone();
                            for pt in spx1.points() {
                                spx1_replacement.extend_by_point(pt.clone());
                            }
                            let mut spx2_replacement = overlap.clone();
                            for pt in spx2.points() {
                                spx2_replacement.extend_by_point(pt.clone());
                            }

                            //remove any pairs containing spx1 or spx2
                            let mut spx1_paired = vec![];
                            for spx in pairs_todo.get(&spx1).unwrap_or(&HashSet::new()).clone() {
                                debug_assert_ne!(spx, spx1);
                                debug_assert_ne!(spx, spx2);
                                debug_assert!(self.simplexes.contains(&spx));
                                pairs_todo
                                    .get_mut(&spx)
                                    .unwrap_or(&mut HashSet::new())
                                    .remove(&spx1);
                                spx1_paired.push(spx);
                            }
                            pairs_todo.remove(&spx1);

                            let mut spx2_paired = vec![];
                            for spx in pairs_todo.get(&spx2).unwrap_or(&HashSet::new()).clone() {
                                debug_assert_ne!(spx, spx1);
                                debug_assert_ne!(spx, spx2);
                                debug_assert!(self.simplexes.contains(&spx));
                                pairs_todo
                                    .get_mut(&spx)
                                    .unwrap_or(&mut HashSet::new())
                                    .remove(&spx2);
                                spx2_paired.push(spx);
                            }
                            pairs_todo.remove(&spx2);

                            // //pairs should now be in a valid state again
                            // for (a, bs) in &pairs_todo {
                            //     debug_assert!(self.simplexes.contains(a));
                            //     for b in bs {
                            //         debug_assert!(self.simplexes.contains(b));
                            //     }
                            // }

                            //remove spx1 and spx2
                            self.simplexes.remove(&spx1);
                            self.simplexes.remove(&spx2);

                            //add the refinements of spx1 and spx2 and update the pairs todo
                            for spx1_repl in spx1_replacement
                                .as_simplicial_complex()
                                .labelled_subset(&InteriorBoundaryLabel::Interior)
                                .into_simplexes()
                            {
                                for spx in &spx1_paired {
                                    pairs_todo
                                        .entry(spx1_repl.clone())
                                        .or_insert(HashSet::new())
                                        .insert(spx.clone());
                                    pairs_todo
                                        .entry(spx.clone())
                                        .or_insert(HashSet::new())
                                        .insert(spx1_repl.clone());
                                }
                                self.simplexes.insert(spx1_repl);
                            }

                            for spx2_repl in spx2_replacement
                                .as_simplicial_complex()
                                .labelled_subset(&InteriorBoundaryLabel::Interior)
                                .into_simplexes()
                            {
                                for spx in &spx2_paired {
                                    pairs_todo
                                        .entry(spx2_repl.clone())
                                        .or_insert(HashSet::new())
                                        .insert(spx.clone());
                                    pairs_todo
                                        .entry(spx.clone())
                                        .or_insert(HashSet::new())
                                        .insert(spx2_repl.clone());
                                }
                                self.simplexes.insert(spx2_repl);
                            }
                        }
                    }
                }
            }
        }

        PartialSimplicialComplex::new(ambient_space, self.into_simplexes()).unwrap()
    }

    // pub fn closure_as_simplicial_complex(self) -> SimplicialComplex<FS, SP> {
    //     self.refine_to_partial_simplicial_complex()
    //         .closure_as_simplicial_complex()
    // }
}
