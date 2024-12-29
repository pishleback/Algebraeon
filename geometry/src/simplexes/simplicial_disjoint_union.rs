use std::collections::{HashMap, HashSet};

use super::*;

#[derive(Clone)]
pub struct LabelledSimplicialDisjointUnion<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> where
    FS::Set: Hash,
{
    ambient_space: SP,
    simplexes: HashMap<Simplex<FS, SP>, T>,
}

pub type SimplicialDisjointUnion<FS, SP> = LabelledSimplicialDisjointUnion<FS, SP, ()>;

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > From<&LabelledSimplicialComplex<FS, SP, T>> for LabelledSimplicialDisjointUnion<FS, SP, T>
where
    FS::Set: Hash,
{
    fn from(sc: &LabelledSimplicialComplex<FS, SP, T>) -> Self {
        Self {
            ambient_space: sc.ambient_space(),
            simplexes: sc
                .labelled_simplexes()
                .into_iter()
                .map(|(spx, label)| (spx.clone(), label.clone()))
                .collect(),
        }
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > From<&LabelledPartialSimplicialComplex<FS, SP, T>>
    for LabelledSimplicialDisjointUnion<FS, SP, T>
where
    FS::Set: Hash,
{
    fn from(sc: &LabelledPartialSimplicialComplex<FS, SP, T>) -> Self {
        Self {
            ambient_space: sc.ambient_space(),
            simplexes: sc
                .labelled_simplexes()
                .into_iter()
                .map(|(s, label)| (s.clone(), label.clone()))
                .collect(),
        }
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledSimplexCollection<FS, SP, T> for LabelledSimplicialDisjointUnion<FS, SP, T>
where
    FS::Set: Hash,
{
    type WithLabel<S: Eq + Clone> = LabelledSimplicialDisjointUnion<FS, SP, S>;
    type SubsetType = LabelledSimplicialDisjointUnion<FS, SP, T>;

    fn new_labelled(
        ambient_space: SP,
        simplexes: HashMap<Simplex<FS, SP>, T>,
    ) -> Result<Self, &'static str> {
        //todo: check simplexes are disjoint
        Ok(Self {
            ambient_space,
            simplexes,
        })
    }

    fn new_labelled_unchecked(ambient_space: SP, simplexes: HashMap<Simplex<FS, SP>, T>) -> Self {
        Self {
            ambient_space,
            simplexes,
        }
    }

    fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    fn labelled_simplexes(&self) -> HashMap<&Simplex<FS, SP>, &T> {
        self.simplexes.iter().collect()
    }

    fn into_labelled_simplexes(self) -> HashMap<Simplex<FS, SP>, T> {
        self.simplexes
    }

    fn into_partial_simplicial_complex(self) -> LabelledPartialSimplicialComplex<FS, SP, T> {
        self.refine_to_partial_simplicial_complex()
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledSimplicialDisjointUnion<FS, SP, T>
where
    FS::Set: Hash,
{
    pub(super) fn check(&self) {
        for (spx_a, _label_a) in &self.simplexes {
            for (spx_b, _label_b) in &self.simplexes {
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

    pub fn refine_to_partial_simplicial_complex(
        mut self,
    ) -> LabelledPartialSimplicialComplex<FS, SP, T> {
        let ambient_space = self.ambient_space();

        //maintain a list of pairs of simplexes which may intersect on their boundary
        let mut pairs_todo: HashMap<Simplex<FS, SP>, HashSet<Simplex<FS, SP>>> = HashMap::new();
        let simplexes = self.simplexes().into_iter().collect::<Vec<_>>();
        for i in 0..simplexes.len() {
            for j in 0..simplexes.len() {
                if i != j {
                    let spx_i = simplexes[i];
                    let spx_j = simplexes[j];
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
                    debug_assert!(self.simplexes.contains_key(&spx1));
                    debug_assert!(self.simplexes.contains_key(&spx2));

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
                                debug_assert!(self.simplexes.contains_key(&spx));
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
                                debug_assert!(self.simplexes.contains_key(&spx));
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

                            let spx1_label = self.simplexes.get(&spx1).unwrap().clone();
                            let spx2_label = self.simplexes.get(&spx2).unwrap().clone();

                            //remove spx1 and spx2
                            self.simplexes.remove(&spx1);
                            self.simplexes.remove(&spx2);

                            //add the refinements of spx1 and spx2 and update the pairs todo
                            for spx1_repl in spx1_replacement
                                .as_simplicial_complex()
                                .subset_by_label(&InteriorBoundaryLabel::Interior)
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
                                self.simplexes.insert(spx1_repl, spx1_label.clone());
                            }

                            for spx2_repl in spx2_replacement
                                .as_simplicial_complex()
                                .subset_by_label(&InteriorBoundaryLabel::Interior)
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
                                self.simplexes.insert(spx2_repl, spx2_label.clone());
                            }
                        }
                    }
                }
            }
        }

        LabelledPartialSimplicialComplex::new_labelled_unchecked(ambient_space, self.simplexes)
    }
}
