use std::{
    collections::{HashMap, HashSet},
    marker::PhantomData,
};

use itertools::Itertools;

use crate::rings::ring_structure::structure::RealToFloatStructure;

use super::*;

pub struct VennResult<W, X> {
    pub left: W,
    pub middle: X,
    pub right: W,
}

#[derive(Clone)]
pub struct SCSpxInfo<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
{
    inv_bdry: HashSet<Simplex<FS, SP>>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> std::fmt::Debug
    for SCSpxInfo<FS, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SCSpxInfo")
            .field("inv_bdry", &self.inv_bdry)
            .finish()
    }
}

#[derive(Clone)]
pub struct SimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> {
    ambient_space: SP,
    simplexes: HashMap<Simplex<FS, SP>, SCSpxInfo<FS, SP>>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> std::fmt::Debug
    for SimplicialComplex<FS, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialComplex<FS, SP>
where
    FS::Set: Hash,
{
    fn check(&self) {
        let mut inv_bdry_map = HashMap::new();
        for spx in self.simplexes.keys() {
            inv_bdry_map.insert(spx.clone(), HashSet::new());
        }

        for (spx, info) in &self.simplexes {
            for bdry_spx in spx.proper_sub_simplices_not_null() {
                assert!(self.simplexes.contains_key(&bdry_spx));
                inv_bdry_map.get_mut(&bdry_spx).unwrap().insert(spx.clone());
            }
        }

        for (spx, info) in &self.simplexes {
            assert_eq!(&info.inv_bdry, inv_bdry_map.get(spx).unwrap());
        }

        //check that every pair of distinct simplexes intersect in the empty set
        SimplicialDisjointUnion::new_unchecked(
            self.ambient_space(),
            self.simplexes.keys().map(|s| s.clone()).collect(),
        )
        .check();
    }

    pub fn new(
        ambient_space: SP,
        simplexes: HashSet<Simplex<FS, SP>>,
    ) -> Result<Self, &'static str> {
        for simplex in &simplexes {
            assert_eq!(simplex.ambient_space().borrow(), ambient_space.borrow());
            if simplex.points().len() == 0 {
                return Err("Simplicial complex musn't contain the null simplex");
            }
        }

        let mut simplexes = simplexes
            .into_iter()
            .map(|spx| {
                (
                    spx,
                    SCSpxInfo {
                        inv_bdry: HashSet::new(),
                    },
                )
            })
            .collect::<HashMap<_, _>>();

        for simplex in simplexes.keys().map(|s| s.clone()).collect::<Vec<_>>() {
            for bdry_spx in simplex.proper_sub_simplices_not_null() {
                match simplexes.get_mut(&bdry_spx) {
                    Some(entry) => {
                        entry.inv_bdry.insert(simplex.clone());
                    }
                    None => {
                        return Err("Simplicial complex must be closed under taking boundaries");
                    }
                }
            }
        }

        Ok(Self {
            ambient_space,
            simplexes,
        })
    }

    pub fn simplexes(&self) -> Vec<&Simplex<FS, SP>> {
        self.simplexes.keys().collect()
    }

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    //remove simplexes and remove them from the inverse boundary of any others
    //self may not be in a valid state after this operation
    fn remove_simplexes(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
        for spx in &simplexes {
            for bdry_spx in spx.proper_sub_simplices_not_null() {
                match self.simplexes.get_mut(&bdry_spx) {
                    Some(info) => {
                        info.inv_bdry.remove(spx);
                    }
                    None => {}
                }
            }
        }
        for spx in &simplexes {
            self.simplexes.remove(spx);
        }
    }

    //add the given simplexes and add them to the inverse boundary map on anything on their boundaries
    //must be added together to cover the case where there are mutual boundary relations
    //self may not be in a valid state after this operation
    fn add_simplexes(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
        for spx in &simplexes {
            self.simplexes.insert(
                spx.clone(),
                SCSpxInfo {
                    inv_bdry: HashSet::new(),
                },
            );
        }
        for spx in &simplexes {
            for bdry_spx in spx.proper_sub_simplices_not_null() {
                self.simplexes
                    .get_mut(&bdry_spx)
                    .unwrap()
                    .inv_bdry
                    .insert(spx.clone());
            }
        }
    }

    //refine self by replacing a simplex with a refinement
    fn refine_simplex(&mut self, spx: &Simplex<FS, SP>, replacement: &SimplicialComplex<FS, SP>) {
        let spx_points = spx.points().iter().collect::<HashSet<_>>();

        //for each simplex in the replacement calculate which boundary simplex of spx is belongs to
        //classify these into a list
        //[belongs to a point, belongs to an edge, belongs to a face, ..., belongs to a ride, belongs to a facet, belongs to the interior]
        let sub_spx_info = {
            let mut sub_spx_info = (0..spx.n())
                .map(|i| HashMap::new())
                .collect::<Vec<HashMap<Vec<usize>, Vec<Simplex<FS, SP>>>>>();
            for rep_spx in replacement.simplexes() {
                //within which intersection of facets of spx does the refined simplex live
                let mut spx_facets = (0..spx.n()).collect::<HashSet<_>>();
                for pt in rep_spx.points() {
                    for i in spx_facets.clone() {
                        match EmbeddedAffineSubspace::<_, _, Rc<_>>::new_affine_span(
                            self.ambient_space(),
                            spx.facet(i).into_points(),
                        )
                        .unwrap()
                        .0
                        .unembed_point(&pt)
                        {
                            Some(_) => {}
                            None => {
                                spx_facets.remove(&i);
                            }
                        }
                    }
                }
                let mut spx_points = vec![];
                for i in 0..spx.n() {
                    if !spx_facets.contains(&i) {
                        spx_points.push(i);
                    }
                }
                debug_assert!(spx_points.len() >= 1);
                debug_assert!(spx_points.len() <= spx.n());
                spx_points.sort_unstable();
                sub_spx_info
                    .get_mut(spx_points.len() - 1)
                    .unwrap()
                    .entry(spx_points)
                    .or_insert(vec![])
                    .push(rep_spx.clone());
            }
            sub_spx_info
        };

        for (sub_n, sub_replacements) in sub_spx_info.into_iter().enumerate() {
            let sub_n = sub_n + 1;
            for (key, sub_replacement) in &sub_replacements {
                debug_assert_eq!(key.len(), sub_n);
                let sub_spx = Simplex::new(
                    self.ambient_space(),
                    key.into_iter().map(|i| spx.point(*i).clone()).collect(),
                )
                .unwrap();

                let mut spxs_to_remove = vec![];
                let mut spxs_to_add = vec![];
                spxs_to_remove.push(sub_spx.clone());
                for sub_repl_spx in sub_replacement {
                    spxs_to_add.push(sub_repl_spx.clone());
                }
                for ext_points in self
                    .simplexes
                    .get(&sub_spx)
                    .unwrap()
                    .inv_bdry
                    .iter()
                    .map(|sub_spx_inv_bdry_spx| {
                        let mut points = sub_spx_inv_bdry_spx
                            .points()
                            .iter()
                            .map(|pt| pt.clone())
                            .collect::<HashSet<_>>();
                        for pt in sub_spx.points() {
                            points.remove(pt);
                        }
                        points
                    })
                    .filter(|points| !points.iter().all(|pt| spx_points.contains(pt)))
                    .collect::<Vec<_>>()
                {
                    //adding points to sub_spx.points gives a simplex with sub_spx on its boundary which is not some larger facet of spx with sub_spx on its boundary
                    //we shall refine this simplex by joining the refinement of sub_spx to the extended points and update the inverse boundary maps as we go

                    //old_ext_spx = Simplex(ext_points + sub_spx.points)
                    let old_ext_spx = Simplex::new(self.ambient_space(), {
                        let mut points = ext_points.clone().into_iter().collect::<Vec<_>>();
                        for pt in sub_spx.points() {
                            points.push(pt.clone());
                        }
                        points
                    })
                    .unwrap();
                    spxs_to_remove.push(old_ext_spx);
                    for sub_spx_repl in sub_replacement {
                        let new_ext_spx = Simplex::new(self.ambient_space(), {
                            let mut points = ext_points.clone().into_iter().collect::<Vec<_>>();
                            for pt in sub_spx_repl.points() {
                                points.push(pt.clone());
                            }
                            points
                        })
                        .unwrap();
                        spxs_to_add.push(new_ext_spx);
                    }
                }

                debug_assert_eq!(
                    spxs_to_add.len(),
                    spxs_to_add.iter().collect::<HashSet<_>>().len()
                );

                if spxs_to_remove.len() == spxs_to_add.len() {
                    debug_assert_eq!(
                        spxs_to_remove.iter().collect::<HashSet<_>>(),
                        spxs_to_add.iter().collect::<HashSet<_>>()
                    );
                }

                self.remove_simplexes(spxs_to_remove);
                self.add_simplexes(spxs_to_add);
            }
        }

        #[cfg(debug_assertions)]
        self.check();
    }

    /*
    pub fn venn(&self, other: &Self) -> (Self, Self)
    where
        FS: RealToFloatStructure,
    {
        let mut a = self.clone();
        let mut b = other.clone();
        println!("sc venn");
        // println!("a = {:?}", a);
        // println!("b = {:?}", b);

        'MAIN_LOOP: loop {
            'PAIR_SEARCH: {
                let mut spx_pairs = a
                    .simplexes()
                    .into_iter()
                    .map(|s| s.clone())
                    .cartesian_product(b.simplexes().into_iter().map(|s| s.clone()))
                    .collect::<Vec<_>>();

                spx_pairs.sort_unstable_by_key(|(x, y)| x.n() + y.n());

                for (a_spx, b_spx) in spx_pairs.into_iter().rev() {
                    let a_bdry = a_spx
                        .sub_simplices_not_null()
                        .into_iter()
                        .collect::<HashSet<_>>();
                    let b_bdry = b_spx
                        .sub_simplices_not_null()
                        .into_iter()
                        .collect::<HashSet<_>>();

                    if a_bdry.contains(&b_spx) || b_bdry.contains(&a_spx) {
                        // equal
                    } else {
                        let a_spx_ch = ConvexHull::from_simplex(a_spx.clone());
                        let b_spx_ch = ConvexHull::from_simplex(b_spx.clone());
                        let overlap = ConvexHull::intersect(&a_spx_ch, &b_spx_ch);
                        if overlap.affine_span_dimension() < a_spx.n()
                            && overlap.affine_span_dimension() < b_spx.n()
                        {
                            // disjoint
                        } else {
                            println!("scary {} {} {} {}", a_spx.n(), b_spx.n(), a.simplexes.len(), b.simplexes.len());

                            // if b.simplexes.len() > 300 {
                            //     <crate::drawing::canvas2d::Diagram2dCanvas as crate::drawing::Canvas>::run(
                            //         |canvas| {
                            //             canvas.draw(
                            //                 &a,
                            //                 (1.0, 0.0, 0.0),
                            //             );
                            //             canvas.draw(
                            //                 &b,
                            //                 (0.0, 1.0, 0.0),
                            //             );
                            //         },
                            //     );
                            // }

                            println!("    {:?}", a_spx);
                            println!("    {:?}", b_spx);

                            let mut refined_a_spx = overlap.clone();
                            let mut refined_b_spx = overlap.clone();
                            for pt in a_spx.points() {
                                refined_a_spx.extend_by_point(pt.clone());
                            }
                            for pt in b_spx.points() {
                                refined_b_spx.extend_by_point(pt.clone());
                            }

                            a.refine_simplex(&a_spx, &refined_a_spx.as_simplicial_complex().entire);
                            b.refine_simplex(&b_spx, &refined_b_spx.as_simplicial_complex().entire);

                            break 'PAIR_SEARCH;

                            // todo!();
                        }
                    }
                }
                break 'MAIN_LOOP;
            }
        }
        println!("Done :OOOO");

        (a, b)
    }
    */
}

#[derive(Clone)]
pub struct PartialSimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> {
    ambient_space: SP,
    simplexes: HashSet<Simplex<FS, SP>>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> std::fmt::Debug
    for PartialSimplicialComplex<FS, SP>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("PartialSimplicialComplex")
            .field("simplexes", &self.simplexes)
            .finish()
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    PartialSimplicialComplex<FS, SP>
where
    FS::Set: Hash,
{
    pub fn new(
        ambient_space: SP,
        simplexes: HashSet<Simplex<FS, SP>>,
    ) -> Result<Self, &'static str> {
        Ok(Self {
            ambient_space,
            simplexes,
        })
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
}

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

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        SC: Borrow<SimplicialComplex<FS, SP>>,
    > From<SC> for SimplicialDisjointUnion<FS, SP>
where
    FS::Set: Hash,
{
    fn from(sc: SC) -> Self {
        let sc = sc.borrow();
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
    fn check(&self) {
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
                    assert!(
                        overlap.affine_span_dimension() < spx_a.n()
                            && overlap.affine_span_dimension() < spx_b.n()
                    );
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

    pub fn subtract(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        Self::new_unchecked(ambient_space.clone(), {
            let mut simplexes = HashSet::new();
            for self_spx in &self.simplexes {
                let mut self_leftover = HashSet::from([self_spx.clone()]);
                for other_spx in &other.simplexes {
                    self_leftover = self_leftover
                        .into_iter()
                        .map(|self_leftover_spx| {
                            Simplex::venn(&self_leftover_spx, other_spx)
                                .left
                                .into_simplexes()
                        })
                        .flatten()
                        .collect();
                }
                for spx in self_leftover {
                    simplexes.insert(spx);
                }
            }
            simplexes
        })
    }

    pub fn venn(
        &self,
        other: &Self,
    ) -> VennResult<SimplicialDisjointUnion<FS, SP>, SimplicialDisjointUnion<FS, SP>> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        VennResult {
            left: Self::subtract(self, other),
            middle: Self::new_unchecked(ambient_space.clone(), {
                let mut simplexes = HashSet::new();
                for self_spx in &self.simplexes {
                    for other_spx in &other.simplexes {
                        for spx in Simplex::venn(self_spx, other_spx).middle.into_simplexes() {
                            simplexes.insert(spx);
                        }
                    }
                }
                simplexes
            }),
            right: Self::subtract(other, self),
        }
    }

    pub fn intersect(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        Self::venn(self, other).middle
    }

    pub fn union(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let VennResult {
            left,
            middle,
            right,
        } = Self::venn(&self, other);
        let mut simplexes = HashSet::new();
        for spx in left.into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in middle.into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in right.into_simplexes() {
            simplexes.insert(spx);
        }
        return Self::new_unchecked(ambient_space, simplexes);
    }

    pub fn refine_to_simplicial_complex(mut self) -> PartialSimplicialComplex<FS, SP> {
        let ambient_space = self.ambient_space();

        //maintain a list of pairs of simplexes which may intersect on their boundary
        let mut pairs_todo: HashMap<Simplex<FS, SP>, HashSet<Simplex<FS, SP>>> = HashMap::new();
        let simplexes = self
            .simplexes()
            .iter()
            .map(|spx| spx.clone())
            .collect::<Vec<_>>();
        for i in (0..simplexes.len()) {
            for j in (0..simplexes.len()) {
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
                            println!("sad {:?} {:?}", spx1, spx2);

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
                                .interior
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
                                .interior
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
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Simplex<FS, SP>
where
    FS::Set: Hash,
{
    pub fn venn(
        &self,
        other: &Self,
    ) -> VennResult<SimplicialDisjointUnion<FS, SP>, SimplicialDisjointUnion<FS, SP>> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        let overlap = ConvexHull::intersect(
            &ConvexHull::from_simplex(self.clone()),
            &ConvexHull::from_simplex(other.clone()),
        );

        let mut self_ext = overlap.clone();
        for pt in self.points() {
            self_ext.extend_by_point(pt.clone());
        }
        let self_parts = self_ext
            .as_simplicial_complex()
            .interior
            .simplexes()
            .into_iter()
            .map(|s| s.clone())
            .collect::<HashSet<_>>();

        let mut other_ext = overlap.clone();
        for pt in other.points() {
            other_ext.extend_by_point(pt.clone());
        }
        let other_parts = other_ext
            .as_simplicial_complex()
            .interior
            .simplexes()
            .into_iter()
            .map(|s| s.clone())
            .collect::<HashSet<_>>();

        let mut left = self_parts.clone();
        let mut middle = HashSet::new();
        let mut right = HashSet::new();
        for other_part in other_parts {
            if left.contains(&other_part) {
                left.remove(&other_part);
                middle.insert(other_part);
            } else {
                right.insert(other_part);
            }
        }

        VennResult {
            left: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), left),
            middle: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), middle),
            right: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), right),
        }
    }
}
