use std::collections::{HashMap, HashSet};

use super::*;

#[derive(Clone)]
pub struct SCSpxInfo<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> {
    inv_bdry: HashSet<Simplex<FS, SP>>,
    label: T,
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > std::fmt::Debug for SCSpxInfo<FS, SP, T>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("SCSpxInfo")
            .field("inv_bdry", &self.inv_bdry)
            .finish()
    }
}

#[derive(Clone)]
pub struct LabelledSimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    T: Eq + Clone,
> {
    ambient_space: SP,
    simplexes: HashMap<Simplex<FS, SP>, SCSpxInfo<FS, SP, T>>,
}

pub type SimplicialComplex<FS, SP> = LabelledSimplicialComplex<FS, SP, ()>;

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > std::fmt::Debug for LabelledSimplicialComplex<FS, SP, T>
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
    pub fn new(
        ambient_space: SP,
        simplexes: HashSet<Simplex<FS, SP>>,
    ) -> Result<Self, &'static str> {
        Self::new_labelled(
            ambient_space,
            simplexes.into_iter().map(|spx| (spx, ())).collect(),
        )
    }

    pub fn interior_and_boundary(
        &self,
    ) -> LabelledSimplicialComplex<FS, SP, InteriorBoundaryLabel> {
        /*
        let n be the dimension of the space self is living in
         - every simplex of rank n is part of the interior
         - a simplex of rank n-1 is the facet of at most 2 simplices of rank n, and is part of the interior if and only if it is the facet of exactly 2 simplices of rank n
         - a simplex of rank less or equal to n-2 is part of the interior iff it is in the boundary of some strictly higher rank simplex AND every strictly higher rank simplex containing it as part of the boundary is part of the interior
        */

        let n = self.ambient_space().borrow().affine_dimension();

        let mut simplexes = HashMap::new();

        let mut all = self.simplexes.keys().cloned().collect::<Vec<_>>();
        all.sort_unstable_by_key(|s| std::cmp::Reverse(s.n())); //so that we process largest rank first
        for simplex in all {
            let r = simplex.n();
            if r == n {
                //rank n simplex is always part of the interior
                simplexes.insert(simplex, InteriorBoundaryLabel::Interior);
            } else {
                let inv_bdry = &self.simplexes.get(&simplex).unwrap().inv_bdry;
                if r == n - 1 {
                    //rank n-1 simplex is part of the boundary of at most 2 rank n simplices
                    //it is part of the boundary of the simplicial complex iff it is the boundary of exactly 2
                    match inv_bdry.len() {
                        0 | 1 => {
                            simplexes.insert(simplex, InteriorBoundaryLabel::Boundary);
                        }
                        2 => {
                            simplexes.insert(simplex, InteriorBoundaryLabel::Interior);
                        }
                        _ => panic!(
                        "rank n-1 simplex should be in the boundary of at most 2 rank n simplices"
                    ),
                    }
                } else {
                    //rank < n-1 simplex is part of the interior iff it is part of the boundary of at least one simplex and every such simplex is part of the interior
                    debug_assert!(r < n - 1);
                    if inv_bdry.is_empty() {
                        simplexes.insert(simplex, InteriorBoundaryLabel::Boundary);
                    } else {
                        if inv_bdry
                            .iter()
                            .all(|b| simplexes.get(b).unwrap() == &InteriorBoundaryLabel::Interior)
                        {
                            simplexes.insert(simplex, InteriorBoundaryLabel::Interior);
                        } else {
                            simplexes.insert(simplex, InteriorBoundaryLabel::Boundary);
                        }
                    }
                }
            }
        }

        LabelledSimplicialComplex::new_labelled(self.ambient_space(), simplexes).unwrap()
    }

    pub fn interior(&self) -> PartialSimplicialComplex<FS, SP> {
        self.interior_and_boundary()
            .labelled_subset(&InteriorBoundaryLabel::Interior)
    }

    pub fn boundary(&self) -> SimplicialComplex<FS, SP> {
        self.interior_and_boundary()
            .labelled_subset(&InteriorBoundaryLabel::Boundary)
            .try_as_simplicial_complex()
            .unwrap()
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledSimplicialComplex<FS, SP, T>
where
    FS::Set: Hash,
{
    fn check(&self) {
        let mut inv_bdry_map = HashMap::new();
        for spx in self.simplexes.keys() {
            inv_bdry_map.insert(spx.clone(), HashSet::new());
        }

        for (spx, _info) in &self.simplexes {
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

    pub fn new_labelled(
        ambient_space: SP,
        simplexes: HashMap<Simplex<FS, SP>, T>,
    ) -> Result<Self, &'static str> {
        for simplex in simplexes.keys() {
            assert_eq!(simplex.ambient_space().borrow(), ambient_space.borrow());
            if simplex.points().len() == 0 {
                return Err("Simplicial complex musn't contain the null simplex");
            }
        }

        let mut simplexes = simplexes
            .into_iter()
            .map(|(spx, label)| {
                (
                    spx,
                    SCSpxInfo {
                        inv_bdry: HashSet::new(),
                        label,
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

    pub fn common_label<'a>(
        &'a self,
        simplexes: impl Iterator<Item = &'a Simplex<FS, SP>>,
    ) -> Option<&'a T> {
        let mut label = None;
        for spx in simplexes {
            let spx_label = &self.simplexes.get(&spx).unwrap().label;
            match label {
                Some(label) => {
                    if label != spx_label {
                        return None;
                    }
                }
                None => {
                    label = Some(spx_label);
                }
            }
        }
        label
    }

    pub fn labelled_subset(&self, label: &T) -> PartialSimplicialComplex<FS, SP> {
        PartialSimplicialComplex::new_unchecked(
            self.ambient_space().clone(),
            self.simplexes
                .iter()
                .filter(|(_, info)| &info.label == label)
                .map(|(spx, _)| spx.clone())
                .collect(),
        )
    }

    pub fn forget_labels(self) -> SimplicialComplex<FS, SP> {
        SimplicialComplex {
            ambient_space: self.ambient_space,
            simplexes: self
                .simplexes
                .into_iter()
                .map(|(spx, info)| {
                    (
                        spx,
                        SCSpxInfo {
                            inv_bdry: info.inv_bdry,
                            label: (),
                        },
                    )
                })
                .collect(),
        }
    }
}

/*
Input:
    A list of oriented simplicies which join to form a closed region of space
    with negative side inside and positive side outside which join

Output:
    Figure out whether the region can be filled in by fanning out from some point
    If it can't: return None
    If it can: return the simplicies to use to fill in the interior region

*/
fn simplify_in_region<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    space: SP,
    boundary_facets: Vec<OrientedSimplex<FS, SP>>,
) -> Option<Vec<Simplex<FS, SP>>>
where
    FS::Set: Hash,
{
    for spx in &boundary_facets {
        debug_assert_eq!(spx.ambient_space().borrow(), space.borrow());
    }

    let mut boundary_points: HashMap<Vector<FS, SP>, Vec<usize>> = HashMap::new();
    for (idx, spx) in boundary_facets.iter().enumerate() {
        for pt in spx.simplex().points() {
            if boundary_points.contains_key(pt) {
                boundary_points.get_mut(pt).unwrap().push(idx);
            } else {
                boundary_points.insert(pt.clone(), vec![idx]);
            }
        }
    }

    for (boundary_point, adjacent_facets) in boundary_points {
        let mut nonadjacent_facets = (0..boundary_facets.len()).collect::<HashSet<_>>();
        for idx in &adjacent_facets {
            nonadjacent_facets.remove(idx);
        }

        if nonadjacent_facets.iter().all(|idx| {
            match boundary_facets[*idx].classify_point(&boundary_point) {
                OrientationSide::Positive => false,
                OrientationSide::Neutral => false,
                OrientationSide::Negative => true,
            }
        }) {
            let mut nonadjacent_simplexes = HashSet::new();
            for idx in nonadjacent_facets {
                for spx in boundary_facets[idx].simplex().sub_simplices_not_null() {
                    nonadjacent_simplexes.insert(spx);
                }
            }
            for idx in adjacent_facets {
                for spx in boundary_facets[idx].simplex().sub_simplices_not_null() {
                    nonadjacent_simplexes.remove(&spx);
                }
            }

            let filler = nonadjacent_simplexes
                .into_iter()
                .map(|spx| {
                    let mut points = spx.points().clone();
                    points.push(boundary_point.clone());
                    Simplex::new(spx.ambient_space().clone(), points).unwrap()
                })
                .collect();

            return Some(filler);
        }
    }
    None
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        T: Eq + Clone,
    > LabelledSimplicialComplex<FS, SP, T>
where
    FS::Set: Hash,
{
    pub fn simplify(mut self) -> Self {
        //go through each point
        //compute its star
        //enter the affine subspace spanned by its star
        //compute its link as oriented simplices
        //if point is part of the link then no simplification can be made, so move on
        //check whether any point of the link is in the interior with respect to every boundary of the link
        //if such a point exists, it can be used to fill in the star in a simpler way by fanning out

        let mut pts_todo = HashSet::new();
        for simplex in self.simplexes.keys() {
            for pt in simplex.points() {
                pts_todo.insert(pt.clone());
            }
        }

        while !pts_todo.is_empty() {
            let pt = {
                let mut pts_todo_iter = pts_todo.into_iter();
                let pt = pts_todo_iter.next().unwrap();
                pts_todo = pts_todo_iter.collect();
                pt
            };

            let pt_spx = Simplex::new(self.ambient_space(), vec![pt.clone()]).unwrap();
            let (star, link) = {
                let mut star = self.simplexes.get(&pt_spx).unwrap().inv_bdry.clone();
                star.insert(pt_spx.clone());

                let mut nbd = HashSet::new();
                for spx in &star {
                    for bdry in spx.sub_simplices_not_null() {
                        nbd.insert(bdry);
                    }
                }

                let mut link = nbd.clone();
                for spx in &star {
                    link.remove(spx);
                }

                debug_assert_eq!(star.len() + link.len(), nbd.len());
                for link_spx in &link {
                    debug_assert!(self.simplexes.contains_key(&link_spx));
                }

                (star, link)
            };

            let link_points = {
                let mut link_points: Vec<Vector<FS, SP>> = vec![];
                for spx in &link {
                    for p in spx.points() {
                        link_points.push(p.clone());
                    }
                }
                link_points
            };

            let nbd_points = {
                let mut nbd_points = link_points.clone();
                nbd_points.push(pt.clone());
                nbd_points
            };

            let nbd_affine_subspace = EmbeddedAffineSubspace::new_affine_span_linearly_dependent(
                self.ambient_space(),
                nbd_points.iter().collect(),
            );

            let nbd_points_img = nbd_points
                .iter()
                .map(|pt| nbd_affine_subspace.unembed_point(pt).unwrap())
                .collect::<Vec<_>>();
            let pt_img = nbd_affine_subspace.unembed_point(&pt).unwrap();
            let pt_img_spx =
                Simplex::new(nbd_affine_subspace.embedded_space(), vec![pt_img.clone()]).unwrap();
            let star_img = star
                .iter()
                .map(|s| nbd_affine_subspace.unembed_simplex(s).unwrap())
                .collect::<HashSet<_>>();
            let link_img = link
                .iter()
                .map(|s| nbd_affine_subspace.unembed_simplex(s).unwrap())
                .collect::<HashSet<_>>();

            let nbd = LabelledSimplicialComplex::new(nbd_affine_subspace.embedded_space(), {
                let mut simplexes = HashSet::new();
                simplexes.extend(star_img.clone());
                simplexes.extend(link_img.clone());
                simplexes
            })
            .unwrap();
            let nbd_interior = nbd.interior().into_simplexes();

            if !nbd_interior.contains(&pt_img_spx) {
                /*
                pt is on the boundary of nbd and nbd looks something like this

                    l     l      l
                     +----------+
                    / \\__       \
                l  /   \  \__ s   \ l
                  /   s \    \__   \
                 /       \      \__ \
                +---------+---------+
                r    b    p    b    r

                where
                    p = point
                    b = boundary = intersection of star and ndb.boundary
                    r = rim = closure of boundary minus boundary
                    l = semilink = link minus rim
                    s = semistar = star minus boundary

                so the idea here is
                    simplify the boundary by filling in the rim to replace the boundary
                    then simplify the rest by filling in the new boundary and the link to replace the star
                */

                let boundary_img = star_img
                    .iter()
                    .filter(|spx| !nbd_interior.contains(&spx))
                    .collect::<HashSet<_>>();

                let boundary = boundary_img
                    .iter()
                    .map(|spx| nbd_affine_subspace.embed_simplex(spx))
                    .collect::<Vec<_>>();
                let semistar = {
                    let mut semistar = star.clone();
                    for spx in &boundary {
                        semistar.remove(spx);
                    }
                    semistar
                };
                debug_assert_eq!(boundary.len() + semistar.len(), star.len());
                if let (Some(boundary_label), Some(semistar_label)) = (
                    self.common_label(boundary.iter()).cloned(),
                    self.common_label(semistar.iter()).cloned(),
                ) {
                    let mut boundary_img_points = HashSet::new();
                    for spx in &boundary_img {
                        for p in spx.points() {
                            boundary_img_points.insert(p);
                        }
                    }
                    let nbd_boundary_affine_subspace =
                        EmbeddedAffineSubspace::new_affine_span_linearly_dependent(
                            nbd_affine_subspace.embedded_space(),
                            boundary_img_points.into_iter().collect(),
                        );
                    debug_assert!(
                        nbd_boundary_affine_subspace
                            .embedded_space()
                            .affine_dimension()
                            <= nbd_affine_subspace.embedded_space().affine_dimension()
                    );
                    if nbd_boundary_affine_subspace
                        .embedded_space()
                        .affine_dimension()
                        + 1
                        == nbd_affine_subspace.embedded_space().affine_dimension()
                    {
                        let ref_point_img = {
                            let mut ref_point_img = None;
                            for pt in &nbd_points_img {
                                if nbd_boundary_affine_subspace.unembed_point(pt).is_none() {
                                    ref_point_img = Some(pt.clone());
                                    break;
                                }
                            }
                            ref_point_img.unwrap()
                        };
                        let oriented_hyperplane = OrientedSimplex::new_with_negative_point(
                            nbd_boundary_affine_subspace.ambient_space(),
                            nbd_boundary_affine_subspace.get_embedding_points().clone(),
                            &ref_point_img,
                        )
                        .unwrap();
                        for pt in &nbd_points_img {
                            debug_assert!(
                                oriented_hyperplane.classify_point(pt) != OrientationSide::Positive
                            );
                        }

                        let rim_img = {
                            let mut rim_img = HashSet::new();
                            for spx in &boundary_img {
                                for bspx in spx.sub_simplices_not_null() {
                                    rim_img.insert(bspx);
                                }
                            }
                            for spx in &boundary_img {
                                rim_img.remove(spx);
                            }
                            rim_img
                        };

                        let pt_img_img =
                            nbd_boundary_affine_subspace.unembed_point(&pt_img).unwrap();

                        let rim_img_img = rim_img
                            .iter()
                            .map(|spx| {
                                OrientedSimplex::new_with_negative_point(
                                    nbd_boundary_affine_subspace.embedded_space(),
                                    nbd_boundary_affine_subspace
                                        .unembed_simplex(spx)
                                        .unwrap()
                                        .points()
                                        .clone(),
                                    &pt_img_img,
                                )
                                .unwrap()
                            })
                            .collect::<Vec<_>>();

                        if let Some(new_boundary_img_img) = simplify_in_region(
                            nbd_boundary_affine_subspace.embedded_space(),
                            rim_img_img,
                        ) {
                            let new_boundary_img = new_boundary_img_img
                                .iter()
                                .map(|spx| nbd_boundary_affine_subspace.embed_simplex(spx))
                                .collect::<Vec<_>>();

                            let sphere_img = {
                                let mut sphere_img = vec![];
                                for spx in &new_boundary_img {
                                    if spx.n() + 1
                                        == nbd_affine_subspace.embedded_space().affine_dimension()
                                    {
                                        sphere_img.push(
                                            OrientedSimplex::new_with_negative_point(
                                                nbd_affine_subspace.embedded_space(),
                                                spx.points().clone(),
                                                &ref_point_img,
                                            )
                                            .unwrap(),
                                        );
                                    }
                                }
                                for spx in &link_img {
                                    if spx.n() + 1
                                        == nbd_affine_subspace.embedded_space().affine_dimension()
                                    {
                                        sphere_img.push(
                                            OrientedSimplex::new_with_negative_point(
                                                nbd_affine_subspace.embedded_space(),
                                                spx.points().clone(),
                                                &pt_img,
                                            )
                                            .unwrap(),
                                        );
                                    }
                                }
                                sphere_img
                            };

                            if let Some(new_star_img) =
                                simplify_in_region(nbd_affine_subspace.embedded_space(), sphere_img)
                            {
                                self.remove_simplexes_unchecked(star.into_iter().collect());
                                self.add_simplexes_unchecked(
                                    new_boundary_img
                                        .into_iter()
                                        .map(|spx_img| nbd_affine_subspace.embed_simplex(&spx_img))
                                        .collect(),
                                    &boundary_label,
                                );
                                self.add_simplexes_unchecked(
                                    new_star_img
                                        .into_iter()
                                        .map(|spx_img| nbd_affine_subspace.embed_simplex(&spx_img))
                                        .collect(),
                                    &semistar_label,
                                );
                                pts_todo.extend(link_points);
                            }
                        }
                    }
                }
            } else {
                /*
                pt is in the interior and nbd looks something like

                          l
                     +---------+
                    / \       / \
                l  /   \  s  /   \ l
                  /  s  \   /  s  \
                 /       \ /       \
                +---------p---------+
                 \       / \       /
                  \  s  /   \  s  /
                l  \   /  s  \   / l
                    \ /       \ /
                     +---------+
                          l

                where
                    p = point
                    s = star
                    l = link
                */

                if let Some(star_label) = self.common_label(star.iter()).cloned() {
                    let boundary = link_img
                        .iter()
                        .filter(|spx| {
                            let n = spx.n();
                            let a = nbd_affine_subspace.embedded_space().affine_dimension();
                            if n >= a {
                                unreachable!()
                            } else if n + 1 == a {
                                true
                            } else {
                                debug_assert!(n + 1 < a);
                                false
                            }
                        })
                        .map(|spx| {
                            OrientedSimplex::new_with_negative_point(
                                nbd_affine_subspace.embedded_space(),
                                spx.points().clone(),
                                &pt_img,
                            )
                            .unwrap()
                        })
                        .collect();

                    if let Some(new_star_img) =
                        simplify_in_region(nbd_affine_subspace.embedded_space(), boundary)
                    {
                        self.remove_simplexes_unchecked(star.into_iter().collect());
                        self.add_simplexes_unchecked(
                            new_star_img
                                .into_iter()
                                .map(|spx_img| nbd_affine_subspace.embed_simplex(&spx_img))
                                .collect(),
                            &star_label,
                        );
                        pts_todo.extend(link_points);
                    }
                }
            }
        }

        #[cfg(debug_assertions)]
        self.check();

        self
    }

    // #[deprecated]
    // pub fn frick(self) -> LabelledSimplicialComplex<FS, SP, usize> {
    //     let mut count = 0;
    //     LabelledSimplicialComplex::new_labelled(
    //         self.ambient_space().clone(),
    //         self.simplexes
    //             .iter()
    //             .map(|(spx, _)| {
    //                 count += 1;
    //                 (spx.clone(), count)
    //             })
    //             .collect(),
    //     )
    //     .unwrap()
    // }

    //remove simplexes and remove them from the inverse boundary of any others
    //self may not be in a valid state after this operation
    fn remove_simplexes_unchecked(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
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

    fn remove_simplexes(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
        self.remove_simplexes_unchecked(simplexes);
        #[cfg(debug_assertions)]
        self.check();
    }

    //add the given simplexes and add them to the inverse boundary map on anything on their boundaries
    //must be added together to cover the case where there are mutual boundary relations
    //self may not be in a valid state after this operation
    fn add_simplexes_unchecked(&mut self, simplexes: Vec<Simplex<FS, SP>>, label: &T) {
        for spx in &simplexes {
            self.simplexes.insert(
                spx.clone(),
                SCSpxInfo {
                    inv_bdry: HashSet::new(),
                    label: label.clone(),
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

    fn add_simplexes(&mut self, simplexes: Vec<Simplex<FS, SP>>, label: &T) {
        self.add_simplexes_unchecked(simplexes, label);
        #[cfg(debug_assertions)]
        self.check();
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ClosureLabel {
    Origional,
    Closure,
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

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    pub fn simplexes(&self) -> &HashSet<Simplex<FS, SP>> {
        &self.simplexes
    }

    pub fn into_simplexes(self) -> HashSet<Simplex<FS, SP>> {
        self.simplexes
    }

    pub fn try_as_simplicial_complex(self) -> Result<SimplicialComplex<FS, SP>, &'static str> {
        SimplicialComplex::new(self.ambient_space, self.simplexes)
    }

    pub fn closure(&self) -> LabelledSimplicialComplex<FS, SP, ClosureLabel> {
        let mut simplexes = HashSet::new();
        for spx in &self.simplexes {
            for bdry in spx.sub_simplices_not_null() {
                simplexes.insert(bdry);
            }
        }
        LabelledSimplicialComplex::new_labelled(
            self.ambient_space(),
            simplexes
                .into_iter()
                .map(|spx| {
                    let label = match self.simplexes.contains(&spx) {
                        true => ClosureLabel::Origional,
                        false => ClosureLabel::Closure,
                    };
                    (spx, label)
                })
                .collect(),
        )
        .unwrap()
    }

    pub fn simplify(&self) -> Self {
        self.closure()
            .simplify()
            .labelled_subset(&simplexes::ClosureLabel::Origional)
    }
}
