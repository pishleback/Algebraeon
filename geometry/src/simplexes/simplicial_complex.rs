use std::collections::{HashMap, HashSet};

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

    pub fn interior_and_boundary(
        &self,
    ) -> (PartialSimplicialComplex<FS, SP>, SimplicialComplex<FS, SP>) {
        /*
        let n be the dimension of the space self is living in
         - every simplex of rank n is part of the interior
         - a simplex of rank n-1 is the facet of at most 2 simplices of rank n, and is part of the interior if and only if it is the facet of exactly 2 simplices of rank n
         - a simplex of rank less or equal to n-2 is part of the interior iff it is in the boundary of some strictly higher rank simplex AND every strictly higher rank simplex containing it as part of the boundary is part of the interior
        */

        let n = self.ambient_space().borrow().affine_dimension();

        let mut interior: HashSet<&Simplex<FS, SP>> = HashSet::new();
        let mut boundary: HashSet<&Simplex<FS, SP>> = HashSet::new();

        let mut all = self.simplexes.keys().collect::<Vec<_>>();
        all.sort_unstable_by_key(|s| std::cmp::Reverse(s.n())); //so that we process largest rank first
        for simplex in all {
            let r = simplex.n();
            if r == n {
                //rank n simplex is always part of the interior
                interior.insert(simplex);
            } else {
                let inv_bdry = &self.simplexes.get(simplex).unwrap().inv_bdry;
                if r == n - 1 {
                    //rank n-1 simplex is part of the boundary of at most 2 rank n simplices
                    //it is part of the boundary of the simplicial complex iff it is the boundary of exactly 2
                    match inv_bdry.len() {
                        0 | 1 => {
                            boundary.insert(simplex);
                        }
                        2 => {
                            interior.insert(simplex);
                        }
                        _ => panic!(
                        "rank n-1 simplex should be in the boundary of at most 2 rank n simplices"
                    ),
                    }
                } else {
                    //rank < n-1 simplex is part of the interior iff it is part of the boundary of at least one simplex and every such simplex is part of the interior
                    debug_assert!(r < n - 1);
                    if inv_bdry.is_empty() {
                        boundary.insert(simplex);
                    } else {
                        if inv_bdry.iter().all(|b| interior.contains(b)) {
                            interior.insert(simplex);
                        } else {
                            boundary.insert(simplex);
                        }
                    }
                }
            }
        }

        (
            PartialSimplicialComplex::new_unchecked(
                self.ambient_space(),
                interior.into_iter().map(|s| s.clone()).collect(),
            ),
            SimplicialComplex::new(
                self.ambient_space(),
                boundary.into_iter().map(|s| s.clone()).collect(),
            )
            .unwrap(),
        )
    }

    pub fn interior(&self) -> PartialSimplicialComplex<FS, SP> {
        self.interior_and_boundary().0
    }

    pub fn boundary(&self) -> SimplicialComplex<FS, SP> {
        self.interior_and_boundary().1
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

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialComplex<FS, SP>
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

            let nbd = SimplicialComplex::new(nbd_affine_subspace.embedded_space(), {
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

                    let pt_img_img = nbd_boundary_affine_subspace.unembed_point(&pt_img).unwrap();

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
                            );
                            self.add_simplexes_unchecked(
                                new_star_img
                                    .into_iter()
                                    .map(|spx_img| nbd_affine_subspace.embed_simplex(&spx_img))
                                    .collect(),
                            );
                            pts_todo.extend(link_points);
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
                    );
                    pts_todo.extend(link_points);
                }
            }
        }

        #[cfg(debug_assertions)]
        self.check();

        self
    }

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
    fn add_simplexes_unchecked(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
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

    fn add_simplexes(&mut self, simplexes: Vec<Simplex<FS, SP>>) {
        self.add_simplexes_unchecked(simplexes);
        #[cfg(debug_assertions)]
        self.check();
    }

    pub fn union(&self, other: &Self) -> Self {
        SimplicialDisjointUnion::union(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .closure_as_simplicial_complex()
    }

    pub fn intersection(&self, other: &Self) -> Self {
        SimplicialDisjointUnion::intersection(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .closure_as_simplicial_complex()
    }
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

    pub fn closure_as_simplicial_complex(&self) -> SimplicialComplex<FS, SP> {
        let mut simplexes = HashSet::new();
        for spx in &self.simplexes {
            for bdry in spx.sub_simplices_not_null() {
                simplexes.insert(bdry);
            }
        }
        SimplicialComplex::new(self.ambient_space(), simplexes).unwrap()
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
        // let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        VennResult {
            left: Self::subtract(self, other),
            middle: Self::intersection(self, other),
            right: Self::subtract(other, self),
        }
    }

    pub fn intersection(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        Self::new_unchecked(ambient_space.clone(), {
            let mut simplexes = HashSet::new();
            for self_spx in &self.simplexes {
                for other_spx in &other.simplexes {
                    for spx in Simplex::venn(self_spx, other_spx).middle.into_simplexes() {
                        simplexes.insert(spx);
                    }
                }
            }
            simplexes
        })
    }

    pub fn union(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let mut simplexes = HashSet::new();
        for spx in Self::subtract(other, self).into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in self.simplexes() {
            simplexes.insert(spx.clone());
        }
        return Self::new_unchecked(ambient_space, simplexes);
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

    pub fn closure_as_simplicial_complex(self) -> SimplicialComplex<FS, SP> {
        self.refine_to_partial_simplicial_complex()
            .closure_as_simplicial_complex()
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

        if middle.len() == 0 {
            VennResult {
                left: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::from([self.clone()]),
                ),
                middle: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::new(),
                ),
                right: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::from([other.clone()]),
                ),
            }
        } else {
            VennResult {
                left: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), left),
                middle: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), middle),
                right: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), right),
            }
        }
    }
}
