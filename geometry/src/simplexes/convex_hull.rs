use super::*;
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct ConvexHull<
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
> where
    FS::Set: Hash,
{
    // the space in which this convex hull lives
    ambient_space: SP,
    // the affine subspace spanned by this convex hull
    // so that this convex hull is "full" in this embedded subspace
    subspace: EmbeddedAffineSubspace<FS, SP, AffineSpace<FS>>,
    // oriented facets belonging to the embedded subspace such that
    // the positive side of each facet is on the interior of the convex hull and
    // the negative side of each facet is on the outside of the convex hull.
    // These facets should form a simplicial complex
    facets: Vec<OrientedSimplex<FS, AffineSpace<FS>>>,
    // These interior simplicies are the full-dimensional simplicies in the embedded subspace forming the interior of the convex hull
    interior: Vec<Simplex<FS, AffineSpace<FS>>>,
    /*
    Consider the case of a convex hull given by a simplex in dimension d:

    d    | #facets | #interior
    3    | 4 tri   | 1 tet
    2    | 3 line  | 1 tri
    1    | 2 pt    | 1 line
    0    | 1 null  | 1 pt
    null | 0 n/a   | 1 null

    This highlights what the behavour should been in the case where the dimension is null and 0
    */
}

impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone> std::fmt::Debug
    for ConvexHull<FS, SP>
where
    FS::Set: std::hash::Hash + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ConvexHull")
            .field(
                "embedded_in_dim",
                &self.subspace.embedded_space().affine_dimension(),
            )
            .field("facets", &self.facets)
            .field("interior", &self.interior)
            .finish()
    }
}

// #[derive(Debug, Clone)]
// pub struct ConvexHullAsSimplicialComplexResult<
//     FS: OrderedRingStructure + FieldStructure,
//     SP: Borrow<AffineSpace<FS>> + Clone,
// > {
//     pub entire: SimplicialComplex<FS, SP>,
//     pub boundary: SimplicialComplex<FS, SP>,
//     pub interior: PartialSimplicialComplex<FS, SP>,
// }

impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone>
    ConvexHull<FS, SP>
where
    FS::Set: Hash,
{
    fn check(&self) -> Result<(), &'static str> {
        assert_eq!(
            self.subspace.ambient_space().borrow(),
            self.ambient_space.borrow()
        );

        {
            for facet in &self.facets {
                if facet.ambient_space() != self.subspace.embedded_space() {
                    return Err("Facet must belong to the embedded subspace");
                }
            }
            //interior simplicies must have dimenion equal to self.subspace
            for spx in &self.interior {
                if spx.ambient_space() != self.subspace.embedded_space() {
                    return Err("Interior simplex must belong to the embedded subspace");
                }
                if spx.n() != self.subspace.embedded_space().affine_dimension() {
                    return Err("Interior simplex must span the embedded subspace");
                }
            }
        }

        match self.subspace.borrow().embedded_space().affine_dimension() {
            0 => {
                if self.facets.len() != 0 {
                    return Err("Empty convex hull should have no facets");
                }
                if self.interior
                    != vec![Simplex::new(self.subspace.embedded_space(), vec![]).unwrap()]
                {
                    return Err(
                        "Empty convex hull should have a single null simplex for its interior",
                    );
                }
            }
            1 => {
                if self.facets.len() != 1 {
                    return Err("0D convex hull should have one null facet");
                }
                if self.interior
                    != vec![
                        Simplex::new(
                            self.subspace.embedded_space(),
                            vec![Vector::construct(
                                self.subspace.embedded_space(),
                                |_| unreachable!(),
                            )],
                        )
                        .unwrap(),
                    ]
                {
                    return Err("0D convex hull should have one point for its interior");
                }
            }
            _ => {}
        }

        //facets should be non-empty whenenver self.subspace has dimension >= 1
        match self.subspace.borrow().embedded_space().linear_dimension() {
            Some(_) => {
                if self.facets.is_empty() {
                    return Err("Facets should be non-empty whenenver the subspace is non-empty");
                }
            }
            None => {}
        }

        //check that facets each share exactly one ridge
        {
            let mut ridges_count = HashMap::new();
            for facet in &self.facets {
                for ridge in facet.simplex().facets() {
                    if !ridges_count.contains_key(&ridge) {
                        ridges_count.insert(ridge.clone(), 0);
                    }
                    *ridges_count.get_mut(&ridge).unwrap() += 1;
                }
            }
            if !ridges_count.into_iter().all(|(_ridge, count)| count == 2) {
                return Err("Ridges of facets should each be shared between exactly two facets");
            }
        }

        //check that facets are convex with negative side inwards by checking that no point is on the positive side of a facet
        {
            let mut all_pts = HashSet::new();
            for facet in &self.facets {
                for pt in facet.simplex().points() {
                    all_pts.insert(pt);
                }
            }
            for spx in &self.interior {
                for pt in spx.points() {
                    all_pts.insert(pt);
                }
            }
            for facet in &self.facets {
                for pt in &all_pts {
                    match facet.classify_point(pt) {
                        OrientationSide::Negative => {
                            return Err(
                                "Every point must be on the positive or neutral side of every facet",
                            );
                        }
                        OrientationSide::Neutral | OrientationSide::Positive => {}
                    }
                }
            }
        }

        Ok(())
    }

    pub fn new_empty(ambient_space: SP) -> Self {
        let subspace = EmbeddedAffineSubspace::new_empty(ambient_space.clone());
        Self {
            ambient_space: ambient_space.clone(),
            subspace: subspace.clone(),
            facets: vec![],
            interior: vec![Simplex::new(subspace.embedded_space(), vec![]).unwrap()],
        }
    }

    pub fn new(ambient_space: SP, points: Vec<Vector<FS, SP>>) -> Self {
        let mut ch = Self::new_empty(ambient_space);
        for point in points {
            ch.extend_by_point(point);
        }
        ch
    }

    pub fn affine_span_dimension(&self) -> usize {
        self.subspace.embedded_space().affine_dimension()
    }

    pub fn from_simplex(spx: Simplex<FS, SP>) -> Self {
        let ambient_space = spx.ambient_space();
        let (subspace, embedded_pts) =
            EmbeddedAffineSubspace::new_affine_span(ambient_space.clone(), spx.into_points())
                .unwrap();
        let embedded_spx = Simplex::new(subspace.embedded_space(), embedded_pts).unwrap();
        Self {
            ambient_space,
            subspace,
            facets: embedded_spx.oriented_facets(),
            interior: vec![embedded_spx],
        }
    }

    pub fn ambient_space(&self) -> SP {
        self.subspace.ambient_space()
    }

    pub fn is_empty(&self) -> bool {
        self.subspace.borrow().embedded_space().affine_dimension() == 0
    }

    pub fn defining_points(&self) -> HashSet<Vector<FS, SP>> {
        match self.subspace.borrow().embedded_space().affine_dimension() {
            0 => HashSet::new(),
            1 => {
                //need to handle affine_dim = 1 case seperately from higher dimensional cases
                //because the facet in the 1d case is the null simplex with no points
                let (root, span) = self.subspace.get_root_and_span().unwrap();
                debug_assert_eq!(span.len(), 0);
                HashSet::from([root])
            }
            _ => {
                let mut points = HashSet::new();
                for facet in &self.facets {
                    for point in facet.simplex().points() {
                        points.insert(self.subspace.embed_point(point));
                    }
                }
                points
            }
        }
    }

    pub fn extend_by_point(&mut self, pt: Vector<FS, SP>) {
        assert_eq!(pt.ambient_space().borrow(), self.ambient_space.borrow());
        #[cfg(debug_assertions)]
        self.check().unwrap();

        // println!(
        //     "subspace_dim={:?}",
        //     self.subspace.embedded_space().linear_dimension()
        // );

        match self.subspace.unembed_point(&pt) {
            Some(subsp_pt) => {
                //Partition into visible / hidden facets
                //Compute horizon as ridges where visible meets hidden
                //Delete visible facets and extend horizion ridges to the new point for replacement facets

                let mut visible = vec![];
                let mut hidden = vec![];
                for facet in &self.facets {
                    match facet.classify_point(&subsp_pt) {
                        OrientationSide::Negative => {
                            visible.push(facet);
                        }
                        OrientationSide::Neutral | OrientationSide::Positive => {
                            hidden.push(facet);
                        }
                    }
                }

                let mut horizon = HashMap::new(); //horizon ridge -> (its adjacent facet, the other vertex)
                for facet in &visible {
                    let facet_simplex = facet.simplex();
                    for i in 0..facet_simplex.n() {
                        let (ridge, pt) = (facet_simplex.facet(i), &facet_simplex.points()[i]);
                        match horizon.contains_key(&ridge) {
                            true => {
                                horizon.remove(&ridge);
                            }
                            false => {
                                horizon.insert(ridge, (facet, pt));
                            }
                        }
                    }
                }
                #[cfg(debug_assertions)]
                {
                    let mut horizon_alt = HashSet::new();
                    for facet in &hidden {
                        for ridge in facet.simplex().facets() {
                            match horizon_alt.contains(&ridge) {
                                true => {
                                    horizon_alt.remove(&ridge);
                                }
                                false => {
                                    horizon_alt.insert(ridge);
                                }
                            }
                        }
                    }
                    assert_eq!(
                        horizon.keys().map(|r| r.clone()).collect::<HashSet<_>>(),
                        horizon_alt
                    );
                }

                // println!(
                //     "#visible={:?} #hidden={:?} #horizon={:?}",
                //     visible.len(),
                //     hidden.len(),
                //     horizon.len()
                // );

                (self.facets, self.interior) = (
                    horizon
                        .iter()
                        .map(|(ridge, (_facet, facet_pt))| {
                            OrientedSimplex::new_with_positive_point(
                                self.subspace.embedded_space(),
                                {
                                    let mut points = ridge.points().clone();
                                    points.push(subsp_pt.clone());
                                    points
                                },
                                *facet_pt,
                            )
                            .unwrap()
                        })
                        .chain(hidden.into_iter().map(|f| f.clone()))
                        .collect::<Vec<_>>(),
                    visible
                        .into_iter()
                        .map(|facet| {
                            Simplex::new(self.subspace.embedded_space(), {
                                let mut points = facet.simplex().points().clone();
                                points.push(subsp_pt.clone());
                                points
                            })
                            .unwrap()
                        })
                        .chain(self.interior.iter().map(|s| s.clone()))
                        .collect::<Vec<_>>(),
                );
            }
            None => {
                (self.subspace, self.facets, self.interior) = {
                    //The new point is outside the current embedded affine subspace
                    //The new facets are given by the old interior union the old facets extended into the new dimension
                    //The new interior is given by the old interior extended into the new dimension
                    {
                        let (iota, new_subspace_embedding, pt_in_new_subspace) =
                            self.subspace.extend_dimension_by_point_unsafe(pt.clone());

                        let new_subspace = iota.ambient_space();
                        debug_assert_eq!(new_subspace, new_subspace_embedding.embedded_space());

                        //new_facets <- old_interior & old_facets
                        //new_interior <- old_interior
                        let mut new_facets = vec![];
                        let mut new_interior = vec![];
                        for old_facet in &self.facets {
                            debug_assert_eq!(
                                old_facet.ambient_space(),
                                self.subspace.embedded_space()
                            );
                            new_facets.push(
                                OrientedSimplex::new_with_positive_point(
                                    new_subspace.clone(),
                                    {
                                        let mut points = vec![];
                                        for pt in old_facet.simplex().points() {
                                            points.push(iota.embed_point(pt));
                                        }
                                        points.push(pt_in_new_subspace.clone());
                                        points
                                    },
                                    &{
                                        //If old_facet is null living inside 0D space then take the iota-embedding of the unique point in the 0D space as the reference
                                        //If old_facet is not null then it comes equiped with a reference point which we embed via iota and use
                                        match old_facet.positive_point() {
                                            Some(pos_pt) => iota.embed_point(&pos_pt),
                                            None => {
                                                debug_assert_eq!(
                                                    self.subspace
                                                        .embedded_space()
                                                        .affine_dimension(),
                                                    1
                                                );
                                                iota.embed_point(
                                                    &self
                                                        .subspace
                                                        .embedded_space()
                                                        .origin()
                                                        .unwrap(),
                                                )
                                            }
                                        }
                                    },
                                )
                                .unwrap(),
                            );
                        }
                        for old_interior in &self.interior {
                            debug_assert_eq!(
                                old_interior.ambient_space(),
                                self.subspace.embedded_space()
                            );
                            new_facets.push(
                                OrientedSimplex::new_with_positive_point(
                                    new_subspace.clone(),
                                    old_interior
                                        .points()
                                        .iter()
                                        .map(|pt| iota.embed_point(pt))
                                        .collect(),
                                    &pt_in_new_subspace,
                                )
                                .unwrap(),
                            );
                            new_interior.push(
                                Simplex::new(new_subspace.clone(), {
                                    let mut points = old_interior
                                        .points()
                                        .iter()
                                        .map(|pt| iota.embed_point(pt))
                                        .collect::<Vec<_>>();
                                    points.push(pt_in_new_subspace.clone());
                                    points
                                })
                                .unwrap(),
                            );
                        }

                        (new_subspace_embedding, new_facets, new_interior)
                    }
                }
            }
        };

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    fn embedded_interior_simplexes(&self) -> Vec<Simplex<FS, SP>> {
        self.interior
            .iter()
            .map(|spx| {
                Simplex::new(
                    self.ambient_space.clone(),
                    spx.points()
                        .iter()
                        .map(|pt| self.subspace.embed_point(pt))
                        .collect(),
                )
                .unwrap()
            })
            .collect()
    }

    fn embedded_facet_simplexes(&self) -> Vec<Simplex<FS, SP>> {
        self.facets
            .iter()
            .map(|ospx| {
                Simplex::new(
                    self.ambient_space.clone(),
                    ospx.simplex()
                        .points()
                        .iter()
                        .map(|pt| self.subspace.embed_point(pt))
                        .collect(),
                )
                .unwrap()
            })
            .collect()
    }

    pub fn as_simplicial_complex(
        &self,
    ) -> LabelledSimplicialComplex<FS, SP, InteriorBoundaryLabel> {
        let boundary_simplexes = self
            .embedded_facet_simplexes()
            .into_iter()
            .map(|spx| spx.sub_simplices_not_null())
            .flatten()
            .collect::<HashSet<_>>();

        let all_simplexes = boundary_simplexes
            .clone()
            .into_iter()
            .chain(
                self.embedded_interior_simplexes()
                    .into_iter()
                    .map(|spx| spx.sub_simplices_not_null())
                    .flatten(),
            )
            .collect::<HashSet<_>>();

        LabelledSimplicialComplex::new_labelled(
            self.ambient_space().clone(),
            all_simplexes
                .into_iter()
                .map(|spx| {
                    let label = match boundary_simplexes.contains(&spx) {
                        false => InteriorBoundaryLabel::Interior,
                        true => InteriorBoundaryLabel::Boundary,
                    };
                    (spx, label)
                })
                .collect(),
        )
        .unwrap()

        // let mut interior_simplexes = HashSet::new();
        // for spx in &all_simplexes {
        //     if !boundary_simplexes.contains(spx) {
        //         interior_simplexes.insert(spx.clone());
        //     }
        // }

        // let entire =
        //     LabelledSimplicialComplex::new(self.ambient_space.clone(), all_simplexes).unwrap();
        // ConvexHullAsSimplicialComplexResult {
        //     entire: entire.clone(),
        //     boundary: LabelledSimplicialComplex::new(
        //         self.ambient_space.clone(),
        //         boundary_simplexes.into_iter().collect(),
        //     )
        //     .unwrap(),
        //     interior: PartialSimplicialComplex::new_unchecked(
        //         self.ambient_space.clone(),
        //         interior_simplexes,
        //     ),
        // }
    }
}

#[derive(Clone)]
struct ConvexHullWireframe<
    FS: OrderedRingSignature + FieldSignature,
    SP: Borrow<AffineSpace<FS>> + Clone,
> {
    ambient_space: SP,
    points: Vec<Vector<FS, SP>>,
    edges: Vec<(Vector<FS, SP>, Vector<FS, SP>)>,
}

impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone>
    ConvexHullWireframe<FS, SP>
where
    FS::Set: Hash,
{
    fn from_convex_hull(ch: &ConvexHull<FS, SP>) -> Self {
        let mut outer_points = HashSet::new();
        let mut outer_edges = HashSet::new();
        for facet in &ch.facets {
            for point in facet.simplex().points() {
                outer_points.insert(point.clone());
            }
            for edge in facet.simplex().edges() {
                outer_edges.insert(edge);
            }
        }

        // type cast
        let mut outer_points = outer_points.into_iter().collect::<Vec<_>>();
        let mut outer_edges = outer_edges
            .into_iter()
            .map(|line| {
                let mut pts = line.into_points().into_iter();
                (pts.next().unwrap(), pts.next().unwrap())
            })
            .collect::<Vec<_>>();

        // Note that in linear dimension 0 we need to explicitly add the unique point to outer points.
        if ch.subspace.embedded_space().affine_dimension() == 1 {
            let (_root, span) = ch.subspace.get_root_and_span().unwrap();
            debug_assert_eq!(span.len(), 0);
            debug_assert_eq!(outer_points.len(), 0);
            outer_points.push(Vector::new(ch.subspace.embedded_space(), vec![]));
        }

        // Note that in linear dimension 1 this doesnt quite work for outer edges since the boundary is not connected. Instead, outer edges should just be the edge between the two points.
        if ch.subspace.embedded_space().affine_dimension() == 2 {
            debug_assert_eq!(outer_points.len(), 2);
            debug_assert_eq!(outer_edges.len(), 0);
            let mut pts = outer_points.iter();
            outer_edges.push((
                (*pts.next().unwrap()).clone(),
                (*pts.next().unwrap()).clone(),
            ));
        }

        Self {
            ambient_space: ch.ambient_space(),
            points: outer_points
                .into_iter()
                .map(|pt| ch.subspace.embed_point(&pt))
                .collect(),
            edges: outer_edges
                .into_iter()
                .map(|(pt1, pt2)| (ch.subspace.embed_point(&pt1), ch.subspace.embed_point(&pt2)))
                .collect(),
        }
    }

    fn into_convex_hull(self) -> ConvexHull<FS, SP> {
        ConvexHull::new(self.ambient_space, self.points)
    }

    pub fn intersect_with_oriented_hyperplane(
        &self,
        hyperplane: &OrientedHyperplane<FS, SP>,
        region: OrientationSide, //TODO: make this const generic once rust has const generic enums
    ) -> Self {
        /*
        Find positive_points, middle_points, negative_points such that
        ConvexHull(middle_points) = self intersect hyperplane
        ConvexHull(positive_points union middle_points) = self intersect (hyperplane union hyperplane_posside)
        ConvexHull(negative_points union middle_points) = self intersect (hyperplane union hyperplane_negside)
        */
        let mut positive_points = HashSet::new();
        let mut middle_points = HashSet::new();
        let mut negative_points = HashSet::new();
        for point in &self.points {
            match hyperplane.classify_point(&point) {
                OrientationSide::Positive => {
                    positive_points.insert(point.clone());
                }
                OrientationSide::Neutral => {
                    middle_points.insert(point.clone());
                }
                OrientationSide::Negative => {
                    negative_points.insert(point.clone());
                }
            }
        }
        let mut positive_edges = HashSet::new();
        let mut negative_edges = HashSet::new();
        for (a, b) in &self.edges {
            match hyperplane.intersect_line(&a, &b) {
                OrientedHyperplaneIntersectLineSegmentResult::PositivePositive
                | OrientedHyperplaneIntersectLineSegmentResult::PositiveNeutral
                | OrientedHyperplaneIntersectLineSegmentResult::NeutralPositive => {
                    positive_edges.insert((a.clone(), b.clone()));
                }
                OrientedHyperplaneIntersectLineSegmentResult::NegativeNegative
                | OrientedHyperplaneIntersectLineSegmentResult::NegativeNeutral
                | OrientedHyperplaneIntersectLineSegmentResult::NeutralNegative => {
                    negative_edges.insert((a.clone(), b.clone()));
                }

                OrientedHyperplaneIntersectLineSegmentResult::NeutralNeutral => {
                    //We deal with middle edges later
                }
                OrientedHyperplaneIntersectLineSegmentResult::PositiveNegative {
                    intersection_point: m,
                } => {
                    positive_edges.insert((a.clone(), m.clone()));
                    negative_edges.insert((m.clone(), b.clone()));
                    middle_points.insert(m);
                }
                OrientedHyperplaneIntersectLineSegmentResult::NegativePositive {
                    intersection_point: m,
                } => {
                    positive_edges.insert((b.clone(), m.clone()));
                    negative_edges.insert((m.clone(), a.clone()));
                    middle_points.insert(m);
                }
            }
        }

        let ConvexHullWireframe {
            ambient_space: _,
            points: middle_points,
            edges: middle_edges,
        } = Self::from_convex_hull(&ConvexHull::new(
            self.ambient_space.clone(),
            middle_points.into_iter().collect(),
        ));

        if !positive_points.is_empty() && !negative_points.is_empty() {
            debug_assert!(!middle_points.is_empty());
        }

        let mut points: HashSet<Vector<FS, SP>> = middle_points.into_iter().collect();
        let mut edges: HashSet<(Vector<FS, SP>, Vector<FS, SP>)> =
            middle_edges.into_iter().collect();

        match region {
            OrientationSide::Positive => {
                for pt in positive_points {
                    points.insert(pt);
                }
                for edge in positive_edges {
                    edges.insert(edge);
                }
            }
            OrientationSide::Neutral => {}
            OrientationSide::Negative => {
                for pt in negative_points {
                    points.insert(pt);
                }
                for edge in negative_edges {
                    edges.insert(edge);
                }
            }
        }

        Self {
            ambient_space: self.ambient_space.clone(),
            points: points.into_iter().collect(),
            edges: edges.into_iter().collect(),
        }
    }
}

impl<FS: OrderedRingSignature + FieldSignature, SP: Borrow<AffineSpace<FS>> + Clone>
    ConvexHull<FS, SP>
where
    FS::Set: Hash,
{
    pub fn intersect_with_oriented_hyperplane(
        &self,
        hyperplane: &OrientedHyperplane<FS, SP>,
        region: OrientationSide, //TODO: make this const generic once rust has const generic enums
    ) -> Self {
        ConvexHullWireframe::from_convex_hull(self)
            .intersect_with_oriented_hyperplane(hyperplane, region)
            .into_convex_hull()
    }

    pub fn intersect_mut(&mut self, other: &Self) {
        let ambient_space = self.ambient_space();
        assert_eq!(ambient_space.borrow(), other.ambient_space().borrow());
        match other.subspace.as_hyperplane_intersection() {
            Some(hyperplanes) => {
                // step 1: intersect with the affine space spanned by other
                for hyperplane in hyperplanes {
                    *self = self
                        .intersect_with_oriented_hyperplane(&hyperplane, OrientationSide::Neutral);
                }
                // step2: embed self into other.affine_subspace
                let embedded_self_in_other_subspace = &ConvexHull::new(
                    other.subspace.embedded_space(),
                    self.defining_points()
                        .into_iter()
                        .map(|point| other.subspace.unembed_point(&point).unwrap())
                        .collect(),
                );
                debug_assert_eq!(
                    embedded_self_in_other_subspace
                        .subspace
                        .embedded_space()
                        .affine_dimension(),
                    self.subspace.embedded_space().affine_dimension()
                );
                let mut embedded_self_in_other_subspace =
                    ConvexHullWireframe::from_convex_hull(&embedded_self_in_other_subspace);
                // step3 : intersect self with each oriented facet of other
                for facet in &other.facets {
                    embedded_self_in_other_subspace = embedded_self_in_other_subspace
                        .intersect_with_oriented_hyperplane(
                            &facet.clone().into_oriented_hyperplane(),
                            OrientationSide::Positive,
                        );
                }

                // step4: unembed self back to the ambient space
                *self = ConvexHull::new(
                    self.ambient_space(),
                    embedded_self_in_other_subspace
                        .into_convex_hull()
                        .defining_points()
                        .into_iter()
                        .map(|point| other.subspace.embed_point(&point))
                        .collect(),
                );
            }
            None => {
                //other is empty, so the intersection with other is empty
                *self = Self::new_empty(self.ambient_space.clone());
            }
        }
    }

    pub fn intersect(&self, other: &Self) -> Self {
        let mut self_mut = self.clone();
        self_mut.intersect_mut(other);
        self_mut
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;
    use algebraeon_sets::structure::*;

    #[test]
    fn construct_convex_hull() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let mut ch = ConvexHull::new_empty(&space);
        //add a point to get started
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(1)],
        ));
        //add the same point again
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(1)],
        ));
        //add a point to form a line
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
        ));
        //add a point such that the line is extended to a longer line
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(-1), Rational::from(0)],
        ));
        //add that point again
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(-1), Rational::from(0)],
        ));
        //add the middle point of the line again
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
        ));
        //add a middle point in the middle of one of the two line segments
        ch.extend_by_point(Vector::new(
            &space,
            vec![
                Rational::from(1) / Rational::from(2),
                Rational::from(3) / Rational::from(4),
            ],
        ));
        //add a point to form a double triangle
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(2), Rational::from(-1)],
        ));
        //add a point to extend by one triangle
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(0), Rational::from(-2)],
        ));
        //extend by a triangle such that the boundary ends up with two points on the same straight edge
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(2), Rational::from(0)],
        ));
        //add a point to extend by three triangles
        ch.extend_by_point(Vector::new(
            &space,
            vec![Rational::from(2), Rational::from(2)],
        ));
    }

    #[test]
    fn convex_hull_intersect_hyperplane() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let ch = ConvexHull::new(
            &space,
            vec![
                Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                Vector::new(
                    &space,
                    vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                ),
                Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                Vector::new(
                    &space,
                    vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                ),
                Vector::new(
                    &space,
                    vec![
                        Rational::from(1) / Rational::from(2),
                        Rational::from(3) / Rational::from(4),
                    ],
                ),
                Vector::new(&space, vec![Rational::from(2), Rational::from(-1)]),
                Vector::new(&space, vec![Rational::from(0), Rational::from(-2)]),
                Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
                Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
            ],
        );

        let ohsp = OrientedSimplex::new_with_positive_point(
            &space,
            vec![
                Vector::new(&space, vec![Rational::from(1), Rational::from(4)]),
                Vector::new(&space, vec![Rational::from(1), Rational::from(-4)]),
            ],
            &Vector::new(&space, vec![Rational::from(10), Rational::from(0)]),
        )
        .unwrap()
        .into_oriented_hyperplane();

        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Neutral);
        println!("{:?}", smaller_ch);
        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Positive);
        println!("{:?}", smaller_ch);
        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Negative);
        println!("{:?}", smaller_ch);
    }

    #[test]
    fn convex_hull_intersections() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        //2d intersect 2d
        {
            let ch1 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                    Vector::new(
                        &space,
                        vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                    ),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                    Vector::new(
                        &space,
                        vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                    ),
                    Vector::new(
                        &space,
                        vec![
                            Rational::from(1) / Rational::from(2),
                            Rational::from(3) / Rational::from(4),
                        ],
                    ),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(-1)]),
                    Vector::new(&space, vec![Rational::from(0), Rational::from(-2)]),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
                ],
            );
            let ch2 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(-2), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(3), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(3), Rational::from(2)]),
                    Vector::new(&space, vec![Rational::from(0), Rational::from(-1)]),
                ],
            );
            let ch3 = ch1.intersect(&ch2);
            println!("{:?}", ch3);
        }
        //2d intersect 1d
        {
            let ch1 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                    Vector::new(
                        &space,
                        vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                    ),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                    Vector::new(
                        &space,
                        vec![Rational::from(0), Rational::from(1) / Rational::from(2)],
                    ),
                    Vector::new(
                        &space,
                        vec![
                            Rational::from(1) / Rational::from(2),
                            Rational::from(3) / Rational::from(4),
                        ],
                    ),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(-1)]),
                    Vector::new(&space, vec![Rational::from(0), Rational::from(-2)]),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(2), Rational::from(2)]),
                ],
            );
            let ch2 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(3), Rational::from(2)]),
                    Vector::new(&space, vec![Rational::from(0), Rational::from(-1)]),
                ],
            );
            let ch3 = ch1.intersect(&ch2);
            println!("{:?}", ch3);
        }
        //line misses line
        {
            let ch1 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(-2), Rational::from(-1)]),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(1)]),
                ],
            );
            let ch2 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                    Vector::new(&space, vec![Rational::from(-1), Rational::from(0)]),
                ],
            );
            let ch3 = ch1.intersect(&ch2);
            debug_assert_eq!(ch3.subspace.embedded_space().affine_dimension(), 0);
        }
        //line hits line
        {
            let ch1 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(2), Rational::from(0)]),
                    Vector::new(&space, vec![Rational::from(-2), Rational::from(-1)]),
                ],
            );
            let ch2 = ConvexHull::new(
                &space,
                vec![
                    Vector::new(&space, vec![Rational::from(0), Rational::from(-1)]),
                    Vector::new(&space, vec![Rational::from(1), Rational::from(1)]),
                ],
            );
            let ch3 = ch1.intersect(&ch2);
            debug_assert_eq!(ch3.subspace.embedded_space().affine_dimension(), 1);
        }
    }

    #[test]
    fn convex_hull_from_simplex() {
        let space = AffineSpace::new_linear(Rational::structure(), 3);
        let p1 = Vector::new(
            &space,
            vec![Rational::from(4), Rational::from(2), Rational::from(2)],
        );
        let p2 = Vector::new(
            &space,
            vec![Rational::from(5), Rational::from(-3), Rational::from(3)],
        );
        let p3 = Vector::new(
            &space,
            vec![Rational::from(-5), Rational::from(6), Rational::from(2)],
        );
        let p4 = Vector::new(
            &space,
            vec![Rational::from(8), Rational::from(2), Rational::from(-9)],
        );

        ConvexHull::from_simplex(Simplex::new(&space, vec![]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(Simplex::new(&space, vec![p1.clone()]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(Simplex::new(&space, vec![p1.clone(), p2.clone()]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(
            Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone()]).unwrap(),
        )
        .check()
        .unwrap();
        ConvexHull::from_simplex(
            Simplex::new(&space, vec![p1.clone(), p2.clone(), p3.clone(), p4.clone()]).unwrap(),
        )
        .check()
        .unwrap();
    }
}
