use crate::{
    affine_subspace::EmbeddedAffineSubspace,
    ambient_space::AffineSpace,
    oriented_simplex::{
        OrientationSide, OrientedHyperplane, OrientedHyperplaneIntersectLineSegmentResult,
        OrientedSimplex,
    },
    simplex::Simplex,
    simplex_collection::{InteriorOrBoundary, LabelledSimplexCollection},
    simplicial_complex::LabelledSimplicialComplex,
    vector::Vector,
};

use super::*;
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct ConvexHull<'f, FS: OrderedRingSignature + FieldSignature>
where
    FS::Set: Hash,
{
    // the space in which this convex hull lives
    ambient_space: AffineSpace<'f, FS>,
    // the affine subspace spanned by this convex hull
    // so that this convex hull is "full" in this embedded subspace
    subspace: EmbeddedAffineSubspace<'f, FS>,
    // oriented facets belonging to the embedded subspace such that
    // the positive side of each facet is on the interior of the convex hull and
    // the negative side of each facet is on the outside of the convex hull.
    // These facets should form a simplicial complex
    facets: Vec<OrientedSimplex<'f, FS>>,
    // These interior simplicies are the full-dimensional simplicies in the embedded subspace forming the interior of the convex hull
    interior: Vec<Simplex<'f, FS>>,
    /*
    Consider the case of a convex hull given by a simplex in dimension d:

    d    | #facets | #interior
    3    | 4 tri   | 1 tet
    2    | 3 line  | 1 tri
    1    | 2 pt    | 1 line
    0    | 1 null  | 1 pt
    null | 0 n/a   | 1 null

    This highlights what the behaviour should been in the case where the dimension is null and 0
    */
}

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for ConvexHull<'f, FS>
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

impl<'f, FS: OrderedRingSignature + FieldSignature> ConvexHull<'f, FS>
where
    FS::Set: Hash,
{
    #[allow(unused)]
    fn check(&self) -> Result<(), &'static str> {
        assert_eq!(self.subspace.ambient_space(), self.ambient_space);

        {
            for facet in &self.facets {
                if facet.ambient_space() != self.subspace.embedded_space() {
                    return Err("Facet must belong to the embedded subspace");
                }
            }
            //interior simplicies must have dimension equal to self.subspace
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
                if !self.facets.is_empty() {
                    return Err("Empty convex hull should have no facets");
                }
                if self.interior != vec![self.subspace.embedded_space().simplex(vec![]).unwrap()] {
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
                        self.subspace
                            .embedded_space()
                            .simplex(vec![Vector::construct(
                                self.subspace.embedded_space(),
                                |_| unreachable!(),
                            )])
                            .unwrap(),
                    ]
                {
                    return Err("0D convex hull should have one point for its interior");
                }
            }
            _ => {}
        }

        //facets should be non-empty whenenver self.subspace has dimension >= 1
        if self
            .subspace
            .borrow()
            .embedded_space()
            .linear_dimension()
            .is_some()
            && self.facets.is_empty()
        {
            return Err("Facets should be non-empty whenenver the subspace is non-empty");
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

    pub(crate) fn new_empty(ambient_space: AffineSpace<'f, FS>) -> Self {
        let subspace = EmbeddedAffineSubspace::new_empty(ambient_space);
        Self {
            ambient_space,
            subspace: subspace.clone(),
            facets: vec![],
            interior: vec![subspace.clone().embedded_space().simplex(vec![]).unwrap()],
        }
    }

    pub(crate) fn new(ambient_space: AffineSpace<'f, FS>, points: Vec<Vector<'f, FS>>) -> Self {
        let mut ch = Self::new_empty(ambient_space);
        for point in points {
            ch.extend_by_point(point);
        }
        ch
    }

    pub(crate) fn from_simplex(spx: Simplex<'f, FS>) -> Self {
        let ambient_space = spx.ambient_space();
        let (subspace, embedded_pts) =
            EmbeddedAffineSubspace::new_affine_span(ambient_space, spx.into_points()).unwrap();
        let embedded_spx = subspace.embedded_space().simplex(embedded_pts).unwrap();
        Self {
            ambient_space,
            subspace,
            facets: embedded_spx.oriented_facets(),
            interior: vec![embedded_spx],
        }
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> AffineSpace<'f, FS>
where
    FS::Set: Hash,
{
    pub fn convex_hull(&self, points: Vec<Vector<'f, FS>>) -> ConvexHull<'f, FS> {
        ConvexHull::new(*self, points)
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> ConvexHull<'f, FS>
where
    FS::Set: Hash,
{
    pub fn affine_span_dimension(&self) -> usize {
        self.subspace.embedded_space().affine_dimension()
    }

    pub fn ambient_space(&self) -> AffineSpace<'f, FS> {
        self.subspace.ambient_space()
    }

    pub fn is_empty(&self) -> bool {
        self.subspace.borrow().embedded_space().affine_dimension() == 0
    }

    pub fn defining_points(&self) -> HashSet<Vector<'f, FS>> {
        match self.subspace.borrow().embedded_space().affine_dimension() {
            0 => HashSet::new(),
            1 => {
                //need to handle affine_dim = 1 case separately from higher dimensional cases
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

    pub fn extend_by_point(&mut self, pt: Vector<'f, FS>) {
        assert_eq!(pt.ambient_space(), self.ambient_space);
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
                        match horizon.entry(ridge) {
                            std::collections::hash_map::Entry::Occupied(entry) => {
                                entry.remove();
                            }
                            std::collections::hash_map::Entry::Vacant(entry) => {
                                entry.insert((facet, pt));
                            }
                        }
                    }
                }
                #[cfg(debug_assertions)]
                {
                    let mut horizon_alt = HashSet::new();
                    for facet in &hidden {
                        for ridge in facet.simplex().facets() {
                            if horizon_alt.contains(&ridge) {
                                horizon_alt.remove(&ridge);
                            } else {
                                horizon_alt.insert(ridge);
                            }
                        }
                    }
                    assert_eq!(horizon.keys().cloned().collect::<HashSet<_>>(), horizon_alt);
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
                        .chain(hidden.into_iter().cloned())
                        .collect::<Vec<_>>(),
                    visible
                        .into_iter()
                        .map(|facet| {
                            self.subspace
                                .embedded_space()
                                .simplex({
                                    let mut points = facet.simplex().points().clone();
                                    points.push(subsp_pt.clone());
                                    points
                                })
                                .unwrap()
                        })
                        .chain(self.interior.iter().cloned())
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
                                    new_subspace,
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
                                        //If old_facet is not null then it comes equipped with a reference point which we embed via iota and use
                                        if let Some(pos_pt) = old_facet.positive_point() {
                                            iota.embed_point(&pos_pt)
                                        } else {
                                            debug_assert_eq!(
                                                self.subspace.embedded_space().affine_dimension(),
                                                1
                                            );
                                            iota.embed_point(
                                                &self.subspace.embedded_space().origin().unwrap(),
                                            )
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
                                    new_subspace,
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
                                new_subspace
                                    .simplex({
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
        }

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    fn embedded_interior_simplexes(&self) -> Vec<Simplex<'f, FS>> {
        self.interior
            .iter()
            .map(|spx| {
                self.ambient_space
                    .simplex(
                        spx.points()
                            .iter()
                            .map(|pt| self.subspace.embed_point(pt))
                            .collect(),
                    )
                    .unwrap()
            })
            .collect()
    }

    fn embedded_facet_simplexes(&self) -> Vec<Simplex<'f, FS>> {
        self.facets
            .iter()
            .map(|ospx| {
                self.ambient_space
                    .simplex(
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

    pub fn to_simplicial_complex(&self) -> LabelledSimplicialComplex<'f, FS, InteriorOrBoundary> {
        let boundary_simplexes = self
            .embedded_facet_simplexes()
            .into_iter()
            .flat_map(|spx| spx.sub_simplices_not_null())
            .collect::<HashSet<_>>();

        let all_simplexes = boundary_simplexes
            .clone()
            .into_iter()
            .chain(
                self.embedded_interior_simplexes()
                    .into_iter()
                    .flat_map(|spx| spx.sub_simplices_not_null()),
            )
            .collect::<HashSet<_>>();

        LabelledSimplicialComplex::try_new_labelled(
            self.ambient_space(),
            all_simplexes
                .into_iter()
                .map(|spx| {
                    let label = if boundary_simplexes.contains(&spx) {
                        InteriorOrBoundary::Boundary
                    } else {
                        InteriorOrBoundary::Interior
                    };
                    (spx, label)
                })
                .collect(),
        )
        .unwrap()
    }
}

#[derive(Clone)]
struct ConvexHullWireframe<'f, FS: OrderedRingSignature + FieldSignature> {
    ambient_space: AffineSpace<'f, FS>,
    points: Vec<Vector<'f, FS>>,
    edges: Vec<(Vector<'f, FS>, Vector<'f, FS>)>,
}

impl<'f, FS: OrderedRingSignature + FieldSignature> ConvexHullWireframe<'f, FS>
where
    FS::Set: Hash,
{
    fn from_convex_hull(ch: &ConvexHull<'f, FS>) -> Self {
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
            outer_points.push(ch.subspace.embedded_space().vector(Vec::<FS::Set>::new()));
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

    fn into_convex_hull(self) -> ConvexHull<'f, FS> {
        ConvexHull::new(self.ambient_space, self.points)
    }

    pub fn intersect_with_oriented_hyperplane(
        &self,
        hyperplane: &OrientedHyperplane<'f, FS>,
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
            match hyperplane.classify_point(point) {
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
            match hyperplane.intersect_line(a, b) {
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
            self.ambient_space,
            middle_points.into_iter().collect(),
        ));

        if !positive_points.is_empty() && !negative_points.is_empty() {
            debug_assert!(!middle_points.is_empty());
        }

        let mut points: HashSet<Vector<'f, FS>> = middle_points.into_iter().collect();
        let mut edges: HashSet<(Vector<'f, FS>, Vector<'f, FS>)> =
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
            ambient_space: self.ambient_space,
            points: points.into_iter().collect(),
            edges: edges.into_iter().collect(),
        }
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> ConvexHull<'f, FS>
where
    FS::Set: Hash,
{
    pub fn intersect_with_oriented_hyperplane(
        &self,
        hyperplane: &OrientedHyperplane<'f, FS>,
        region: OrientationSide, //TODO: make this const generic once rust has const generic enums
    ) -> Self {
        ConvexHullWireframe::from_convex_hull(self)
            .intersect_with_oriented_hyperplane(hyperplane, region)
            .into_convex_hull()
    }

    pub fn intersect_mut(&mut self, other: &Self) {
        let ambient_space = self.ambient_space();
        assert_eq!(ambient_space, other.ambient_space());
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
                    ConvexHullWireframe::from_convex_hull(embedded_self_in_other_subspace);

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
                *self = Self::new_empty(self.ambient_space);
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

    #[test]
    fn construct_convex_hull() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let mut ch = ConvexHull::new_empty(space);
        //add a point to get started
        ch.extend_by_point(space.vector([1, 1]));
        //add the same point again
        ch.extend_by_point(space.vector([1, 1]));
        //add a point to form a line
        ch.extend_by_point(
            space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
        );
        //add a point such that the line is extended to a longer line
        ch.extend_by_point(space.vector([-1, 0]));
        //add that point again
        ch.extend_by_point(space.vector([-1, 0]));
        //add the middle point of the line again
        ch.extend_by_point(
            space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
        );
        //add a middle point in the middle of one of the two line segments
        ch.extend_by_point(space.vector([
            Rational::from(1) / Rational::from(2),
            Rational::from(3) / Rational::from(4),
        ]));
        //add a point to form a double triangle
        ch.extend_by_point(space.vector([2, -1]));
        //add a point to extend by one triangle
        ch.extend_by_point(space.vector([0, -2]));
        //extend by a triangle such that the boundary ends up with two points on the same straight edge
        ch.extend_by_point(space.vector([2, 0]));
        //add a point to extend by three triangles
        ch.extend_by_point(space.vector([2, 2]));
    }

    #[test]
    fn convex_hull_intersect_hyperplane() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let ch = ConvexHull::new(
            space,
            vec![
                space.vector([1, 1]),
                space.vector([1, 1]),
                space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                space.vector([-1, 0]),
                space.vector([-1, 0]),
                space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                space.vector([
                    Rational::from(1) / Rational::from(2),
                    Rational::from(3) / Rational::from(4),
                ]),
                space.vector([2, -1]),
                space.vector([0, -2]),
                space.vector([2, 0]),
                space.vector([2, 2]),
            ],
        );

        let ohsp = OrientedSimplex::new_with_positive_point(
            space,
            vec![
                space.vector([Rational::from(1), Rational::from(4)]),
                space.vector([Rational::from(1), Rational::from(-4)]),
            ],
            &space.vector(vec![Rational::from(10), Rational::from(0)]),
        )
        .unwrap()
        .into_oriented_hyperplane();

        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Neutral);
        println!("{smaller_ch:?}");
        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Positive);
        println!("{smaller_ch:?}");
        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Negative);
        println!("{smaller_ch:?}");
    }

    #[test]
    fn convex_hull_intersections() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        //2d intersect 2d
        {
            let ch1 = ConvexHull::new(
                space,
                vec![
                    space.vector([1, 1]),
                    space.vector([1, 1]),
                    space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                    space.vector([-1, 0]),
                    space.vector([-1, 0]),
                    space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                    space.vector([
                        Rational::from(1) / Rational::from(2),
                        Rational::from(3) / Rational::from(4),
                    ]),
                    space.vector([2, 1]),
                    space.vector([0, -2]),
                    space.vector([2, 0]),
                    space.vector([2, 2]),
                ],
            );
            let ch2 = ConvexHull::new(
                space,
                vec![
                    space.vector([-2, 0]),
                    space.vector([3, 0]),
                    space.vector([3, 2]),
                    space.vector([0, -1]),
                ],
            );
            let ch3 = ch1.intersect(&ch2);
            println!("{ch3:?}");
        }
        //2d intersect 1d
        {
            let ch1 = ConvexHull::new(
                space,
                vec![
                    space.vector([1, 1]),
                    space.vector([1, 1]),
                    space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                    space.vector([-1, 0]),
                    space.vector([-1, 0]),
                    space.vector([Rational::from(0), Rational::from(1) / Rational::from(2)]),
                    space.vector([
                        Rational::from(1) / Rational::from(2),
                        Rational::from(3) / Rational::from(4),
                    ]),
                    space.vector([2, 1]),
                    space.vector([0, -2]),
                    space.vector([2, 0]),
                    space.vector([2, 2]),
                ],
            );
            let ch2 = ConvexHull::new(space, vec![space.vector([3, 2]), space.vector([0, -1])]);
            let ch3 = ch1.intersect(&ch2);
            println!("{ch3:?}");
        }
        //line misses line
        {
            let ch1 = ConvexHull::new(space, vec![space.vector([-2, -1]), space.vector([-1, 1])]);
            let ch2 = ConvexHull::new(space, vec![space.vector([1, 1]), space.vector([-1, 0])]);
            let ch3 = ch1.intersect(&ch2);
            debug_assert_eq!(ch3.subspace.embedded_space().affine_dimension(), 0);
        }
        //line hits line
        {
            let ch1 = ConvexHull::new(space, vec![space.vector([2, 0]), space.vector([-2, -1])]);
            let ch2 = ConvexHull::new(space, vec![space.vector([0, -1]), space.vector([1, 1])]);
            let ch3 = ch1.intersect(&ch2);
            debug_assert_eq!(ch3.subspace.embedded_space().affine_dimension(), 1);
        }
    }

    #[test]
    fn convex_hull_from_simplex() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 3);
        let p1 = space.vector([4, 2, 2]);
        let p2 = space.vector([5, -3, 3]);
        let p3 = space.vector([-5, 6, 2]);
        let p4 = space.vector([8, 2, -9]);

        ConvexHull::from_simplex(space.simplex(vec![]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(space.simplex(vec![p1.clone()]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(space.simplex(vec![p1.clone(), p2.clone()]).unwrap())
            .check()
            .unwrap();
        ConvexHull::from_simplex(
            space
                .simplex(vec![p1.clone(), p2.clone(), p3.clone()])
                .unwrap(),
        )
        .check()
        .unwrap();
        ConvexHull::from_simplex(
            space
                .simplex(vec![p1.clone(), p2.clone(), p3.clone(), p4.clone()])
                .unwrap(),
        )
        .check()
        .unwrap();
    }
}
