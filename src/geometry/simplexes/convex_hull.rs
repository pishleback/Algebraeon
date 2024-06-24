use std::collections::{HashMap, HashSet};

use crate::geometry_old::simplicial_complex::SubSimplicialComplex;

use super::*;

#[derive(Debug, Clone)]
pub struct ConvexHull<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> where
    FS::Set: Hash,
{
    ambient_space: SP,
    subspace: EmbeddedAffineSubspace<FS, SP, Rc<AffineSpace<FS>>>,
    // oriented facets belonging to subspace such
    // the positive side of each facet is on the interior of the convex hull
    // the facets form a simplicial complex
    facets: Vec<OrientedSimplex<FS, Rc<AffineSpace<FS>>>>,
    interior: Vec<Simplex<FS, Rc<AffineSpace<FS>>>>,
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

#[derive(Debug, Clone)]
pub struct ConvexHullAsSimplicialComplexResult<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> {
    pub entire: Rc<SimplicialComplex<FS, SP>>,
    pub boundary: FullSubSimplicialComplex<FS, SP, Rc<SimplicialComplex<FS, SP>>>,
    pub interior: PartialSubSimplicialComplex<FS, SP, Rc<SimplicialComplex<FS, SP>>>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
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
                    != vec![Simplex::new(
                        self.subspace.embedded_space(),
                        vec![Vector::construct(
                            self.subspace.embedded_space(),
                            |i| unreachable!(),
                        )],
                    )
                    .unwrap()]
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
            if !ridges_count.into_iter().all(|(ridge, count)| count == 2) {
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
                            return Err("Every point must be on the positive or neutral side of every facet");
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

                let mut horizon = HashMap::new();
                for facet in &visible {
                    for ridge in facet.simplex().facets() {
                        match horizon.contains_key(&ridge) {
                            true => {
                                horizon.remove(&ridge);
                            }
                            false => {
                                horizon.insert(ridge, facet);
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
                        .into_iter()
                        .map(|(ridge, facet)| {
                            OrientedSimplex::new_with_positive_point(
                                self.subspace.embedded_space(),
                                {
                                    let mut points = ridge.points().clone();
                                    points.push(subsp_pt.clone());
                                    points
                                },
                                &{
                                    match facet.positive_point() {
                                        Some(pos_pt) => pos_pt,
                                        None => {
                                            debug_assert_eq!(
                                                self.subspace.embedded_space().affine_dimension(),
                                                1
                                            );
                                            /*
                                            unreachable because the embedded subspace has linear dimension 0
                                            therefore the only way we are here is if we start with a convex hull equal to a point and add that same point
                                            in this case hidden={the null facet} and visible={} so horizon={}
                                            this is contained in iteration over the horizon, so doesn't run
                                            */
                                            unreachable!()
                                        }
                                    }
                                },
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

    pub fn as_simplicial_complex(&self) -> ConvexHullAsSimplicialComplexResult<FS, SP> {
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
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        let mut boundary_indexes = HashSet::new();
        let mut interior_indexes = HashSet::new();
        for (idx, spx) in all_simplexes.iter().enumerate() {
            match boundary_simplexes.contains(spx) {
                true => {
                    boundary_indexes.insert(idx);
                }
                false => {
                    interior_indexes.insert(idx);
                }
            }
        }

        let entire = Rc::new(SimplicialComplex::new_unchecked(
            self.ambient_space.clone(),
            all_simplexes,
        ));
        ConvexHullAsSimplicialComplexResult {
            entire: entire.clone(),
            boundary: FullSubSimplicialComplex::new_unchecked(entire.clone(), boundary_indexes),
            interior: PartialSubSimplicialComplex::new_unchecked(entire.clone(), interior_indexes),
        }
    }

    pub fn intersect_with_oriented_hyperplane(
        &self,
        hyperplane: &OrientedHyperplane<FS, SP>,
        region: OrientationSide, //TODO: make this const generic once rust has const generic enums
    ) -> Self {
        //find outer_points and outer_edges giving the 1-skeleton of the surface of the convex hull
        let mut outer_points = HashSet::new();
        let mut outer_edges = HashSet::new();
        for facet in &self.facets {
            for point in facet.simplex().points() {
                outer_points.insert(point);
            }
            for edge in facet.simplex().edges() {
                outer_edges.insert(edge);
            }
        }

        /*
        Find positive_points, middle_points, negative_points such that
        ConvexHull(middle_points) = self intersect hyperplane
        ConvexHull(positive_points union middle_points) = self intersect (hyperplane union hyperplane_posside)
        ConvexHull(negative_points union middle_points) = self intersect (hyperplane union hyperplane_negside)
        */
        let mut positive_points = HashSet::new();
        let mut middle_points = HashSet::new();
        let mut negative_points = HashSet::new();
        for subsp_point in outer_points {
            let point = self.subspace.embed_point(subsp_point);
            match hyperplane.classify_point(&point) {
                OrientationSide::Positive => {
                    positive_points.insert(point);
                }
                OrientationSide::Neutral => {
                    middle_points.insert(point);
                }
                OrientationSide::Negative => {
                    negative_points.insert(point);
                }
            }
        }
        for edge in &outer_edges {
            debug_assert_eq!(edge.points().len(), 2);
            let a = self.subspace.embed_point(&edge.points()[0]);
            let b = self.subspace.embed_point(&edge.points()[1]);
            match hyperplane.intersect_line(&a, &b) {
                Some(m) => {
                    middle_points.insert(m);
                }
                None => {}
            }
        }

        Self::new(self.ambient_space.clone(), {
            match region {
                OrientationSide::Positive => middle_points
                    .union(&positive_points)
                    .map(|pt| pt.clone())
                    .collect(),
                OrientationSide::Neutral => {
                    middle_points.into_iter().map(|pt| pt.clone()).collect()
                }
                OrientationSide::Negative => middle_points
                    .union(&negative_points)
                    .map(|pt| pt.clone())
                    .collect(),
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use malachite_base::num::conversion::traits::FromSciString;
    use malachite_q::Rational;

    use crate::rings::structure::StructuredType;

    use super::*;

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
        let mut ch = ConvexHull::new(
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

        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Positive);

        let smaller_ch = ch.intersect_with_oriented_hyperplane(&ohsp, OrientationSide::Negative);
    }
}
