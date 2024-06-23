use super::*;

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> PartialOrd
    for Vector<FS, SP>
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let space = common_space(
            self.ambient_space().borrow(),
            other.ambient_space().borrow(),
        )?;
        for i in 0..space.linear_dimension().unwrap() {
            match space
                .ordered_field()
                .ring_cmp(self.coordinate(i), other.coordinate(i))
            {
                std::cmp::Ordering::Less => {
                    return Some(std::cmp::Ordering::Less);
                }
                std::cmp::Ordering::Equal => {}
                std::cmp::Ordering::Greater => {
                    return Some(std::cmp::Ordering::Greater);
                }
            }
        }
        Some(std::cmp::Ordering::Equal)
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Ord
    for Vector<FS, SP>
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.partial_cmp(other) {
            Some(ans) => ans,
            None => panic!(),
        }
    }
}

mod simplex {
    use itertools::Itertools;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct Simplex<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
    > {
        ambient_space: SP,
        // points must be ordered w.r.t the ordering on vectors
        // points must be non-degenerate
        // points must belong to the ambient_space
        points: Vec<Vector<FS, SP>>,
    }

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> PartialEq
        for Simplex<FS, SP>
    {
        fn eq(&self, other: &Self) -> bool {
            self.ambient_space.borrow() == other.ambient_space.borrow()
                && self.points == other.points
        }
    }

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Eq
        for Simplex<FS, SP>
    {
    }

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Hash
        for Simplex<FS, SP>
    where
        FS::Set: Hash,
    {
        fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
            self.points.hash(state);
        }
    }

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Simplex<FS, SP>
    where
        SP: Clone,
    {
        pub fn new(
            ambient_space: SP,
            mut points: Vec<Vector<FS, SP>>,
        ) -> Result<Self, &'static str> {
            for point in &points {
                assert_eq!(ambient_space.borrow(), point.ambient_space().borrow());
            }
            points.sort_unstable();
            if !ambient_space
                .borrow()
                .are_points_affine_independent(points.iter().collect())
            {
                Err("Can't make a simplex using degenerate points")
            } else {
                Ok(Self {
                    ambient_space,
                    points,
                })
            }
        }

        pub fn ambient_space(&self) -> SP {
            self.ambient_space.clone()
        }

        pub fn n(&self) -> usize {
            self.points.len()
        }

        pub fn points(&self) -> &Vec<Vector<FS, SP>> {
            &self.points
        }

        pub fn skeleton(&self, skel_n: isize) -> Vec<Self> {
            if skel_n < 0 {
                vec![]
            } else {
                let skel_n = skel_n as usize;
                let mut parts = vec![];
                for subset in (0..self.points.len()).combinations(skel_n) {
                    let part = Self::new(
                        self.ambient_space.clone(),
                        subset.into_iter().map(|i| self.points[i].clone()).collect(),
                    )
                    .unwrap();
                    parts.push(part);
                }
                parts
            }
        }

        pub fn vertices(&self) -> Vec<Self> {
            self.skeleton(1)
        }

        pub fn edges(&self) -> Vec<Self> {
            self.skeleton(2)
        }

        pub fn faces(&self) -> Vec<Self> {
            self.skeleton(3)
        }

        pub fn ridges(&self) -> Vec<Self> {
            self.skeleton(self.points.len() as isize - 2)
        }

        pub fn facets(&self) -> Vec<Self> {
            self.skeleton(self.points.len() as isize - 1)
        }

        pub fn facet(&self, k: usize) -> Self {
            assert!(k <= self.points.len());
            let facet = Self::new(self.ambient_space.clone(), {
                let mut facet_points = self.points.clone();
                facet_points.remove(k);
                facet_points
            })
            .unwrap();
            facet
        }

        pub fn sub_simplices(&self) -> Vec<Self> {
            self.points()
                .clone()
                .into_iter()
                .powerset()
                .map(|sub_points| Self::new(self.ambient_space.clone(), sub_points).unwrap())
                .collect()
        }

        pub fn oriented_facet(&self, k: usize) -> OrientedSimplex<FS, SP> {
            //return the oriented facet of self with positive side on the outside and negative side on the inside
            assert!(k <= self.points.len());
            let mut facet_points = self.points.clone();
            let other_pt = facet_points.remove(k);
            OrientedSimplex::new_with_negative_point(
                self.ambient_space.clone(),
                facet_points,
                &other_pt,
            )
            .unwrap()
        }

        pub fn oriented_facets(&self) -> Vec<OrientedSimplex<FS, SP>> {
            assert_eq!(self.ambient_space.borrow().affine_dimension(), self.n());
            (0..self.n()).map(|k| self.oriented_facet(k)).collect()
        }
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum OrientationSide {
        Positive,
        Neutral,
        Negative,
    }

    /*
    a simplex spanning a hyperplane with a positive side and a negative side
    in 3d space it is a triangle
    in 2d space it is a line
    in 1d space it is a point
    in 0d space it is a null simplex. The orientation looses meaning but it is nice to still count this case.
    */

    #[derive(Debug, Clone)]
    struct OrientedSimplexOrientation<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
    > {
        flip: bool,                     // flip the orientation if necessary
        positive_point: Vector<FS, SP>, // a point on the positive side
    }

    #[derive(Debug, Clone)]
    pub struct OrientedSimplex<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
    > {
        simplex: Simplex<FS, SP>,
        orientation: Option<OrientedSimplexOrientation<FS, SP>>, //None iff simplex is null
    }

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
        OrientedSimplex<FS, SP>
    {
        pub fn new_with_positive_point(
            ambient_space: SP,
            mut points: Vec<Vector<FS, SP>>,
            ref_point: &Vector<FS, SP>,
        ) -> Result<Self, &'static str> {
            assert_eq!(ref_point.ambient_space().borrow(), ambient_space.borrow());
            if points.len() != ambient_space.borrow().linear_dimension().unwrap() {
                return Err("Oriented simplex must have dimension one less than the ambient space");
            }
            let n = points.len();
            let simplex = Simplex::new(ambient_space, points)?;
            if n == 0 {
                Ok(Self {
                    simplex,
                    orientation: None,
                })
            } else {
                let mut guess = Self {
                    simplex,
                    orientation: Some(OrientedSimplexOrientation {
                        flip: false,
                        positive_point: ref_point.clone(),
                    }),
                };
                match guess.classify_point(ref_point) {
                    OrientationSide::Positive => Ok(guess),
                    OrientationSide::Neutral => Err("The reference point lies on the hyperplane"),
                    OrientationSide::Negative => {
                        let orientation = guess.orientation.as_mut().unwrap();
                        orientation.flip = !orientation.flip;
                        Ok(guess)
                    }
                }
            }
        }

        pub fn new_with_negative_point(
            ambient_space: SP,
            mut points: Vec<Vector<FS, SP>>,
            ref_point: &Vector<FS, SP>,
        ) -> Result<Self, &'static str> {
            Ok(Self::new_with_positive_point(ambient_space, points, ref_point)?.flip())
        }

        pub fn positive_point(&self) -> Option<Vector<FS, SP>> {
            Some(self.orientation.as_ref()?.positive_point.clone())
        }

        pub fn negative_point(&self) -> Option<Vector<FS, SP>> {
            let positive_point = self.positive_point()?;
            Some({
                let pt: &Vector<FS, SP> = &self.simplex.points()[0];
                &(pt + pt) - &positive_point
            })
        }

        pub fn ambient_space(&self) -> SP {
            self.simplex().ambient_space()
        }

        pub fn simplex(&self) -> &Simplex<FS, SP> {
            &self.simplex
        }

        pub fn flip(mut self) -> Self {
            let negative_point = self.negative_point();
            match &mut self.orientation {
                Some(OrientedSimplexOrientation {
                    flip,
                    positive_point,
                }) => {
                    (*flip, *positive_point) = (!*flip, negative_point.unwrap());
                }
                None => {}
            }
            match &self.orientation {
                Some(OrientedSimplexOrientation {
                    flip,
                    positive_point,
                }) => {
                    debug_assert_eq!(
                        self.classify_point(positive_point),
                        OrientationSide::Positive
                    );
                }
                None => {}
            }
            self
        }

        pub fn classify_point(&self, point: &Vector<FS, SP>) -> OrientationSide {
            let space = self.ambient_space();
            let field = space.borrow().ordered_field();
            match &self.orientation {
                Some(OrientedSimplexOrientation {
                    flip,
                    positive_point,
                }) => match space.borrow().linear_dimension().unwrap() {
                    0 => unreachable!(),
                    d => {
                        let root = &self.simplex.points[0];
                        let mut vecs = (1..d)
                            .map(|i| &self.simplex.points[i] - root)
                            .collect::<Vec<_>>();
                        vecs.push(point - root);
                        let det = space.borrow().determinant(vecs.iter().collect());
                        match flip {
                            false => match field.ring_cmp(&det, &field.zero()) {
                                std::cmp::Ordering::Less => OrientationSide::Negative,
                                std::cmp::Ordering::Equal => OrientationSide::Neutral,
                                std::cmp::Ordering::Greater => OrientationSide::Positive,
                            },
                            true => match field.ring_cmp(&det, &field.zero()) {
                                std::cmp::Ordering::Less => OrientationSide::Positive,
                                std::cmp::Ordering::Equal => OrientationSide::Neutral,
                                std::cmp::Ordering::Greater => OrientationSide::Negative,
                            },
                        }
                    }
                },
                None => OrientationSide::Neutral,
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use malachite_q::Rational;

        use crate::rings::structure::StructuredType;

        use super::*;

        #[test]
        fn make_simplex() {
            let space = AffineSpace::new_linear(Rational::structure(), 2);
            let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(1)]);
            let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
            let v3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
            let s = Simplex::new(&space, vec![v1, v2, v3]);
            assert!(s.is_ok());

            let space = AffineSpace::new_linear(Rational::structure(), 2);
            let v1 = Vector::new(&space, vec![Rational::from(0), Rational::from(0)]);
            let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
            let v3 = Vector::new(&space, vec![Rational::from(2), Rational::from(0)]);
            let s = Simplex::new(&space, vec![v1, v2, v3]);
            assert!(s.is_err());
        }

        #[test]
        fn simplex_skeleton() {
            let space = AffineSpace::new_linear(Rational::structure(), 2);
            let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(1)]);
            let v2 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
            let v3 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
            let s = Simplex::new(&space, vec![v1, v2, v3]).unwrap();

            assert_eq!(s.skeleton(-2).len(), 0);
            assert_eq!(s.skeleton(-1).len(), 0);
            assert_eq!(s.skeleton(0).len(), 1);
            assert_eq!(s.vertices().len(), 3);
            assert_eq!(s.edges().len(), 3);
            assert_eq!(s.faces().len(), 1);
            assert_eq!(s.skeleton(4).len(), 0);
        }

        #[test]
        fn make_oriented_simplex() {
            let space = AffineSpace::new_linear(Rational::structure(), 2);
            let v1 = Vector::new(&space, vec![Rational::from(1), Rational::from(0)]);
            let v2 = Vector::new(&space, vec![Rational::from(0), Rational::from(1)]);
            let v3 = Vector::new(&space, vec![Rational::from(2), Rational::from(3)]);
            let s_pos =
                OrientedSimplex::new_with_positive_point(&space, vec![v1.clone(), v2.clone()], &v3)
                    .unwrap();
            let s_neg =
                OrientedSimplex::new_with_negative_point(&space, vec![v1, v2], &v3).unwrap();

            assert_ne!(
                s_pos.orientation.unwrap().flip,
                s_neg.orientation.unwrap().flip
            );
        }
    }
}
pub use simplex::*;

mod simplicial_complex {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct SimplicialComplex {}

    #[derive(Debug, Clone)]
    pub struct SubSimplicialComplex<SC: Borrow<SimplicialComplex>> {
        simplicial_complex: SC,
    }
}
pub use simplicial_complex::*;

mod partial_simplicial_complex {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct PartialSimplicialComplex {}

    impl PartialSimplicialComplex {
        pub fn overlap_components(
            a: &PartialSimplicialComplex,
            b: &PartialSimplicialComplex,
        ) -> (
            Rc<SimplicialComplex>,                       // a union b
            SubSimplicialComplex<Rc<SimplicialComplex>>, // a intersect b
            SubSimplicialComplex<Rc<SimplicialComplex>>, // a without b
            SubSimplicialComplex<Rc<SimplicialComplex>>, // b without a
        ) {
            todo!()
        }
    }
}
pub use partial_simplicial_complex::*;

mod convex_hull {
    use std::collections::{HashMap, HashSet};

    use super::*;
    #[derive(Debug, Clone)]
    pub struct ConvexHull<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
    >
    where
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
                        return Err(
                            "Facets should be non-empty whenenver the subspace is non-empty",
                        );
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
                    return Err(
                        "Ridges of facets should each be shared between exactly two facets",
                    );
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

        fn extend_by_point(&mut self, pt: Vector<FS, SP>)
        where
            SP: std::fmt::Debug,
        {
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
                                                    self.subspace
                                                        .embedded_space()
                                                        .affine_dimension(),
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
            //add a point to extend by three triangles
            ch.extend_by_point(Vector::new(
                &space,
                vec![Rational::from(2), Rational::from(2)],
            ));
        }
    }
}
pub use convex_hull::*;

// fn split<
//     FS: OrderedRingStructure + FieldStructure,
//     SP: Borrow<LinearSpace<FS>> + Clone,
// >(
//     a: &Simplex<FS, SP>,
//     b: &OrientedSimplex<FS, SP>,
// ) -> (
//     Vec<Simplex<FS, SP>>,
//     Vec<Simplex<FS, SP>>,
//     Vec<Simplex<FS, SP>>,
// ) {
//     todo!()
// }

// fn intersect<
//     FS: OrderedRingStructure + FieldStructure,
//     SP: Borrow<LinearSpace<FS>> + Clone,
// >(
//     a: &Simplex<FS, SP>,
//     b: &Simplex<FS, SP>,
// ) -> Vec<Simplex<FS, SP>> {
//     todo!()
// }

/*
Plan

TODO: SC venn SC:
    For each simplex in one paired with a simplex in the other subdivide each and subdivide the SCs as induced such that we have a venn diagram

TODO: Refine SC given refinement of some simplex:
    Delete the interior of the simplex
    Get induced refinement of all facets, ridges, ...
    Replace facets, riges, ... with their refinement and induce refinement of adjacent simplexes in the SC
    Reinsert the refined interior of the simplex

TODO: Intersect a pair of simplexes A, B
    Express the subspace spanned by B as the intersection of hyperplanes
    Intersect A other with each hyperplane
    Restrict to the affine subspace spanned by B
    Intersect A with the interior side of each oriented facet of B
    Return to the ambient space
    Compute the convex hull as an SC

TODO: Intersect a simplex with an oriented hyperplane
    Compute the convex hull of any neutral points and any edges with one positive and one negative endpoint

TODO: Get the component of a simplex on the positive side of a hyperplane
    induct on the dimension of the simplex
    0d: the simplex is a point and it is trivial
    nd: take the interior of the convex hull of the positive component of each facet

*/
