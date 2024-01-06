use std::collections::{HashMap, HashSet};

use itertools::Itertools;

use crate::geometry::vector::are_points_nondegenerage;

use super::{
    affine_coordinate_system::affine_span_of_simplices, oriented_simplex::OrientedSimplex,
    shape::Shape, simplex::Simplex, vector::Point,
};

#[derive(Debug, Clone)]
struct SimplicialComplexSimplexData {
    boundary_of: HashSet<Simplex>,
}

#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    dim: usize,
    simplices_and_data: HashMap<Simplex, SimplicialComplexSimplexData>,
}

impl SimplicialComplex {
    pub fn check(&self) -> Result<(), &'static str> {
        self.clone().as_shape().check()?; //simplices are disjoint

        for (simplex, data) in &self.simplices_and_data {
            for b in simplex.boundary().simplices_ref() {
                if !self.simplices_and_data.contains_key(b) {
                    return Err("simplicial complex is missing the boundary of some simplex");
                }
            }

            for inv_b in &data.boundary_of {
                if !self.simplices_and_data.contains_key(inv_b) {
                    return Err(
                        "inverse boundary contains a simplex not present in the simplicial complex",
                    );
                }

                if !inv_b.boundary_simplices().iter().any(|b| simplex == b) {
                    return Err(
                        "inverse boundary of simplex does not contain simplex in its boundary",
                    );
                }

                for (possible_b_inv, _) in &self.simplices_and_data {
                    if possible_b_inv
                        .boundary_simplices()
                        .iter()
                        .any(|b| simplex == b)
                    {
                        if !data.boundary_of.contains(possible_b_inv) {
                            return Err(
                                "inverse boundary of simplex is missing a simplex which contains simplex in its boundary",
                            );
                        }
                    }
                }
            }
        }

        Ok(())
    }

    pub fn empty(dim: usize) -> Self {
        Self {
            dim,
            simplices_and_data: HashMap::new(),
        }
    }

    pub fn new(dim: usize, simplices: Vec<Simplex>) -> Self {
        let mut simplices_and_data = HashMap::new();

        for simplex in &simplices {
            simplices_and_data.insert(
                simplex.clone(),
                SimplicialComplexSimplexData {
                    boundary_of: HashSet::new(),
                },
            );
        }

        for simplex in simplices {
            for bdry in simplex.boundary_simplices() {
                simplices_and_data
                    .get_mut(&bdry)
                    .unwrap()
                    .boundary_of
                    .insert(simplex.clone());
            }
        }

        let ans = Self {
            dim,
            simplices_and_data,
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn simplices(self) -> Vec<Simplex> {
        self.simplices_and_data
            .into_iter()
            .map(|(k, v)| k)
            .collect()
    }

    pub fn simplices_ref(&self) -> Vec<&Simplex> {
        self.simplices_and_data.iter().map(|(k, v)| k).collect()
    }

    pub fn as_shape(self) -> Shape {
        Shape::new(
            self.dim,
            self.simplices_and_data
                .into_iter()
                .map(|(k, v)| k)
                .collect(),
        )
    }

    pub fn is_empty(&self) -> bool {
        self.simplices_ref().is_empty()
    }

    pub fn point_boundary_of(&self, point: Point) -> &HashSet<Simplex> {
        self.simplex_boundary_of(Simplex::new(point.dim(), vec![point]))
    }

    pub fn simplex_boundary_of(&self, simplex: Simplex) -> &HashSet<Simplex> {
        &self.simplices_and_data.get(&simplex).unwrap().boundary_of
    }

    // fn subset(&self) -> SubSimplicialComplex {
    //     SubSimplicialComplex {
    //         simplicial_complex: self,
    //         subset: vec![&self.simplices[0]],
    //     }
    // }

    pub fn interior_and_boundary(&self) -> (SubSimplicialComplex, SubSimplicialComplex) {
        /*
        let n be the dimension of the space self is living in
         - every simplex of rank n is part of the interior
         - a simplex of rank n-1 is the facet of at most 2 simplices of rank n, and is part of the interior if and only if it is the facet of exactly 2 simplices of rank n
         - a simplex of rank less or equal to n-2 is part of the interior iff it is in the boundary of some strictly higher rank simplex AND every strictly higher rank simplex containing it as part of the boundary is part of the interior
        */

        //each simplex x points at all simplicies y such that x is part of the boundary of y
        let mut inverse_boundary_lookup: HashMap<&Simplex, HashSet<&Simplex>> = self
            .simplices_ref()
            .into_iter()
            .map(|s| (s, HashSet::new()))
            .collect();
        for s in self.simplices_ref() {
            for b in s.boundary_simplices() {
                inverse_boundary_lookup.get_mut(&b).unwrap().insert(s);
            }
        }

        let n = self.dim();

        let mut interior: HashSet<&Simplex> = HashSet::new();
        let mut boundary: HashSet<&Simplex> = HashSet::new();

        let mut all = self.simplices_ref();
        all.sort_by_key(|s| std::cmp::Reverse(s.rank().unwrap())); //so that we process largest rank first
        for simplex in all {
            let r = simplex.rank().unwrap();
            if r == n {
                //rank n simplex is always part of the interior
                interior.insert(simplex);
            } else if r == n - 1 {
                //rank n-1 simplex is part of the boundary of at most 2 rank n simplices
                //it is part of the boundary of the simplicial complex iff it is the boundary of exactly 2
                match inverse_boundary_lookup.get(simplex).unwrap().len() {
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
                let bdry = inverse_boundary_lookup.get(simplex).unwrap();
                if bdry.is_empty() {
                    boundary.insert(simplex);
                } else {
                    if bdry.iter().all(|b| interior.contains(*b)) {
                        interior.insert(simplex);
                    } else {
                        boundary.insert(simplex);
                    }
                }
            }
        }

        (
            SubSimplicialComplex {
                simplicial_complex: self,
                subset: interior,
            },
            SubSimplicialComplex {
                simplicial_complex: self,
                subset: boundary,
            },
        )
    }

    pub fn interior(&self) -> SubSimplicialComplex {
        self.interior_and_boundary().0
    }

    pub fn boundary(&self) -> SubSimplicialComplex {
        self.interior_and_boundary().1
    }

    pub fn simplify(mut self) -> Self {
        //go through each point
        //compute its star
        //enter the affine subspace spanned by its star
        //compute its link as oriented simplices
        //if point is part of the link then no simplification can be made, so more on
        //check whether any point of the link is in the interior with respect to every boundary of the link
        //if such a point exists, it can be used to fill in the star in a simpler way by fanning out

        'main: while true {
            for (simplex, data) in &self.simplices_and_data {
                if simplex.rank().unwrap() == 0 {
                    //loop through every point
                    let point = simplex.points().pop().unwrap();

                    //compute its star, link and nbd := star \sqcup link
                    let star = &data.boundary_of;
                    let mut nbd = HashSet::new();
                    nbd.insert(simplex.clone());
                    for s in star {
                        nbd.insert(s.clone());
                        for b in s.boundary_simplices() {
                            nbd.insert(b);
                        }
                    }
                    let nbd = nbd;
                    let mut link = nbd.clone();
                    for s in star {
                        link.remove(s);
                    }
                    link.remove(simplex);
                    let link = link;

                    let nbdsc =
                        SimplicialComplex::new(self.dim(), nbd.iter().map(|s| s.clone()).collect());

                    //compute a coordinate system for the affine subspace spanned by the star/link/nbd of the point
                    let nbd_coords =
                        affine_span_of_simplices(self.dim(), nbd.iter().collect()).unwrap();

                    let point_img = nbd_coords.point_preimage(&point).unwrap();
                    let simplex_img = nbd_coords.simplex_preimage(&simplex).unwrap();
                    simplex_img.check().unwrap();
                    let star_img: HashSet<_> = star
                        .iter()
                        .map(|s| nbd_coords.simplex_preimage(s).unwrap())
                        .collect();
                    let link_img: HashSet<_> = link
                        .iter()
                        .map(|s| nbd_coords.simplex_preimage(s).unwrap())
                        .collect();
                    let nbdsc_img = nbd_coords.simplicial_complex_preimage(&nbdsc).unwrap();
                    nbdsc_img.check().unwrap();

                    //proceed only if the point is locally an interior point
                    if nbdsc_img.interior().contains(&simplex_img) {
                        //compute orientations for the full rank simplicies of the link
                        let oriented_link_img = link_img
                            .iter()
                            .filter(|s| s.rank().unwrap() + 1 == nbd_coords.rank())
                            .map(|s| OrientedSimplex::from_simplex(s.clone(), &point_img))
                            .collect_vec();

                        //search for a link point which is on the interior with respect to every oriented facet of the link boundary
                        for link_point_img in link_img {
                            if link_point_img.rank().unwrap() == 0 {
                                let link_point_img = link_point_img.points().pop().unwrap();

                                if oriented_link_img.iter().all(|link_facet_img| {
                                    link_facet_img.sign_point(&link_point_img).is_le()
                                }) {
                                    //at this point we've found a point on the link which can be used as a replacement for point
                                    let new_point = nbd_coords.point_image(&link_point_img);
                                    println!("new_point = {:?}", new_point);

                                    let mut replacement_star_simplices = HashSet::new();
                                    for s in self.simplices_ref() {
                                        replacement_star_simplices.insert(s.clone());
                                    }
                                    replacement_star_simplices.remove(simplex);
                                    for s in star {
                                        replacement_star_simplices.remove(s);
                                    }
                                    for s in &link {
                                        if !s.points().contains(&new_point) {
                                            let mut i_points = s.points();
                                            i_points.push(new_point.clone());
                                            if are_points_nondegenerage(
                                                self.dim(),
                                                i_points.iter().collect(),
                                            ) {
                                                let i = Simplex::new_unsafe(self.dim(), i_points);
                                                replacement_star_simplices.insert(i);
                                            }
                                        }
                                    }

                                    self = Self::new(
                                        self.dim(),
                                        replacement_star_simplices.into_iter().collect(),
                                    );

                                    continue 'main;
                                }
                            }
                        }
                    }
                }
            }
            break;
        }
        self
    }
}

#[derive(Debug, Clone)]
pub struct SubSimplicialComplex<'a> {
    simplicial_complex: &'a SimplicialComplex,
    subset: HashSet<&'a Simplex>,
}

impl<'a> SubSimplicialComplex<'a> {
    pub fn check(&self) -> Result<(), &'static str> {
        self.simplicial_complex.check()?;

        for simplex in &self.subset {
            if !self
                .simplicial_complex
                .simplices_ref()
                .into_iter()
                .any(|sc_simplex| std::ptr::eq(*simplex, sc_simplex))
            {
                return Err("simplices in the subset should reference simplices in the simplicial complex only");
            }
        }

        Ok(())
    }

    pub fn dim(&self) -> usize {
        self.simplicial_complex.dim()
    }

    pub fn as_shape(&self) -> Shape {
        Shape::new(
            self.dim(),
            self.subset.iter().map(|s| (*s).clone()).collect(),
        )
    }

    pub fn contains(&self, simplex: &Simplex) -> bool {
        self.subset.contains(simplex)
    }
}
