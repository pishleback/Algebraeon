use std::collections::{HashMap, HashSet};

use malachite_q::Rational;

use crate::{
    geometry::vector::Vector,
    rings::{
        lattice::{AffineLatticeElements, QQ_AFFLAT, QQ_LINLAT},
        matrix::{Matrix, QQ_MAT},
    },
};

use super::{
    oriented_simplex::OrientedSimplex, simplex::Simplex, simplicial_complex::SimplicialComplex,
};

//disjoint union of simplices
#[derive(Debug, Clone)]
pub struct Shape {
    dim: usize,
    simplices: Vec<Simplex>,
}

pub fn shape_disjoint_union(dim: usize, shapes: Vec<Shape>) -> Shape {
    let mut simplices = vec![];
    for mut shape in shapes {
        assert_eq!(shape.dim, dim);
        simplices.append(&mut shape.simplices);
    }
    let ans = Shape { dim, simplices };
    debug_assert!(ans.check().is_ok());
    ans
}

impl Shape {
    pub fn check(&self) -> Result<(), &'static str> {
        for simplex in &self.simplices {
            if !simplex.rank().is_some() {
                return Err("simplex in a shape musn't be the empty simplex");
            }

            if simplex.dim() != self.dim {
                return Err("Simplex dim does not match shape dim");
            }
            simplex.check()?;
        }
        //TODO: check that the simplices are disjoint
        Ok(())
    }

    pub fn new(dim: usize, simplices: Vec<Simplex>) -> Self {
        let ans = Self { dim, simplices };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn dim(&self) -> usize {
        self.dim
    }

    pub fn empty(dim: usize) -> Self {
        Self {
            dim,
            simplices: vec![],
        }
    }

    pub fn is_empty(&self) -> bool {
        self.simplices.is_empty()
    }

    pub fn simplex(simplex: Simplex) -> Self {
        Self {
            dim: simplex.dim(),
            simplices: vec![simplex],
        }
    }

    pub fn simplices(self) -> Vec<Simplex> {
        self.simplices
    }

    pub fn simplices_ref(&self) -> Vec<&Simplex> {
        self.simplices.iter().collect()
    }

    #[deprecated]
    fn transform(self, new_dim: usize, f: &dyn Fn(Vector) -> Vector) -> Self {
        Self::new(
            new_dim,
            self.simplices
                .into_iter()
                .map(|s| s.transform(new_dim, f))
                .collect(),
        )
    }

    fn points(&self) -> Vec<Vector> {
        let mut all_points = vec![];
        for s in &self.simplices {
            for p in s.points().iter() {
                all_points.push(p.clone());
            }
        }
        all_points
    }

    pub fn is_complete(&self) -> bool {
        let all_simplices: HashSet<_> = self.simplices.iter().collect();
        for s in &self.simplices {
            if s.n() >= 2 {
                for f in s.facets() {
                    if !all_simplices.contains(&f) {
                        return false;
                    }
                }
            }
        }
        true
    }

    pub fn completion(self) -> Shape {
        //does not always produce a valid shape
        let mut all_simplices = HashSet::new();
        for s in self.simplices {
            for b in s.boundary().simplices_ref() {
                all_simplices.insert(b.clone());
            }
            all_simplices.insert(s);
        }
        let ans = Shape {
            dim: self.dim,
            simplices: all_simplices.into_iter().collect(),
        };
        debug_assert!(ans.check().is_ok());
        ans
    }

    pub fn boundary(&self) -> Shape {
        let mut facets = HashSet::new();
        for simplex in &self.simplices {
            for facet in simplex.facets() {
                if facets.contains(&facet) {
                    facets.remove(&facet);
                } else {
                    facets.insert(facet);
                }
            }
        }
        let shape = Shape {
            dim: self.dim,
            simplices: facets.into_iter().collect(),
        };
        debug_assert!(shape.check().is_ok());
        shape
    }

    pub fn equal(&self, other: &Self) -> bool {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        self.symmetric_difference_nosimp(other).is_empty()
    }

    fn simplify(&self) -> Self {
        // 1) turn self into a subcomplex of a simplicial complex
        // 2) simplify
        todo!()
        // self.clone() //TODO
    }

    fn intersect_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        let mut parts = vec![];
        for s in &self.simplices {
            for t in &other.simplices {
                let (inner, _outer) = cut_simplex_by_simplex(&s, &t);
                parts.push(inner);
            }
        }
        shape_disjoint_union(dim, parts)
    }

    fn union_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        shape_disjoint_union(
            dim,
            vec![
                self.intersect_nosimp(other),
                self.symmetric_difference_nosimp(other),
            ],
        )
        // shape_union(dim, vec![self.subtract_nosimp(other), other.clone()])
    }

    fn subtract_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        let mut ans = self.clone();
        for s in &other.simplices {
            let (_inner, outer) = cut_shape_by_simplex(&s, &ans);
            ans = outer;
        }
        ans
    }

    fn symmetric_difference_nosimp(&self, other: &Self) -> Self {
        let dim = self.dim;
        assert_eq!(dim, other.dim);
        shape_disjoint_union(
            dim,
            vec![self.subtract_nosimp(other), other.subtract_nosimp(self)],
        )
    }

    pub fn intersect(&self, other: &Self) -> Self {
        self.intersect_nosimp(other).simplify()
    }
    pub fn union(&self, other: &Self) -> Self {
        self.union_nosimp(other).simplify()
    }
    pub fn subtract(&self, other: &Self) -> Self {
        self.subtract_nosimp(other).simplify()
    }
    pub fn symmetric_difference(&self, other: &Self) -> Self {
        self.symmetric_difference_nosimp(other).simplify()
    }
}

fn interior_of_convex_shell(shape: &Shape) -> Shape {
    debug_assert!(shape.is_complete());

    //shape should be the full boundary of a convex polytope

    if shape.simplices.is_empty() {
        return Shape::empty(shape.dim);
    }

    let n = shape.simplices.iter().map(|s| s.n()).max().unwrap();
    //n=1: shape is the shell of a line i.e. two points
    //n=2: shape is the shell of a polygon
    //n=3: shape is the shell of a solid
    //...

    if n == 1 {
        for s in &shape.simplices {
            debug_assert_eq!(s.n(), 1);
        }

        //every 1-dimensional convex shape must either be empty, a point, or a line
        match shape.simplices.len() {
            0 => Shape::empty(shape.dim), //no points has empty interior
            1 => Shape::empty(shape.dim), //single point has empty interior
            2 => Shape::simplex(Simplex::new(
                shape.dim,
                vec![
                    shape
                        .simplices_ref()
                        .iter()
                        .nth(0)
                        .unwrap()
                        .points()
                        .iter()
                        .nth(0)
                        .unwrap()
                        .clone(),
                    shape
                        .simplices_ref()
                        .iter()
                        .nth(1)
                        .unwrap()
                        .points()
                        .iter()
                        .nth(0)
                        .unwrap()
                        .clone(),
                ],
            )), //pair of points has interior of a line
            _ => {
                panic!()
            }
        }
    } else {
        debug_assert!(n >= 2);

        //n = 2: interior of a polygon
        //n = 3: interior of a solid
        //n = 4: interior of a 4d solid
        //n = 5: ...

        //the idea is to pick a vertex and add simplexes in a fan pattern from that vertex to disjoint bits of the boundary
        //e.g. for a pentagon: add 3 triangles, and 2 edges
        //e.g. for an icosahedron: add 15 tetrahedra, 25 triangles, and 6 edges

        if shape.simplices.len() == 0 {
            Shape::empty(shape.dim)
        } else {
            // println!("CUT SIMPLEX");

            //1) choose a base point
            let root = shape
                .simplices_ref()
                .into_iter()
                .nth(0)
                .unwrap()
                .points_ref()
                .into_iter()
                .nth(0)
                .unwrap();
            // println!("root = {:?}", root);
            // for s in &shape.simplices {
            //     println!("shell: {:?}", s);
            // }

            //2) find all cells adjacent to root
            let adj_cell_spans: Vec<_> = shape
                .simplices
                .iter()
                .filter(|s| s.has_vertex(&root))
                .map(|s| s.affine_span())
                .collect();

            //3) find which adjacent faces are vertex belongs to
            let mut point_degen = HashMap::new();
            for s in &shape.simplices {
                if s.n() == 1 {
                    let pt = s.points_ref().into_iter().nth(0).unwrap();
                    debug_assert!(!point_degen.contains_key(pt));
                    let adj_cells: HashSet<usize> = adj_cell_spans
                        .iter()
                        .enumerate()
                        .filter(|(_idx, adj_cell_span)| {
                            QQ_AFFLAT.contains_point(adj_cell_span, pt.as_matrix())
                        })
                        .map(|(idx, _adj_cell_span)| idx)
                        .collect();
                    point_degen.insert(pt, adj_cells);
                }
            }
            //check that every vertex of every part of the shell is present - it should be since the shell is complete
            for s in &shape.simplices {
                for p in s.points_ref() {
                    debug_assert!(point_degen.contains_key(p));
                }
            }
            // println!("{:?}", point_degen);

            //4) for each shell simplex, if not all vertices lie in some adjacent face span, fill it in
            let mut interior = Shape::empty(shape.dim);
            for s in &shape.simplices {
                debug_assert!(s.points().len() >= 1);
                let mut common = point_degen[&s.points().iter().nth(0).unwrap()].clone();
                common = common
                    .into_iter()
                    .filter(|idx| {
                        (1..s.points().len())
                            .map(|i| s.points_ref()[i])
                            .all(|pt| point_degen[pt].contains(idx))
                    })
                    .collect();
                if common.len() == 0 {
                    let filler = Simplex::new(shape.dim, {
                        let mut filler_pts = vec![root.clone()];
                        filler_pts.append(&mut s.points().clone());
                        // println!("filler_pts = {:?}", filler_pts);
                        filler_pts
                    });
                    // println!("filler simplex: {:?}", filler);
                    interior =
                        shape_disjoint_union(shape.dim, vec![interior, Shape::simplex(filler)]);
                }
            }
            interior
        }
    }
}

fn cut_simplex_by_plane(cut_plane: &OrientedSimplex, simplex: &Simplex) -> (Shape, Shape, Shape) {
    let dim = cut_plane.dim();
    assert_eq!(dim, simplex.dim());
    debug_assert!(simplex.n() <= dim + 1);

    match simplex.n() {
        0 => (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim)),
        1 => match cut_plane.sign_point(&simplex.points().iter().nth(0).unwrap()) {
            std::cmp::Ordering::Greater => (
                Shape::simplex(simplex.clone()),
                Shape::empty(dim),
                Shape::empty(dim),
            ),
            std::cmp::Ordering::Equal => (
                Shape::empty(dim),
                Shape::simplex(simplex.clone()),
                Shape::empty(dim),
            ),
            std::cmp::Ordering::Less => (
                Shape::empty(dim),
                Shape::empty(dim),
                Shape::simplex(simplex.clone()),
            ),
        },
        2 => {
            let (p, q) = (
                simplex.points_ref().into_iter().nth(0).unwrap(),
                simplex.points_ref().into_iter().nth(1).unwrap(),
            );
            match (cut_plane.sign_point(p), cut_plane.sign_point(q)) {
                (std::cmp::Ordering::Greater, std::cmp::Ordering::Greater)
                | (std::cmp::Ordering::Equal, std::cmp::Ordering::Greater)
                | (std::cmp::Ordering::Greater, std::cmp::Ordering::Equal) => (
                    Shape::simplex(simplex.clone()),
                    Shape::empty(dim),
                    Shape::empty(dim),
                ),
                (std::cmp::Ordering::Equal, std::cmp::Ordering::Equal) => (
                    Shape::empty(dim),
                    Shape::simplex(simplex.clone()),
                    Shape::empty(dim),
                ),
                (std::cmp::Ordering::Less, std::cmp::Ordering::Less)
                | (std::cmp::Ordering::Equal, std::cmp::Ordering::Less)
                | (std::cmp::Ordering::Less, std::cmp::Ordering::Equal) => (
                    Shape::empty(dim),
                    Shape::empty(dim),
                    Shape::simplex(simplex.clone()),
                ),
                (p_sign, q_sign) => {
                    let t = -cut_plane.det_point(&p) / cut_plane.det_vector(&(q - p));
                    let r = p + &(&t * &(q - p));
                    debug_assert!(0 < t && t < 1);
                    debug_assert_eq!(cut_plane.det_point(&r), Rational::from(0));
                    match (p_sign, q_sign) {
                        (std::cmp::Ordering::Greater, std::cmp::Ordering::Less) => (
                            Shape::simplex(Simplex::new(dim, vec![p.clone(), r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r, q.clone()])),
                        ),
                        (std::cmp::Ordering::Less, std::cmp::Ordering::Greater) => (
                            Shape::simplex(Simplex::new(dim, vec![q.clone(), r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r.clone()])),
                            Shape::simplex(Simplex::new(dim, vec![r, p.clone()])),
                        ),
                        _ => panic!(),
                    }
                }
            }
        }
        n => {
            debug_assert!(n >= 3);
            let (open_left, hollow_middle, open_right) =
                cut_shape_by_plane(&cut_plane, &simplex.boundary().as_shape());

            let interior_middle = interior_of_convex_shell(&hollow_middle);
            let hollow_left = shape_disjoint_union(
                dim,
                vec![
                    open_left.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let hollow_right = shape_disjoint_union(
                dim,
                vec![
                    open_right.clone(),
                    hollow_middle.clone(),
                    interior_middle.clone(),
                ],
            );
            let interior_left = interior_of_convex_shell(&hollow_left);
            let interior_right = interior_of_convex_shell(&hollow_right);
            (interior_left, interior_middle, interior_right)
        }
    }
}

fn cut_shape_by_plane(cut_plane: &OrientedSimplex, shape: &Shape) -> (Shape, Shape, Shape) {
    let dim = cut_plane.dim();
    assert_eq!(dim, shape.dim);
    let (mut left, mut middle, mut right) =
        (Shape::empty(dim), Shape::empty(dim), Shape::empty(dim));
    for simplex in &shape.simplices {
        let (s_left, s_middle, s_right) = cut_simplex_by_plane(&cut_plane, &simplex);
        left = shape_disjoint_union(dim, vec![left, s_left]);
        middle = shape_disjoint_union(dim, vec![middle, s_middle]);
        right = shape_disjoint_union(dim, vec![right, s_right]);
    }
    (left, middle, right)
}

//return (part of simplex inside cut_simplex, part of simplex outside cut_simplex)
pub fn cut_simplex_by_simplex(cut_simplex: &Simplex, simplex: &Simplex) -> (Shape, Shape) {
    let dim = cut_simplex.dim();
    assert!(dim >= 1);
    assert_eq!(dim, simplex.dim());

    if !cut_simplex
        .bounding_box()
        .intersects(&simplex.bounding_box())
    {
        return (Shape::empty(dim), Shape::simplex(simplex.clone()));
    }

    //cut in some hyperplanes to restrict to the affine subspace of cut_simplex
    match cut_simplex.affine_span().elems() {
        AffineLatticeElements::Empty() => panic!(),
        AffineLatticeElements::NonEmpty { offset, linlat } => {
            let flat_dim = QQ_LINLAT.rank(&linlat);
            let offset_pt = Vector::from_matrix(&offset);
            let offset_vec = Vector::from_matrix(&offset);

            let linhyperplanes = QQ_LINLAT.as_hyperplane_intersection(&linlat);

            let mut inside = Shape::simplex(simplex.clone());
            let mut outside = Shape::empty(dim);

            for linhyperplane in linhyperplanes {
                //construct an oriented simplex which spans linhyperplane
                let mut pts = vec![];
                pts.push(offset_pt.clone());
                for c in 0..dim - 1 {
                    pts.push(
                        &offset_pt
                            + &Vector::from_matrix(&QQ_LINLAT.basis_matrix(&linhyperplane, c)),
                    );
                }

                let cut_plane = OrientedSimplex::new(Simplex::new(dim, pts));
                debug_assert!(cut_plane.check().is_ok());

                let (inside_left, inside_mid, inside_right) =
                    cut_shape_by_plane(&cut_plane, &inside);
                // let (outside_left, outside_mid, outside_right) =
                //     cut_shape_by_plane(&cut_plane, &outside);

                inside = inside_mid;
                outside = shape_disjoint_union(dim, vec![inside_left, inside_right, outside]);
            }

            //restrict to the span of cut_simplex and intersect in there using the facets as oriented hyperplanes in the subspace
            let mat = Matrix::join_cols(
                dim,
                (0..flat_dim)
                    .map(|i| QQ_LINLAT.basis_matrix(&linlat, i))
                    .collect(),
            );

            if flat_dim != 0 {
                //apply so that cut simplex has full rank in the affine subspace
                // mat(pt - offset)
                //apply to put back into this space
                // (mat(!Ǝ)=x) + offset

                let flat_transform = &|pt: Vector| {
                    Vector::from_matrix(
                        &QQ_MAT
                            .col_solve(&mat, &(&pt - &offset_vec).as_matrix())
                            .unwrap(),
                    )
                };
                let flat_cut_simplex = cut_simplex.clone().transform(flat_dim, flat_transform);
                let mut flat_inside = inside.transform(flat_dim, flat_transform);
                let mut flat_outside = Shape::empty(flat_dim);
                for k in 0..flat_cut_simplex.n() {
                    let (flat_inside_out, flat_inside_bdry, flat_inside_in) =
                        cut_shape_by_plane(&flat_cut_simplex.oriented_facet(k), &flat_inside);
                    // let (flat_outside_out, flat_outside_bdry, flat_outside_in) =
                    //     cut_shape_by_plane(&simplex.oriented_facet(k), &flat_outside);

                    flat_inside = flat_inside_in;
                    flat_outside = shape_disjoint_union(
                        flat_dim,
                        vec![flat_inside_out, flat_inside_bdry, flat_outside],
                    );
                }

                let unflat_transform = &|pt: Vector| {
                    &Vector::from_matrix(&QQ_MAT.mul_refs(&mat, &pt.as_matrix()).unwrap())
                        + &offset_vec
                };
                let unflat_inside = flat_inside.transform(dim, unflat_transform);
                let unflat_outside = flat_outside.transform(dim, unflat_transform);

                (inside, outside) = (
                    unflat_inside,
                    shape_disjoint_union(dim, vec![unflat_outside, outside]),
                );
            }

            if inside.is_empty() {
                (inside, Shape::simplex(simplex.clone()))
            } else if outside.is_empty() {
                (Shape::simplex(simplex.clone()), outside)
            } else {
                (inside, outside)
            }
        }
    }
}

//return (part of shape inside cut_simplex, part of shape outside cut_simplex)
pub fn cut_shape_by_simplex(cut_simplex: &Simplex, shape: &Shape) -> (Shape, Shape) {
    let dim = cut_simplex.dim();
    assert_eq!(dim, shape.dim);
    let (mut inside, mut outside) = (Shape::empty(dim), Shape::empty(dim));
    for simplex in &shape.simplices {
        let (s_inside, s_outside) = cut_simplex_by_simplex(&cut_simplex, &simplex);
        inside = shape_disjoint_union(dim, vec![inside, s_inside]);
        outside = shape_disjoint_union(dim, vec![outside, s_outside]);
    }
    (inside, outside)
}

/*
pub fn convexhull_boundary(dim: usize, points: Vec<Point>) -> (Shape, usize) {
    //return None if the convex hull is flat
    for point in &points {
        assert_eq!(point.dim(), dim);
    }

    if points.len() == 0 {
        panic!("points should be non-empty");
    } else if points.len() == 1 {
        (Shape::simplex(Simplex::new(dim, points)), 0)
    } else {
        let (first, rest) = points.split_first().unwrap();

        let mat = Matrix::construct(dim, rest.len(), |r, c| {
            rest[c].get_coord(r) - first.get_coord(r)
        });
        let (_h, _u, _u_det, pivs) = QQ_MAT.row_hermite_algorithm(mat);
        let rank = pivs.len();

        let starting_simplex = Simplex::new(dim, {
            let mut starting_simplex_points = vec![first.clone()];
            for i in pivs {
                starting_simplex_points.push(rest[i].clone());
            }
            starting_simplex_points
        });
        debug_assert!(starting_simplex.check().is_ok());

        debug_assert!(rank <= dim);
        if rank == dim {
            let middle_point = starting_simplex.centroid();
            let mut facets: HashSet<_> = starting_simplex.oriented_facets().into_iter().collect();
            for pt in points.iter() {
                let visible = facets
                    .iter()
                    .filter(|f| f.sign_point(pt).is_ge())
                    .map(|f| f.clone())
                    .collect_vec();

                let mut horizon = HashSet::new();
                for v in visible {
                    for b in v.simplex.facets() {
                        if horizon.contains(&b) {
                            horizon.remove(&b);
                        } else {
                            horizon.insert(b);
                        }
                    }
                    facets.remove(&v);
                }

                for h in horizon {
                    let mut new_facet_points = h.vertices;
                    new_facet_points.push(pt.clone());
                    let new_facet =
                        OrientedSimplex::from_points(dim, new_facet_points, &middle_point);
                    new_facet.check().unwrap();
                    facets.insert(new_facet);
                }
            }

            (
                Shape {
                    dim,
                    simplices: facets.into_iter().map(|facet| facet.simplex).collect(),
                }
                .completion(),
                rank,
            )
        } else {
            match starting_simplex.affine_span().elems() {
                AffineLatticeElements::Empty() => panic!(),
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    let flat_dim = QQ_LINLAT.rank(linlat);
                    let offset_pt = Point::from_matrix(offset);
                    let offset_vec = Vector::from_matrix(offset);

                    //restrict to the span of starting_simplex and take convex hull there
                    let mat = Matrix::join_cols(
                        dim,
                        (0..flat_dim)
                            .map(|i| QQ_LINLAT.basis_matrix(linlat, i))
                            .collect(),
                    );

                    if flat_dim == 0 {
                        panic!("convex hull of a single point should already be taken care of");
                    } else {
                        //apply so that cut simplex has full rank in the affine subspace
                        // mat(pt - offset)
                        //apply to put back into this space
                        // (mat(!Ǝ)=x) + offset

                        let flat_transform = &|pt: Point| {
                            Point::from_matrix(
                                &QQ_MAT
                                    .col_solve(&mat, &(&pt - &offset_vec).as_matrix())
                                    .unwrap(),
                            )
                        };
                        let flat_points = points
                            .iter()
                            .map(|pt| pt.clone().transform(flat_transform))
                            .collect();
                        let (flat_hull_boundary, rank) = convexhull_boundary(flat_dim, flat_points);
                        debug_assert_eq!(rank, flat_dim);

                        let unflat_transform = &|pt: Point| {
                            &Point::from_matrix(&QQ_MAT.mul_refs(&mat, &pt.as_matrix()).unwrap())
                                + &offset_vec
                        };
                        let unflat_hull_boundary =
                            flat_hull_boundary.transform(dim, unflat_transform);

                        (unflat_hull_boundary, flat_dim)
                    }
                }
            }
        }
    }
}

pub fn convexhull_interior(dim: usize, points: Vec<Point>) -> (Shape, usize) {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    let (shell, rank) = convexhull_boundary(dim, points);
    (interior_of_convex_shell(&shell), rank)
}

pub fn convexhull(dim: usize, points: Vec<Point>) -> Shape {
    for point in &points {
        assert_eq!(point.dim(), dim);
    }
    let (shell, _rank) = convexhull_boundary(dim, points);
    let interior = interior_of_convex_shell(&shell);
    shape_disjoint_union(dim, vec![shell, interior])
}
*/

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_cut_simplex() {
        let s = Simplex::new(
            3,
            vec![
                Vector::new(vec![
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("1").unwrap(),
                ]),
                Vector::new(vec![
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("-1").unwrap(),
                ]),
                Vector::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("2").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Vector::new(vec![
                    Rational::from_str("2").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
            ],
        );

        let h = OrientedSimplex::new(Simplex::new(
            3,
            vec![
                Vector::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Vector::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                ]),
                Vector::new(vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("1").unwrap(),
                ]),
            ],
        ));

        s.check().unwrap();
        h.check().unwrap();
        let (a, b, c) = cut_simplex_by_plane(&h, &s);
        a.check().unwrap();
        b.check().unwrap();
        c.check().unwrap();
        println!(
            "{:?}",
            a.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        println!(
            "{:?}",
            b.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        println!(
            "{:?}",
            c.clone()
                .simplices()
                .iter()
                .map(|s| s.n())
                .collect::<Vec<_>>()
        );
        let d = shape_disjoint_union(3, vec![a, b, c]);
        d.check().unwrap();
    }

    #[test]
    fn shape_equal() {
        let shape1 = Shape {
            dim: 2,
            simplices: vec![
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
            ],
        };

        let shape2 = Shape {
            dim: 2,
            simplices: vec![
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
                Simplex::new(
                    2,
                    vec![
                        Vector::new(vec![
                            Rational::from_str("0").unwrap(),
                            Rational::from_str("0").unwrap(),
                        ]),
                        Vector::new(vec![
                            Rational::from_str("1").unwrap(),
                            Rational::from_str("1").unwrap(),
                        ]),
                    ],
                ),
            ],
        };

        assert!(shape1.equal(&shape2));
    }

    // #[test]
    // fn test_simplex_boundary_as_mesh() {
    //     let a = Simplex::new(
    //         2,
    //         vec![
    //             Point::new(vec![
    //                 Rational::from_str("0").unwrap(),
    //                 Rational::from_str("1").unwrap(),
    //             ]),
    //             Point::new(vec![
    //                 Rational::from_str("2").unwrap(),
    //                 Rational::from_str("1").unwrap(),
    //             ]),
    //             Point::new(vec![
    //                 Rational::from_str("1").unwrap(),
    //                 Rational::from_str("-1").unwrap(),
    //             ]),
    //         ],
    //     );

    //     debug_assert!(a.boundary_as_mesh().unwrap().check().is_ok());
    // }
}
