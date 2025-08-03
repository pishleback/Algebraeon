use super::*;
use crate::{
    ambient_space::AffineSpace,
    oriented_simplex::{OrientedHyperplane, OrientedSimplex},
    simplex::Simplex,
    vector::Vector,
};
use algebraeon_rings::matrix::{Matrix, MatrixStructure};

#[derive(Debug, Clone)]
pub struct EmbeddedAffineSubspace<'f, FS: OrderedRingSignature + FieldSignature> {
    // The ordered_field of ambient_space and subspace must match
    ambient_space: AffineSpace<'f, FS>,
    embedded_space: AffineSpace<'f, FS>,
    /*
    these vectors must be equal in length to the affine dimension of subspace
    they define the embedding of subspace into ambient space
        if they are empty then they define the empty embedding
        if there is one vector then it defines the location of the embedded point
        if there are vectors [v0, v1, v2, ..., vn] then the embedding sends (a1, a2, ..., an) in subspace to v0 + a1*v1, v0 + a2*v2, ..., v0 + an*vn in ambient space
    */
    embedding_points: Vec<Vector<'f, FS>>,
}

impl<'f, FS: OrderedRingSignature + FieldSignature> EmbeddedAffineSubspace<'f, FS> {
    pub fn new_affine_span(
        ambient_space: AffineSpace<'f, FS>,
        points: Vec<Vector<'f, FS>>,
    ) -> Result<(Self, Vec<Vector<'f, FS>>), &'static str> {
        for point in &points {
            debug_assert_eq!(point.ambient_space(), ambient_space);
        }
        if !ambient_space
            .borrow()
            .are_points_affine_independent(points.iter().collect())
        {
            return Err("Affine embedding points must be affine independent");
        }
        let field = ambient_space.borrow().field();
        let embedded_space: AffineSpace<'f, FS> = AffineSpace::new_affine(field, points.len());
        let n = points.len();
        let embedded_pts = (0..n)
            .map(|i| {
                Vector::construct(embedded_space, |j| {
                    if i == 0 {
                        field.borrow().zero()
                    } else if i == j + 1 {
                        field.borrow().one()
                    } else {
                        field.borrow().zero()
                    }
                })
            })
            .collect();
        Ok((
            Self {
                ambient_space,
                embedded_space,
                embedding_points: points,
            },
            embedded_pts,
        ))
    }

    pub fn new_empty(ambient_space: AffineSpace<'f, FS>) -> Self {
        Self::new_affine_span(ambient_space, vec![]).unwrap().0
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> EmbeddedAffineSubspace<'f, FS> {
    pub fn new(
        ambient_space: AffineSpace<'f, FS>,
        root: Vector<'f, FS>,
        span: Vec<Vector<'f, FS>>,
    ) -> Result<(Self, Vec<Vector<'f, FS>>), &'static str> {
        let mut points = vec![root.clone()];
        points.extend(span.iter().map(|vec| &root + vec));
        Self::new_affine_span(ambient_space, points)
    }

    pub fn new_affine_span_linearly_dependent(
        ambient_space: AffineSpace<'f, FS>,
        points: Vec<&Vector<'f, FS>>,
    ) -> Self {
        if points.is_empty() {
            Self::new_empty(ambient_space)
        } else {
            let dim = ambient_space.borrow().linear_dimension().unwrap();
            let field = ambient_space.borrow().field();
            let mut points = points.into_iter();
            let root = points.next().unwrap();
            let span = points.map(|pt| pt - root).collect::<Vec<_>>();
            //matrix whose columns are pt - root for every other pt in points
            let mat = Matrix::construct(dim, span.len(), |r, c| span[c].coordinate(r).clone());
            let (_, _, _, pivs) = MatrixStructure::new(field.clone()).row_hermite_algorithm(mat);
            Self::new(
                ambient_space,
                root.clone(),
                pivs.into_iter().map(|i| span[i].clone()).collect(),
            )
            .unwrap()
            .0
        }
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> EmbeddedAffineSubspace<'f, FS> {
    pub fn field(&self) -> &'f FS {
        self.ambient_space.borrow().field()
    }

    pub fn ambient_space(&self) -> AffineSpace<'f, FS> {
        self.ambient_space
    }

    pub fn embedded_space(&self) -> AffineSpace<'f, FS> {
        self.embedded_space
    }

    //Let A be the affine subspace and let S be its ambient space
    //Find an affine subspace B of S obtained by linearly extending by pt
    //Return the embeddings (f, g, pt) where
    //  f : A -> B
    //  g : B -> S
    //  pt is pt in B
    pub(crate) fn extend_dimension_by_point_unsafe(
        &self,
        pt: Vector<'f, FS>,
    ) -> (
        EmbeddedAffineSubspace<'f, FS>,
        EmbeddedAffineSubspace<'f, FS>,
        Vector<'f, FS>,
    ) {
        debug_assert_eq!(self.ambient_space, pt.ambient_space());
        debug_assert!(self.unembed_point(&pt).is_none());
        let field = self.field();

        let n = self.embedded_space.borrow().affine_dimension();
        let extended_embedded_space = AffineSpace::new_affine(field, n + 1);

        (
            EmbeddedAffineSubspace {
                ambient_space: extended_embedded_space,
                embedded_space: self.embedded_space,
                // 0, e_1, e_2, ..., e_(n-1)
                embedding_points: {
                    (0..n)
                        .map(|k| {
                            Vector::construct(extended_embedded_space, |i| {
                                if k == 0 {
                                    field.zero()
                                } else {
                                    let j = k - 1;
                                    if i == j { field.one() } else { field.zero() }
                                }
                            })
                        })
                        .collect()
                },
            },
            EmbeddedAffineSubspace {
                ambient_space: self.ambient_space,
                embedded_space: extended_embedded_space,
                embedding_points: {
                    let mut pts = self.embedding_points.clone();
                    pts.push(pt);
                    pts
                },
            },
            Vector::construct(extended_embedded_space, |i| {
                if i + 1 == n {
                    self.field().one()
                } else {
                    self.field().zero()
                }
            }),
        )
    }

    pub fn get_root_and_span(&self) -> Option<(Vector<'f, FS>, Vec<Vector<'f, FS>>)> {
        let mut points = self.embedding_points.iter();
        let root = points.next()?;
        let span = points.map(|pt| pt - root).collect::<Vec<_>>();
        Some((root.clone(), span))
    }

    pub fn get_embedding_points(&self) -> &Vec<Vector<'f, FS>> {
        &self.embedding_points
    }

    pub fn embed_point(&self, pt: &Vector<'f, FS>) -> Vector<'f, FS> {
        assert_eq!(pt.ambient_space(), self.embedded_space);
        let (root, span) = self.get_root_and_span().unwrap(); //pt exists in the embedded space, so the embedded space is non-empty, so has a root and span
        let mut total = root.clone();
        for (i, vec) in span.iter().enumerate() {
            total += &vec.scalar_mul(pt.coordinate(i));
        }
        total
    }

    pub fn embed_simplex(&self, spx: &Simplex<'f, FS>) -> Simplex<'f, FS> {
        self.ambient_space()
            .simplex(spx.points().iter().map(|p| self.embed_point(p)).collect())
            .unwrap()
    }

    pub fn unembed_point(&self, pt: &Vector<'f, FS>) -> Option<Vector<'f, FS>> {
        assert_eq!(pt.ambient_space(), self.ambient_space);
        match self.get_root_and_span() {
            Some((root, span)) => {
                //solve root + x * basis = v for x
                let y = (pt - &root).into_coordinates();
                let basis_matrix = self
                    .ambient_space
                    .borrow()
                    .cols_from_vectors(span.iter().collect());
                let x = MatrixStructure::new(self.ambient_space.borrow().field().clone())
                    .col_solve(basis_matrix, &y);
                Some(self.embedded_space().vector(x?))
            }
            None => None,
        }
    }

    pub fn unembed_simplex(&self, spx: &Simplex<'f, FS>) -> Option<Simplex<'f, FS>> {
        let mut pts = vec![];
        for embedded_pt in spx.points() {
            match self.unembed_point(embedded_pt) {
                Some(pt) => {
                    pts.push(pt);
                }
                None => {
                    return None;
                }
            }
        }
        Some(self.embedded_space().simplex(pts).unwrap())
    }

    pub fn as_hyperplane_intersection(&self) -> Option<Vec<OrientedHyperplane<'f, FS>>> {
        let ambient_space = self.ambient_space();
        let field = ambient_space.field();
        match self.get_root_and_span() {
            Some((root, span)) => {
                let dim_amb = ambient_space.linear_dimension().unwrap();
                let dim_ss = span.len();
                // step 1: extend span to a basis by appending a subset of the elementary basis vectors
                //first columns = span, remaining columns = identity matrix
                let mat = Matrix::construct(dim_amb, dim_ss + dim_amb, |r, c| {
                    if c < dim_ss {
                        span[c].coordinate(r).clone()
                    } else {
                        let c = c - dim_ss;
                        if r == c { field.one() } else { field.zero() }
                    }
                });
                let (_, _, _, pivs) =
                    MatrixStructure::new(field.clone()).row_hermite_algorithm(mat);
                debug_assert_eq!(pivs.len(), dim_amb);
                #[allow(clippy::needless_range_loop)]
                for i in 0..dim_ss {
                    debug_assert_eq!(pivs[i], i); //span is linearly independent so we expect this
                }
                let extension_elementary_basis_vectors = (dim_ss..dim_amb)
                    .map(|i| pivs[i] - dim_ss)
                    .collect::<Vec<_>>();

                // step 2: take the hyperplanes formed by removing exactly one of each of the added elementary basis vectors at a time
                let hyperplanes = (0..extension_elementary_basis_vectors.len())
                    .map(|i| {
                        let ref_point = {
                            //root + e_k
                            let k = extension_elementary_basis_vectors[i];
                            Vector::construct(ambient_space, |l| {
                                field.add(
                                    root.coordinate(l),
                                    &if l == k { field.one() } else { field.zero() },
                                )
                            })
                        };
                        OrientedSimplex::new_with_positive_point(
                            ambient_space,
                            {
                                let mut points = vec![root.clone()];
                                for s in &span {
                                    points.push(&root + s);
                                }
                                for (j, &k) in extension_elementary_basis_vectors.iter().enumerate()
                                {
                                    if i != j {
                                        //push root + e_k
                                        points.push(Vector::construct(ambient_space, |l| {
                                            field.add(root.coordinate(l), &{
                                                if l == k { field.one() } else { field.zero() }
                                            })
                                        }));
                                    }
                                }
                                points
                            },
                            &ref_point,
                        )
                        .unwrap()
                        .into_oriented_hyperplane()
                    })
                    .collect::<Vec<_>>();
                debug_assert_eq!(hyperplanes.len(), dim_amb - dim_ss);
                Some(hyperplanes)
            }
            None => None,
        }
    }
}

pub fn compose_affine_embeddings<'f, FS: OrderedRingSignature + FieldSignature>(
    _a_to_b: EmbeddedAffineSubspace<'f, FS>,
    _b_to_c: EmbeddedAffineSubspace<'f, FS>,
) -> EmbeddedAffineSubspace<'f, FS> {
    todo!() // call b_to_c.embed on the defining points of a_to_b
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Rational;

    use super::*;

    #[test]
    fn make_affine_subspace() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 3);
        let v1 = space.vector([1, 1, 1]);
        let v2 = space.vector([1, 0, 0]);
        let v3 = space.vector([0, 1, 0]);
        let s = EmbeddedAffineSubspace::new(space, v1, vec![v2, v3]);
        s.unwrap();

        let space = AffineSpace::new_linear(Rational::structure_ref(), 3);
        let v1 = space.vector([1, 1, 1]);
        let v2 = space.vector([1, 2, 0]);
        let v3 = space.vector([-2, -4, 0]);
        let s = EmbeddedAffineSubspace::new(space, v1, vec![v2, v3]);
        assert!(s.is_err());
    }

    #[test]
    fn affine_subspace_embed_and_unembed() {
        //1d embedded in 2d
        {
            let plane = AffineSpace::new_linear(Rational::structure_ref(), 2);
            //the line x + y = 2
            let (line, _) = EmbeddedAffineSubspace::new(
                plane,
                plane.vector([1, 1]),
                vec![plane.vector([1, -1])],
            )
            .unwrap();

            assert_eq!(
                line.embed_point(&line.embedded_space().vector([-3])),
                plane.vector([-2, 4])
            );

            assert_eq!(
                line.unembed_point(&plane.vector([-1, 3])),
                Some(line.embedded_space().vector([-2]))
            );

            assert_eq!(line.unembed_point(&plane.vector([1, 2])), None);
        }

        //2d embedded in 3d
        {
            let space = AffineSpace::new_linear(Rational::structure_ref(), 3);
            let (plane, _) = EmbeddedAffineSubspace::new(
                space,
                space.vector([3, 1, 2]),
                vec![space.vector([4, 2, 1]), space.vector([1, -1, 2])],
            )
            .unwrap();

            assert_eq!(
                plane.embed_point(&plane.embedded_space().vector([-3, 2])),
                space.vector([-7, -7, 3])
            );

            assert_eq!(
                plane.unembed_point(&space.vector([0, -2, 3],)),
                Some(plane.embedded_space().vector([-1, 1]))
            );

            assert_eq!(plane.unembed_point(&space.vector([1, 2, 1],)), None);
        }
    }

    #[test]
    fn extend_by_point_embedding_composition() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 4);
        let v1 = space.vector([1, 2, 1, 1]);
        let v2 = space.vector([1, -2, 2, 0]);
        let v3 = space.vector([2, 1, 0, 2]);
        let (h, _) = EmbeddedAffineSubspace::new(space, v1, vec![v2, v3]).unwrap();
        let v4 = space.vector([0, 3, -2, 1]);
        let (f, g, v4_inv) = h.extend_dimension_by_point_unsafe(v4.clone());
        assert_eq!(g.embed_point(&v4_inv), v4);

        let x = h
            .embedded_space()
            .vector([Rational::from(5), Rational::from(7)]);
        //check that g(f(x)) = h(x)
        assert_eq!(g.embed_point(&f.embed_point(&x)), h.embed_point(&x));
    }
}
