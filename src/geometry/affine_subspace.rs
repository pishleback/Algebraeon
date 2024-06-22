use crate::rings::linear::matrix::MatrixStructure;

use super::*;

#[derive(Debug, Clone)]
enum AffineSubspaceEmbedding<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    ESP: Borrow<AffineSpace<FS>> + Clone,
> {
    Empty {
        subspace: ESP, //must be empty i.e. affine dimension = 0 i.e. linear dimension = None
    },
    NonEmpty {
        subspace: ESP, //dimension equal to the length of the basis
        root: Vector<FS, SP>,
        basis: Vec<Vector<FS, SP>>,
    },
}

#[derive(Debug, Clone)]
pub struct EmbeddedAffineSubspace<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    ESP: Borrow<AffineSpace<FS>> + Clone,
> {
    //The ordered_field of ambient_space and subspace must match
    ambient_space: SP,
    embedding: AffineSubspaceEmbedding<FS, SP, ESP>,
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        ESP: Borrow<AffineSpace<FS>> + From<AffineSpace<FS>> + Clone,
    > EmbeddedAffineSubspace<FS, SP, ESP>
{
    pub fn new_empty(ambient_space: SP) -> Self {
        let ordered_field = ambient_space.borrow().ordered_field();
        Self {
            ambient_space,
            embedding: AffineSubspaceEmbedding::Empty {
                subspace: AffineSpace::new_empty(ordered_field.clone()).into(),
            },
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    EmbeddedAffineSubspace<FS, SP, AffineSpace<FS>>
{
    pub fn new(
        ambient_space: SP,
        root: Vector<FS, SP>,
        span: Vec<Vector<FS, SP>>,
    ) -> Result<Self, &'static str> {
        assert_eq!(root.ambient_space().borrow(), ambient_space.borrow());
        for vec in &span {
            assert_eq!(vec.ambient_space().borrow(), ambient_space.borrow());
        }
        let ordered_field = ambient_space.borrow().ordered_field();
        if ambient_space.borrow().rank(span.iter().collect()) != span.len() {
            Err("The span must be linearly independent when constructing an affine subspace")
        } else {
            let d = ambient_space.borrow().linear_dimension().unwrap();
            let embedding = AffineSubspaceEmbedding::NonEmpty {
                subspace: AffineSpace::new(ordered_field.clone(), span.len()),
                root: root,
                basis: span,
            };
            Ok(Self {
                ambient_space,
                embedding,
            })
        }
    }
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        ESP: Borrow<AffineSpace<FS>> + Clone,
    > EmbeddedAffineSubspace<FS, SP, ESP>
{
    pub fn ordered_field(&self) -> Rc<FS> {
        self.ambient_space.borrow().ordered_field()
    }

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }

    // pub fn embedded_dimension(&self) -> Option<usize> {
    //     match &self.embedding {
    //         AffineSubspaceEmbedding::Empty => None,
    //         AffineSubspaceEmbedding::NonEmpty {
    //             subspace,
    //             root,
    //             basis,
    //         } => {
    //             debug_assert_eq!(subspace.borrow().dimension(), basis.len());
    //             Some(basis.len())
    //         }
    //     }
    // }

    pub fn embedded_space(&self) -> ESP {
        match &self.embedding {
            AffineSubspaceEmbedding::Empty { subspace } => {
                debug_assert_eq!(subspace.borrow().affine_dimension(), 0);
                subspace.clone()
            }
            AffineSubspaceEmbedding::NonEmpty { subspace, .. } => subspace.clone(),
        }
    }

    //Let A be the affine subspace and let S be its ambient space
    //Find an affine subspace B of S obtained by linearly extending by pt
    //Return the embeddings (f, g) where f : A -> B and g : B -> S
    pub fn extend_dimension_by_point_unsafe(
        &self,
        pt: &Vector<FS, SP>,
    ) -> (
        EmbeddedAffineSubspace<FS, Rc<AffineSpace<FS>>, ESP>,
        EmbeddedAffineSubspace<FS, SP, Rc<AffineSpace<FS>>>,
    ) {
        debug_assert_eq!(self.ambient_space.borrow(), pt.ambient_space().borrow());
        debug_assert!(self.unembed_point(pt).is_none());
        let ordered_field = self.ordered_field();
        match &self.embedding {
            AffineSubspaceEmbedding::Empty { subspace } => {
                let extended_linear_space = Rc::new(AffineSpace::new(ordered_field, 0));
                (
                    EmbeddedAffineSubspace {
                        ambient_space: extended_linear_space.clone(),
                        embedding: AffineSubspaceEmbedding::Empty {
                            subspace: subspace.clone(),
                        },
                    },
                    EmbeddedAffineSubspace {
                        ambient_space: self.ambient_space.clone(),
                        embedding: AffineSubspaceEmbedding::NonEmpty {
                            subspace: extended_linear_space,
                            root: pt.clone(),
                            basis: vec![],
                        },
                    },
                )
            }
            AffineSubspaceEmbedding::NonEmpty {
                subspace,
                root,
                basis,
            } => {
                let n = subspace.borrow().linear_dimension().unwrap();
                let extended_linear_space = Rc::new(AffineSpace::new(ordered_field.clone(), n + 1));
                (
                    EmbeddedAffineSubspace {
                        ambient_space: extended_linear_space.clone(),
                        embedding: AffineSubspaceEmbedding::NonEmpty {
                            subspace: subspace.clone(),
                            root: Vector::zero(extended_linear_space.clone()),
                            basis: (0..n)
                                .map(|k| {
                                    Vector::construct(extended_linear_space.clone(), |i| {
                                        match i == k {
                                            true => ordered_field.one(),
                                            false => ordered_field.zero(),
                                        }
                                    })
                                })
                                .collect(),
                        },
                    },
                    EmbeddedAffineSubspace {
                        ambient_space: self.ambient_space.clone(),
                        embedding: AffineSubspaceEmbedding::NonEmpty {
                            subspace: extended_linear_space.clone(),
                            root: root.clone(),
                            basis: {
                                let mut basis = basis.clone();
                                basis.push(pt - root);
                                basis
                            },
                        },
                    },
                )
            }
        }
    }

    pub fn embed_point(&self, p: &Vector<FS, ESP>) -> Vector<FS, SP> {
        match &self.embedding {
            AffineSubspaceEmbedding::Empty { .. } => panic!(),
            AffineSubspaceEmbedding::NonEmpty {
                subspace,
                root,
                basis,
            } => {
                assert_eq!(p.ambient_space().borrow(), subspace.borrow());
                let mut total = root.clone();
                for (i, b) in basis.iter().enumerate() {
                    total += &b.scalar_mul(p.coordinate(i));
                }
                total
            }
        }
    }

    pub fn unembed_point(&self, p: &Vector<FS, SP>) -> Option<Vector<FS, ESP>> {
        assert_eq!(p.ambient_space().borrow(), self.ambient_space.borrow());
        //solve root + x * basis = v for x
        match &self.embedding {
            AffineSubspaceEmbedding::Empty { .. } => None,
            AffineSubspaceEmbedding::NonEmpty {
                subspace,
                root,
                basis,
            } => {
                let y = (p - root).into_col();
                let basis_matrix = self
                    .ambient_space
                    .borrow()
                    .cols_from_vectors(basis.iter().collect());
                let x = MatrixStructure::new(self.ambient_space.borrow().ordered_field())
                    .col_solve(&basis_matrix, y);
                Some(vector_from_col(self.embedded_space(), &x?))
            }
        }
    }

    pub fn embed_vector(&self, v: &Vector<FS, ESP>) -> Vector<FS, SP> {
        match &self.embedding {
            AffineSubspaceEmbedding::Empty { .. } => panic!(),
            AffineSubspaceEmbedding::NonEmpty {
                subspace,
                root,
                basis,
            } => {
                assert_eq!(v.ambient_space().borrow(), subspace.borrow());
                let mut total = Vector::zero(self.ambient_space.clone());
                for (i, b) in basis.iter().enumerate() {
                    total += &b.scalar_mul(v.coordinate(i));
                }
                total
            }
        }
    }

    pub fn unembed_vector(&self, v: &Vector<FS, SP>) -> Option<Vector<FS, ESP>> {
        assert_eq!(v.ambient_space().borrow(), self.ambient_space.borrow());
        todo!()
    }
}

pub fn compose_affine_embeddings<
    FS: OrderedRingStructure + FieldStructure,
    SPA: Borrow<AffineSpace<FS>> + Clone,
    SPB: Borrow<AffineSpace<FS>> + Clone,
    SPC: Borrow<AffineSpace<FS>> + Clone,
>(
    a_to_b: EmbeddedAffineSubspace<FS, SPB, SPA>,
    b_to_c: EmbeddedAffineSubspace<FS, SPC, SPB>,
) -> EmbeddedAffineSubspace<FS, SPC, SPA> {
    todo!() // call b_to_c.embed on the defining points of a_to_b
}

#[cfg(test)]
mod tests {
    use malachite_q::Rational;

    use crate::rings::structure::StructuredType;

    use super::*;

    #[test]
    fn make_affine_subspace() {
        let space = AffineSpace::new(Rational::structure(), 3);
        let v1 = Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(1), Rational::from(1)],
        );
        let v2 = Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(0), Rational::from(0)],
        );
        let v3 = Vector::new(
            &space,
            vec![Rational::from(0), Rational::from(1), Rational::from(0)],
        );
        let s = EmbeddedAffineSubspace::new(&space, v1, vec![v2, v3]);
        assert!(s.is_ok());

        let space = AffineSpace::new(Rational::structure(), 3);
        let v1 = Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(1), Rational::from(1)],
        );
        let v2 = Vector::new(
            &space,
            vec![Rational::from(1), Rational::from(2), Rational::from(0)],
        );
        let v3 = Vector::new(
            &space,
            vec![Rational::from(-2), Rational::from(-4), Rational::from(0)],
        );
        let s = EmbeddedAffineSubspace::new(&space, v1, vec![v2, v3]);
        assert!(s.is_err());
    }

    #[test]
    fn affine_subspace_embed_and_unembed() {
        //1d embedded in 2d
        {
            let plane = AffineSpace::new(Rational::structure(), 2);
            //the line x + y = 2
            let line = EmbeddedAffineSubspace::new(
                &plane,
                Vector::new(&plane, vec![Rational::from(1), Rational::from(1)]),
                vec![Vector::new(
                    &plane,
                    vec![Rational::from(1), Rational::from(-1)],
                )],
            )
            .unwrap();

            assert_eq!(
                line.embed_point(&Vector::new(
                    line.embedded_space(),
                    vec![Rational::from(-3)],
                )),
                Vector::new(&plane, vec![Rational::from(-2), Rational::from(4)])
            );

            assert_eq!(
                line.unembed_point(&Vector::new(
                    &plane,
                    vec![Rational::from(-1), Rational::from(3)],
                )),
                Some(Vector::new(line.embedded_space(), vec![Rational::from(-2)],))
            );

            assert_eq!(
                line.unembed_point(&Vector::new(
                    &plane,
                    vec![Rational::from(1), Rational::from(2)],
                )),
                None
            );
        }

        //2d embedded in 3d
        {
            let space = AffineSpace::new(Rational::structure(), 3);
            let plane = EmbeddedAffineSubspace::new(
                &space,
                Vector::new(
                    &space,
                    vec![Rational::from(3), Rational::from(1), Rational::from(2)],
                ),
                vec![
                    Vector::new(
                        &space,
                        vec![Rational::from(4), Rational::from(2), Rational::from(1)],
                    ),
                    Vector::new(
                        &space,
                        vec![Rational::from(1), Rational::from(-1), Rational::from(2)],
                    ),
                ],
            )
            .unwrap();

            assert_eq!(
                plane.embed_point(&Vector::new(
                    plane.embedded_space(),
                    vec![Rational::from(-3), Rational::from(2)],
                )),
                Vector::new(
                    &space,
                    vec![Rational::from(-7), Rational::from(-7), Rational::from(3)]
                )
            );

            assert_eq!(
                plane.unembed_point(&Vector::new(
                    &space,
                    vec![Rational::from(0), Rational::from(-2), Rational::from(3)],
                )),
                Some(Vector::new(
                    plane.embedded_space(),
                    vec![Rational::from(-1), Rational::from(1)],
                ))
            );

            assert_eq!(
                plane.unembed_point(&Vector::new(
                    &space,
                    vec![Rational::from(1), Rational::from(2), Rational::from(2)],
                )),
                None
            );
        }
    }
}
