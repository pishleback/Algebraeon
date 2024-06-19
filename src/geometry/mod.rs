use std::borrow::Borrow;
use std::rc::Rc;

use crate::rings::{
    linear::subspace::{
        AffineLattice, AffineLatticeStructure, LinearLattice, LinearLatticeStructure,
    },
    ring_structure::structure::{FieldStructure, OrderedRingStructure},
};

mod linear_space {

    use crate::rings::linear::matrix::Matrix;

    use super::*;
    #[derive(Debug, Clone)]
    pub struct LinearSpace<FS: OrderedRingStructure + FieldStructure> {
        ordered_field: Rc<FS>,
        dimension: usize,
    }

    impl<FS: OrderedRingStructure + FieldStructure> LinearSpace<FS> {
        pub fn new(ordered_field: Rc<FS>, dimension: usize) -> Self {
            Self {
                ordered_field,
                dimension,
            }
        }

        //this is supposed to be an affine space, so no origin
        pub fn origin<'a>(&'a self) -> Vector<FS, &Self> {
            Vector::new(
                self,
                (0..self.dimension)
                    .map(|i| self.ordered_field.zero())
                    .collect(),
            )
        }

        pub fn ordered_field(&self) -> Rc<FS> {
            self.ordered_field.clone()
        }

        pub fn dimension(&self) -> usize {
            self.dimension
        }

        pub fn vectors_from_rows<'a>(&'a self, mat: &Matrix<FS::Set>) -> Vec<Vector<FS, &Self>> {
            assert_eq!(mat.cols(), self.dimension);
            (0..mat.rows())
                .map(|r| {
                    Vector::new(
                        self,
                        (0..mat.cols())
                            .map(|c| mat.at(r, c).unwrap().clone())
                            .collect(),
                    )
                })
                .collect()
        }

        pub fn vectors_from_cols<'a>(&'a self, mat: &Matrix<FS::Set>) -> Vec<Vector<FS, &Self>> {
            assert_eq!(mat.rows(), self.dimension);
            self.vectors_from_rows(&mat.transpose_ref())
        }

        pub fn vector_from_row<'a>(&'a self, mat: &Matrix<FS::Set>) -> Vector<FS, &Self> {
            assert_eq!(mat.rows(), 1);
            assert_eq!(mat.cols(), self.dimension);
            self.vectors_from_rows(mat).pop().unwrap()
        }

        pub fn vector_from_col<'a>(&'a self, mat: &Matrix<FS::Set>) -> Vector<FS, &Self> {
            assert_eq!(mat.rows(), self.dimension);
            assert_eq!(mat.cols(), 1);
            self.vector_from_row(&mat.transpose_ref())
        }
    }

    pub fn common_space<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>>(
        space1: SP,
        space2: SP,
    ) -> Option<SP> {
        if std::ptr::eq(space1.borrow(), space2.borrow()) {
            Some(space1)
        } else {
            None
        }
    }
}
pub use linear_space::*;

mod coordinates;
pub use coordinates::*;

/*
mod linear_subspace {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct LinearSubspace<FS : OrderedRingStructure + FieldStructure, SP:Borrow<LinearSpace<FS>>>  {
        //The ordered_field of abmient_space and subspace must match
        abmient_space: &'a Space<FS>,
        subspace: Space<FS>,
        embedding: LinearLattice<FS::Set>,
    }

    impl<FS : OrderedRingStructure + FieldStructure, SP:Borrow<LinearSpace<FS>>>  LinearSubspace<FS, SP> {
        // pub fn new(&self, abmient_space: &'a Space<FS>) -> Self {
        //     todo!()
        // }

        pub fn ordered_field(&self) -> Rc<FS> {
            let f = self.abmient_space.ordered_field();
            debug_assert_eq!(f, self.subspace.ordered_field());
            f
        }

        pub fn dimension(&self) -> usize {
            LinearLatticeStructure::new(self.ordered_field()).rank(&self.embedding)
        }

        pub fn space(&self) -> &Space<FS> {
            &self.subspace
        }
    }
}
pub use linear_subspace::*;
*/

pub mod hyperplanes {
    use super::*;

    mod affine_subspace {
        use super::*;

        #[derive(Debug, Clone)]
        pub struct AffineSubspace<
            FS: OrderedRingStructure + FieldStructure,
            SP: Borrow<LinearSpace<FS>>,
        > {
            //The ordered_field of abmient_space and subspace must match
            abmient_space: SP,
            subspace: LinearSpace<FS>,
            embedding: AffineLattice<FS::Set>,
        }

        impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>>
            AffineSubspace<FS, SP>
        {
            pub fn ordered_field(&self) -> Rc<FS> {
                let f = self.abmient_space.borrow().ordered_field();
                debug_assert_eq!(f, self.subspace.ordered_field());
                f
            }

            pub fn dimension(&self) -> Option<usize> {
                AffineLatticeStructure::new(self.ordered_field()).rank(&self.embedding)
            }

            pub fn space(&self) -> &LinearSpace<FS> {
                &self.subspace
            }

            pub fn embed(&self, v: &Vector<FS, SP>) -> Vector<FS, SP> {
                todo!()
            }

            pub fn unembed(&self, v: &Vector<FS, SP>) -> Option<Vector<FS, SP>> {
                todo!()
            }
        }
    }
    pub use affine_subspace::*;

    mod oriented_hyperplane {
        #[derive(Debug, Clone)]
        pub struct OrientedHyperplane {}
    }
    pub use oriented_hyperplane::*;

    mod hyperplane {
        #[derive(Debug, Clone)]
        pub struct Hyperplane {}
    }
    pub use hyperplane::*;
}

pub mod triangles {
    use super::*;

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> PartialOrd
        for Vector<FS, SP>
    {
        fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
            let space = common_space(
                self.abmient_space().borrow(),
                other.abmient_space().borrow(),
            )?;
            for i in 0..space.dimension() {
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

    impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> Ord
        for Vector<FS, SP>
    {
        fn cmp(&self, other: &Self) -> std::cmp::Ordering {
            match self.partial_cmp(other) {
                Some(ans) => ans,
                None => panic!(),
            }
        }
    }

    mod oriented_simplex {
        use super::*;

        #[derive(Debug, Clone)]
        pub struct OrientedSimplex {}
    }
    pub use oriented_simplex::*;

    mod simplex {
        use super::*;

        #[derive(Debug, Clone)]
        pub struct Simplex {}
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
        use super::*;

        #[derive(Debug, Clone)]
        pub struct ConvexHull {}
    }
    pub use convex_hull::*;
}

#[cfg(test)]
mod tests {
    use crate::rings::{linear::matrix::Matrix, structure::StructuredType};
    use malachite_q::Rational;

    use super::*;

    #[test]
    fn test_vector_from_mat() {
        let space = LinearSpace::new(Rational::structure(), 2);
        let mat = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(1), Rational::from(2)],
            vec![Rational::from(3), Rational::from(4)],
        ]);

        mat.pprint();

        let mut vecs = space.vectors_from_rows(&mat);
        let v2 = vecs.pop().unwrap();
        let v1 = vecs.pop().unwrap();
        println!("v1 = {:?}", v1);
        println!("v2 = {:?}", v2);

        assert_eq!(
            v1,
            Vector::new(&space, vec![Rational::from(1), Rational::from(2)])
        );
        assert_eq!(
            v2,
            Vector::new(&space, vec![Rational::from(3), Rational::from(4)])
        );
    }
}
