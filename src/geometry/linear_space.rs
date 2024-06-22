use std::sync::atomic::AtomicUsize;

use crate::rings::linear::matrix::{Matrix, MatrixStructure};

use super::*;
#[derive(Debug, Clone)]
pub struct AffineSpace<FS: OrderedRingStructure + FieldStructure> {
    ordered_field: Rc<FS>,
    //affine dimension = dimension + 1. affine dimension = 0 if linear dimension = None
    //dimension = None -> empty set
    dimension: Option<usize>,
    ident: usize,
}

impl<FS: OrderedRingStructure + FieldStructure> PartialEq for AffineSpace<FS> {
    fn eq(&self, other: &Self) -> bool {
        #[cfg(debug_assertions)]
        if self.ident == other.ident {
            assert_eq!(self.dimension, other.dimension);
            assert_eq!(self.ordered_field, other.ordered_field);
        }
        self.ident == other.ident
    }
}

impl<FS: OrderedRingStructure + FieldStructure> Eq for AffineSpace<FS> {}

impl<FS: OrderedRingStructure + FieldStructure + Hash> Hash for AffineSpace<FS> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ordered_field.hash(state);
        self.dimension.hash(state);
    }
}

impl<FS: OrderedRingStructure + FieldStructure> AffineSpace<FS> {
    fn new_impl(ordered_field: Rc<FS>, dimension: Option<usize>) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        Self {
            ordered_field,
            dimension,
            ident: COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        }
    }

    pub fn new_empty(ordered_field: Rc<FS>) -> Self {
        Self::new_impl(ordered_field, None)
    }

    pub fn new(ordered_field: Rc<FS>, dimension: usize) -> Self {
        Self::new_impl(ordered_field, Some(dimension))
    }


    // //this is supposed to be an affine space, so no origin
    // pub fn origin<'a>(&'a self) -> Option<Vector<FS, &Self>> {
    //     Some(Vector::new(
    //         self,
    //         (0..self.dimension?)
    //             .map(|i| self.ordered_field.zero())
    //             .collect(),
    //     ))
    // }

    pub fn ordered_field(&self) -> Rc<FS> {
        self.ordered_field.clone()
    }

    pub fn linear_dimension(&self) -> Option<usize> {
        self.dimension
    }

    pub fn affine_dimension(&self) -> usize {
        match self.dimension {
            Some(d) => d + 1,
            None => 0,
        }
    }

    pub fn rows_from_vectors(&self, vecs: Vec<&Vector<FS, impl Borrow<Self>>>) -> Matrix<FS::Set> {
        for vec in &vecs {
            assert_eq!(self, vec.ambient_space().borrow());
        }
        Matrix::construct(vecs.len(), self.linear_dimension().unwrap(), |r, c| {
            vecs[r].coordinate(c).clone()
        })
    }

    pub fn cols_from_vectors(&self, vecs: Vec<&Vector<FS, impl Borrow<Self>>>) -> Matrix<FS::Set> {
        self.rows_from_vectors(vecs).transpose()
    }

    pub fn determinant(&self, vecs: Vec<&Vector<FS, impl Borrow<Self>>>) -> FS::Set {
        MatrixStructure::new(self.ordered_field())
            .det(self.rows_from_vectors(vecs))
            .unwrap()
    }

    pub fn rank(&self, vecs: Vec<&Vector<FS, impl Borrow<Self>>>) -> usize {
        MatrixStructure::new(self.ordered_field()).rank(self.rows_from_vectors(vecs))
    }

    pub fn are_points_nondegenerage(
        &self,
        points: Vec<&Vector<FS, impl Borrow<Self> + Clone>>,
    ) -> bool {
        for point in &points {
            assert_eq!(self, point.ambient_space().borrow());
        }
        if points.len() != 0 {
            let vecs = (1..points.len())
                .map(|i| points[i] - points[0])
                .collect::<Vec<_>>();
            let mat = self.rows_from_vectors(vecs.iter().collect());
            MatrixStructure::new(self.ordered_field()).rank(mat) == vecs.len()
        } else {
            true
        }
    }
}

pub fn vectors_from_rows<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<FS, SP>> {
    assert_eq!(mat.cols(), sp.borrow().linear_dimension().unwrap());
    (0..mat.rows())
        .map(|r| {
            Vector::new(
                sp.clone(),
                (0..mat.cols())
                    .map(|c| mat.at(r, c).unwrap().clone())
                    .collect(),
            )
        })
        .collect()
}

pub fn vectors_from_cols<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<FS, SP>> {
    assert_eq!(mat.rows(), sp.borrow().linear_dimension().unwrap());
    vectors_from_rows(sp, &mat.transpose_ref())
}

pub fn vector_from_row<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vector<FS, SP> {
    assert_eq!(mat.rows(), 1);
    assert_eq!(mat.cols(), sp.borrow().linear_dimension().unwrap());
    vectors_from_rows(sp, mat).pop().unwrap()
}

pub fn vector_from_col<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vector<FS, SP> {
    assert_eq!(mat.rows(), sp.borrow().linear_dimension().unwrap());
    assert_eq!(mat.cols(), 1);
    vector_from_row(sp, &mat.transpose_ref())
}

pub fn common_space<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
>(
    space1: SP,
    space2: SP,
) -> Option<SP> {
    if space1.borrow() == space2.borrow() {
        Some(space1)
    } else {
        None
    }
}
