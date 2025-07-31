use super::*;
use algebraeon_rings::matrix::{Matrix, MatrixStructure};
use std::sync::atomic::AtomicUsize;

/// An affine space over a field.
/// affine_dimension = 0 => the empty space
/// affine_dimension = 1 => one point space
/// affine_dimension = 2 => a line
/// affine_dimension = 3 => a plane
/// ...
#[derive(Debug, Clone)]
pub struct AffineSpace<'f, FS: FieldSignature> {
    field: &'f FS,
    // linear dimension = affine dimension - 1
    affine_dimension: usize,
    ident: usize,
}

impl<'f, FS: FieldSignature> PartialEq for AffineSpace<'f, FS> {
    fn eq(&self, other: &Self) -> bool {
        #[cfg(debug_assertions)]
        if self.ident == other.ident {
            assert_eq!(self.affine_dimension, other.affine_dimension);
            assert_eq!(self.field(), other.field());
        }
        self.ident == other.ident
    }
}

impl<'f, FS: FieldSignature> Eq for AffineSpace<'f, FS> {}

impl<'f, FS: FieldSignature + Hash> Hash for AffineSpace<'f, FS> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ident.hash(state);
    }
}

impl<'f, FS: FieldSignature> AffineSpace<'f, FS> {
    pub fn new_affine(field: &'f FS, affine_dimension: usize) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        Self {
            field,
            affine_dimension,
            ident: COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
        }
    }

    pub fn new_empty(field: &'f FS) -> Self {
        Self::new_affine(field, 0)
    }

    pub fn new_linear(field: &'f FS, linear_dimension: usize) -> Self {
        Self::new_affine(field, linear_dimension + 1)
    }

    pub fn field(&self) -> &'f FS {
        self.field
    }

    pub fn origin<SP: Borrow<Self> + Clone + From<Self>>(&self) -> Option<Vector<'f, FS, SP>> {
        Some(Vector::construct(self.clone().into(), |_i| {
            self.field().zero()
        }))
    }

    pub fn linear_dimension(&self) -> Option<usize> {
        if self.affine_dimension == 0 {
            None
        } else {
            Some(self.affine_dimension - 1)
        }
    }

    pub fn affine_dimension(&self) -> usize {
        self.affine_dimension
    }

    #[allow(clippy::needless_pass_by_value)]
    pub fn rows_from_vectors(
        &self,
        vecs: Vec<&Vector<'f, FS, impl Borrow<Self>>>,
    ) -> Matrix<FS::Set> {
        for vec in &vecs {
            assert_eq!(self, vec.ambient_space().borrow());
        }
        Matrix::construct(vecs.len(), self.linear_dimension().unwrap(), |r, c| {
            vecs[r].coordinate(c).clone()
        })
    }

    pub fn cols_from_vectors(
        &self,
        vecs: Vec<&Vector<'f, FS, impl Borrow<Self>>>,
    ) -> Matrix<FS::Set> {
        self.rows_from_vectors(vecs).transpose()
    }

    pub fn determinant(&self, vecs: Vec<&Vector<'f, FS, impl Borrow<Self>>>) -> FS::Set {
        MatrixStructure::new(self.field().clone())
            .det(self.rows_from_vectors(vecs))
            .unwrap()
    }

    pub fn rank(&self, vecs: Vec<&Vector<'f, FS, impl Borrow<Self>>>) -> usize {
        MatrixStructure::new(self.field().clone()).rank(self.rows_from_vectors(vecs))
    }

    #[allow(clippy::needless_pass_by_value)]
    pub fn are_points_affine_independent(
        &self,
        points: Vec<&Vector<'f, FS, impl Borrow<Self> + Clone>>,
    ) -> bool {
        for point in &points {
            assert_eq!(self, point.ambient_space().borrow());
        }
        if points.is_empty() {
            true
        } else {
            let vecs = (1..points.len())
                .map(|i| points[i] - points[0])
                .collect::<Vec<_>>();
            let mat = self.rows_from_vectors(vecs.iter().collect());
            // println!("{:?}", mat);
            // println!("{:?} {:?}", vecs.len(), MatrixStructure::new(self.field()).rank(mat.clone()));
            MatrixStructure::new(self.field().clone()).rank(mat) == vecs.len()
        }
    }
}

pub fn vectors_from_rows<'f, FS: FieldSignature + 'f, SP: Borrow<AffineSpace<'f, FS>> + Clone>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<'f, FS, SP>> {
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

pub fn vectors_from_cols<'f, FS: FieldSignature + 'f, SP: Borrow<AffineSpace<'f, FS>> + Clone>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<'f, FS, SP>> {
    assert_eq!(mat.rows(), sp.borrow().linear_dimension().unwrap());
    vectors_from_rows(sp, &mat.transpose_ref())
}

pub fn vector_from_row<'f, FS: FieldSignature + 'f, SP: Borrow<AffineSpace<'f, FS>> + Clone>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vector<'f, FS, SP> {
    assert_eq!(mat.rows(), 1);
    assert_eq!(mat.cols(), sp.borrow().linear_dimension().unwrap());
    vectors_from_rows(sp, mat).pop().unwrap()
}

pub fn vector_from_col<'f, FS: FieldSignature + 'f, SP: Borrow<AffineSpace<'f, FS>> + Clone>(
    sp: SP,
    mat: &Matrix<FS::Set>,
) -> Vector<'f, FS, SP> {
    assert_eq!(mat.rows(), sp.borrow().linear_dimension().unwrap());
    assert_eq!(mat.cols(), 1);
    vector_from_row(sp, &mat.transpose_ref())
}

pub fn common_space<'f, FS: FieldSignature + 'f, SP: Borrow<AffineSpace<'f, FS>> + Clone>(
    space1: SP,
    space2: SP,
) -> Option<SP> {
    if space1.borrow() == space2.borrow() {
        Some(space1)
    } else {
        None
    }
}
