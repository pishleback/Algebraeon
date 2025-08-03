use crate::ambient_space::{AffineSpace, common_space};

use super::*;
use algebraeon_rings::matrix::Matrix;
use std::borrow::Borrow;
use std::hash::Hash;

#[derive(Clone)]
pub struct Vector<'f, FS: FieldSignature + 'f> {
    ambient_space: AffineSpace<'f, FS>,
    coordinates: Vec<FS::Set>, //length equal to ambient_space.dimension()
}

impl<'f, FS: FieldSignature> std::fmt::Debug for Vector<'f, FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Vector")
            .field("coordinates", &self.coordinates)
            .finish()
    }
}

impl<'f, FS: FieldSignature> PartialEq for Vector<'f, FS> {
    fn eq(&self, other: &Self) -> bool {
        match common_space(self.ambient_space, other.ambient_space) {
            Some(space) => {
                let n = space.linear_dimension().unwrap();
                (0..n).all(|i| space.field().equal(self.coordinate(i), other.coordinate(i)))
            }
            None => false,
        }
    }
}

impl<'f, FS: FieldSignature> Eq for Vector<'f, FS> {}

impl<'f, FS: FieldSignature> Hash for Vector<'f, FS>
where
    FS::Set: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // self.ambient_space.borrow().hash(state);
        self.coordinates.hash(state);
    }
}
impl<'f, FS: FieldSignature> Vector<'f, FS> {
    pub fn ambient_space(&self) -> AffineSpace<'f, FS> {
        self.ambient_space
    }

    pub(crate) fn new(
        ambient_space: AffineSpace<'f, FS>,
        coordinates: impl IntoIterator<Item = impl Into<FS::Set>>,
    ) -> Self {
        let coordinates = coordinates
            .into_iter()
            .map(|c| c.into())
            .collect::<Vec<_>>();
        assert_eq!(ambient_space.linear_dimension().unwrap(), coordinates.len());
        Self {
            ambient_space,
            coordinates,
        }
    }

    pub fn construct(
        ambient_space: AffineSpace<'f, FS>,
        coordinate_func: impl FnMut(usize) -> FS::Set,
    ) -> Self {
        let coordinates = (0..ambient_space.borrow().linear_dimension().unwrap())
            .map(coordinate_func)
            .collect();
        Self {
            ambient_space,
            coordinates,
        }
    }

    pub fn zero(ambient_space: AffineSpace<'f, FS>) -> Self {
        let field = ambient_space.borrow().field().clone();
        Self::construct(ambient_space, |_i| field.zero())
    }

    pub fn coordinate(&self, i: usize) -> &FS::Set {
        self.coordinates.get(i).unwrap()
    }

    pub fn coordinate_mut(&mut self, i: usize) -> &mut FS::Set {
        self.coordinates.get_mut(i).unwrap()
    }

    pub fn into_coordinates(self) -> Vec<FS::Set> {
        self.coordinates
    }

    pub fn into_row(&self) -> Matrix<FS::Set> {
        Matrix::construct(
            1,
            self.ambient_space().linear_dimension().unwrap(),
            |_r, c| self.coordinate(c).clone(),
        )
    }

    pub fn into_col(&self) -> Matrix<FS::Set> {
        Matrix::construct(
            self.ambient_space().linear_dimension().unwrap(),
            1,
            |r, _c| self.coordinate(r).clone(),
        )
    }
}

// -&vector
impl<'f, FS: FieldSignature> std::ops::Neg for &Vector<'f, FS> {
    type Output = Vector<'f, FS>;

    fn neg(self) -> Self::Output {
        Vector {
            ambient_space: self.ambient_space,
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ambient_space().field().neg(x))
                .collect(),
        }
    }
}

// &vector + &vector
impl<'f, FS: FieldSignature> std::ops::Add<&Vector<'f, FS>> for &Vector<'f, FS> {
    type Output = Vector<'f, FS>;

    fn add(self, other: &Vector<'f, FS>) -> Self::Output {
        match common_space(self.ambient_space, other.ambient_space) {
            Some(space) => {
                let n = space.linear_dimension().unwrap();
                let coordinates = (0..n)
                    .map(|i| space.field().add(self.coordinate(i), other.coordinate(i)))
                    .collect();
                Vector {
                    ambient_space: space,
                    coordinates,
                }
            }
            None => panic!("Can't add vectors belonging to different spaces"),
        }
    }
}

// mut vector += &vector
impl<'f, FS: FieldSignature> std::ops::AddAssign<&Vector<'f, FS>> for Vector<'f, FS> {
    fn add_assign(&mut self, other: &Vector<'f, FS>) {
        match common_space(self.ambient_space, other.ambient_space) {
            Some(space) => {
                let n = space.borrow().linear_dimension().unwrap();
                for i in 0..n {
                    space
                        .field()
                        .add_mut(self.coordinate_mut(i), other.coordinate(i));
                }
            }
            None => panic!("Can't add vectors belonging to different spaces"),
        }
    }
}

// &vector - &vector
impl<'f, FS: FieldSignature> std::ops::Sub<&Vector<'f, FS>> for &Vector<'f, FS> {
    type Output = Vector<'f, FS>;

    fn sub(self, other: &Vector<'f, FS>) -> Self::Output {
        self + &(-other)
    }
}

// &vector * &scalar
impl<'f, FS: FieldSignature> Vector<'f, FS> {
    pub fn scalar_mul(&self, other: &FS::Set) -> Vector<'f, FS> {
        Vector {
            ambient_space: self.ambient_space,
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ambient_space().field().mul(x, other))
                .collect(),
        }
    }
}

// It is helpful for computational reasons to put an ordering on the vectors
// so that the points of a simplex can be ordered
impl<'f, FS: OrderedRingSignature + FieldSignature> PartialOrd for Vector<'f, FS> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl<'f, FS: OrderedRingSignature + FieldSignature> Ord for Vector<'f, FS> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        for i in 0..space.linear_dimension().unwrap() {
            match space
                .field()
                .ring_cmp(self.coordinate(i), other.coordinate(i))
            {
                std::cmp::Ordering::Less => {
                    return std::cmp::Ordering::Less;
                }
                std::cmp::Ordering::Equal => {}
                std::cmp::Ordering::Greater => {
                    return std::cmp::Ordering::Greater;
                }
            }
        }
        std::cmp::Ordering::Equal
    }
}

pub fn vectors_from_rows<'f, FS: FieldSignature + 'f>(
    sp: AffineSpace<'f, FS>,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<'f, FS>> {
    assert_eq!(mat.cols(), sp.linear_dimension().unwrap());
    (0..mat.rows())
        .map(|r| Vector::new(sp, (0..mat.cols()).map(|c| mat.at(r, c).unwrap().clone())))
        .collect()
}

pub fn vectors_from_cols<'f, FS: FieldSignature + 'f>(
    sp: AffineSpace<'f, FS>,
    mat: &Matrix<FS::Set>,
) -> Vec<Vector<'f, FS>> {
    assert_eq!(mat.rows(), sp.linear_dimension().unwrap());
    vectors_from_rows(sp, &mat.transpose_ref())
}

pub fn vector_from_row<'f, FS: FieldSignature + 'f>(
    sp: AffineSpace<'f, FS>,
    mat: &Matrix<FS::Set>,
) -> Vector<'f, FS> {
    assert_eq!(mat.rows(), 1);
    assert_eq!(mat.cols(), sp.linear_dimension().unwrap());
    vectors_from_rows(sp, mat).pop().unwrap()
}

pub fn vector_from_col<'f, FS: FieldSignature + 'f>(
    sp: AffineSpace<'f, FS>,
    mat: &Matrix<FS::Set>,
) -> Vector<'f, FS> {
    assert_eq!(mat.rows(), sp.linear_dimension().unwrap());
    assert_eq!(mat.cols(), 1);
    vector_from_row(sp, &mat.transpose_ref())
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;

    #[test]
    fn vector_from_mat() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let mat = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(1), Rational::from(2)],
            vec![Rational::from(3), Rational::from(4)],
        ]);

        mat.pprint();

        let mut vecs = vectors_from_rows(space, &mat);
        let v2 = vecs.pop().unwrap();
        let v1 = vecs.pop().unwrap();
        println!("v1 = {v1:?}");
        println!("v2 = {v2:?}");

        assert_eq!(
            v1,
            Vector::new(space, vec![Rational::from(1), Rational::from(2)])
        );
        assert_eq!(
            v2,
            Vector::new(space, vec![Rational::from(3), Rational::from(4)])
        );
    }

    #[test]
    fn det() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space, vec![Rational::from(3), Rational::from(2)]);
        let v2 = Vector::new(space, vec![Rational::from(5), Rational::from(7)]);
        assert_eq!(space.determinant(vec![&v1, &v2]), Rational::from(11));
    }

    #[test]
    fn test_abgroup() {
        let space_ab = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let a = Vector::new(space_ab, vec![Rational::from(1), Rational::from(2)]);
        let b = Vector::new(space_ab, vec![Rational::from(6), Rational::from(3)]);
        let c = Vector::new(space_ab, vec![Rational::from(7), Rational::from(5)]);

        let space_xy = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let x = Vector::new(space_xy, vec![Rational::from(1), Rational::from(2)]);
        let y = Vector::new(space_xy, vec![Rational::from(6), Rational::from(3)]);
        let z = Vector::new(space_xy, vec![Rational::from(7), Rational::from(5)]);
        let w = Vector::new(space_xy, vec![Rational::from(-2), Rational::from(-4)]);

        assert_eq!(c, &a + &b);
        assert_eq!(z, &x + &y);
        assert_eq!(a, a);
        assert_ne!(a, b);
        assert_ne!(a, x); //same coordinates but different space
        assert_ne!(x.scalar_mul(&Rational::from(-2)), z);
        assert_eq!(x.scalar_mul(&Rational::from(-2)), w);
    }
}
