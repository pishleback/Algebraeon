use algebraeon_rings::matrix::Matrix;
use std::borrow::Borrow;
use std::hash::Hash;
use super::*;

#[derive(Clone)]
pub struct Vector<'f, FS: FieldSignature + 'f> {
    ambient_space: AffineSpace<'f, FS>,
    coordinates: Vec<FS::Set>, //length equal to ambient_space.dimension()
}

#[allow(clippy::missing_fields_in_debug)]
impl<'f, FS: FieldSignature> std::fmt::Debug for Vector<'f, FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Vector")
            .field("coordinates", &self.coordinates)
            .finish()
    }
}

impl<'f, FS: FieldSignature> PartialEq for Vector<'f, FS> {
    fn eq(&self, other: &Self) -> bool {
        match common_space(&self.ambient_space, &other.ambient_space) {
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
    pub fn new(ambient_space: AffineSpace<'f, FS>, coordinates: Vec<FS::Set>) -> Self {
        assert_eq!(ambient_space.linear_dimension().unwrap(), coordinates.len());
        Self {
            ambient_space,
            coordinates,
        }
    }

    pub fn construct(
        ambient_space: AffineSpace<'f, FS>,
        mut coordinate_func: impl FnMut(usize) -> FS::Set,
    ) -> Self {
        #[allow(clippy::redundant_closure)]
        let coordinates = (0..ambient_space.borrow().linear_dimension().unwrap())
            .map(|i| coordinate_func(i))
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

    pub fn ambient_space(&self) -> &AffineSpace<'f, FS> {
        &self.ambient_space
    }

    // pub fn ordered_field(&self) -> Rc<FS> {
    //     self.ambient_space.borrow().ordered_field()
    // }

    // pub fn dimension(&self) -> usize {
    //     self.ambient_space.borrow().dimension()
    // }

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
            ambient_space: self.ambient_space.clone(),
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
        match common_space(&self.ambient_space, &other.ambient_space) {
            Some(space) => {
                let n = space.linear_dimension().unwrap();
                let coordinates = (0..n)
                    .map(|i| space.field().add(self.coordinate(i), other.coordinate(i)))
                    .collect();
                Vector {
                    ambient_space: space.clone(),
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
        match common_space(&self.ambient_space, &other.ambient_space) {
            Some(space) => {
                let space = space.clone();
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
            ambient_space: self.ambient_space.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ambient_space().field().mul(x, other))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Rational;

    use super::*;

    #[test]
    fn vector_from_mat() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let mat = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(1), Rational::from(2)],
            vec![Rational::from(3), Rational::from(4)],
        ]);

        mat.pprint();

        let mut vecs = vectors_from_rows(&space, &mat);
        let v2 = vecs.pop().unwrap();
        let v1 = vecs.pop().unwrap();
        println!("v1 = {:?}", v1);
        println!("v2 = {:?}", v2);

        assert_eq!(
            v1,
            Vector::new(space.clone(), vec![Rational::from(1), Rational::from(2)])
        );
        assert_eq!(
            v2,
            Vector::new(space.clone(), vec![Rational::from(3), Rational::from(4)])
        );
    }

    #[test]
    fn det() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space.clone(), vec![Rational::from(3), Rational::from(2)]);
        let v2 = Vector::new(space.clone(), vec![Rational::from(5), Rational::from(7)]);
        assert_eq!(space.determinant(vec![&v1, &v2]), Rational::from(11));
    }

    #[test]
    fn test_abgroup() {
        let space_ab = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let a = Vector::new(space_ab.clone(), vec![Rational::from(1), Rational::from(2)]);
        let b = Vector::new(space_ab.clone(), vec![Rational::from(6), Rational::from(3)]);
        let c = Vector::new(space_ab.clone(), vec![Rational::from(7), Rational::from(5)]);

        let space_xy = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let x = Vector::new(space_xy.clone(), vec![Rational::from(1), Rational::from(2)]);
        let y = Vector::new(space_xy.clone(), vec![Rational::from(6), Rational::from(3)]);
        let z = Vector::new(space_xy.clone(), vec![Rational::from(7), Rational::from(5)]);
        let w = Vector::new(
            space_xy.clone(),
            vec![Rational::from(-2), Rational::from(-4)],
        );

        assert_eq!(c, &a + &b);
        assert_eq!(z, &x + &y);
        assert_eq!(a, a);
        assert_ne!(a, b);
        assert_ne!(a, x); //same coordinates but different space
        assert_ne!(x.scalar_mul(&Rational::from(-2)), z);
        assert_eq!(x.scalar_mul(&Rational::from(-2)), w);
    }
}
