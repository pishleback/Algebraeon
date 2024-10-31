use orthoclase_rings::linear::matrix::Matrix;
use std::borrow::Borrow;
use std::hash::Hash;

use super::*;

#[derive(Clone)]
pub struct Vector<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> {
    ambient_space: SP,
    coordinates: Vec<FS::Set>, //length equal to ambient_space.dimension()
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> std::fmt::Debug for Vector<FS, SP> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Vector").field("coordinates", &self.coordinates).finish()
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> PartialEq
    for Vector<FS, SP>
{
    fn eq(&self, other: &Self) -> bool {
        match common_space(self.ambient_space.borrow(), other.ambient_space.borrow()) {
            Some(space) => {
                let n = space.linear_dimension().unwrap();
                (0..n).all(|i| {
                    space
                        .ordered_field()
                        .equal(self.coordinate(i), other.coordinate(i))
                })
            }
            None => false,
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> Eq for Vector<FS, SP> {}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> Hash for Vector<FS, SP>
where
    FS::Set: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // self.ambient_space.borrow().hash(state);
        self.coordinates.hash(state);
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>>> Vector<FS, SP> {
    pub fn new(ambient_space: SP, coordinates: Vec<FS::Set>) -> Self {
        assert_eq!(
            ambient_space.borrow().linear_dimension().unwrap(),
            coordinates.len()
        );
        Self {
            ambient_space,
            coordinates,
        }
    }

    pub fn construct(ambient_space: SP, mut coordinate_func: impl FnMut(usize) -> FS::Set) -> Self {
        let coordinates = (0..ambient_space.borrow().linear_dimension().unwrap())
            .map(|i| coordinate_func(i))
            .collect();
        Self {
            ambient_space,
            coordinates,
        }
    }

    pub fn zero(ambient_space: SP) -> Self {
        let ordered_field = ambient_space.borrow().ordered_field();
        Self::construct(ambient_space, |i| ordered_field.zero())
    }

    pub fn ambient_space(&self) -> &SP {
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

    pub fn into_row(&self) -> Matrix<FS::Set> {
        Matrix::construct(
            1,
            self.ambient_space().borrow().linear_dimension().unwrap(),
            |_r, c| self.coordinate(c).clone(),
        )
    }

    pub fn into_col(&self) -> Matrix<FS::Set> {
        Matrix::construct(
            self.ambient_space().borrow().linear_dimension().unwrap(),
            1,
            |r, _c| self.coordinate(r).clone(),
        )
    }
}

// -&vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> std::ops::Neg
    for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn neg(self) -> Self::Output {
        Vector {
            ambient_space: self.ambient_space.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ambient_space().borrow().ordered_field().neg(x))
                .collect(),
        }
    }
}

// &vector + &vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    std::ops::Add<&Vector<FS, SP>> for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn add(self, other: &Vector<FS, SP>) -> Self::Output {
        match common_space(self.ambient_space.clone(), other.ambient_space.clone()) {
            Some(space) => {
                let n = space.borrow().linear_dimension().unwrap();
                let coordinates = (0..n)
                    .map(|i| {
                        space
                            .borrow()
                            .ordered_field()
                            .add(self.coordinate(i), other.coordinate(i))
                    })
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
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    std::ops::AddAssign<&Vector<FS, SP>> for Vector<FS, SP>
{
    fn add_assign(&mut self, other: &Vector<FS, SP>) {
        match common_space(self.ambient_space.clone(), other.ambient_space.clone()) {
            Some(space) => {
                let n = space.borrow().linear_dimension().unwrap();
                for i in 0..n {
                    space
                        .borrow()
                        .ordered_field()
                        .add_mut(self.coordinate_mut(i), other.coordinate(i));
                }
            }
            None => panic!("Can't add vectors belonging to different spaces"),
        }
    }
}

// &vector - &vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    std::ops::Sub<&Vector<FS, SP>> for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn sub(self, other: &Vector<FS, SP>) -> Self::Output {
        self + &(-other)
    }
}

// &vector * &scalar
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    Vector<FS, SP>
{
    pub fn scalar_mul(&self, other: &FS::Set) -> Vector<FS, SP> {
        Vector {
            ambient_space: self.ambient_space.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ambient_space().borrow().ordered_field().mul(x, other))
                .collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use malachite_q::Rational;

    use orthoclase_rings::structure::StructuredType;

    use super::*;

    #[test]
    fn vector_from_mat() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
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
            Vector::new(&space, vec![Rational::from(1), Rational::from(2)])
        );
        assert_eq!(
            v2,
            Vector::new(&space, vec![Rational::from(3), Rational::from(4)])
        );
    }

    #[test]
    fn det() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = Vector::new(&space, vec![Rational::from(3), Rational::from(2)]);
        let v2 = Vector::new(&space, vec![Rational::from(5), Rational::from(7)]);
        assert_eq!(space.determinant(vec![&v1, &v2]), Rational::from(11));
    }

    #[test]
    fn test_abgroup() {
        let space_ab = AffineSpace::new_linear(Rational::structure(), 2);
        let a = Vector::new(&space_ab, vec![Rational::from(1), Rational::from(2)]);
        let b = Vector::new(&space_ab, vec![Rational::from(6), Rational::from(3)]);
        let c = Vector::new(&space_ab, vec![Rational::from(7), Rational::from(5)]);

        let space_xy = AffineSpace::new_linear(Rational::structure(), 2);
        let x = Vector::new(&space_xy, vec![Rational::from(1), Rational::from(2)]);
        let y = Vector::new(&space_xy, vec![Rational::from(6), Rational::from(3)]);
        let z = Vector::new(&space_xy, vec![Rational::from(7), Rational::from(5)]);
        let w = Vector::new(&space_xy, vec![Rational::from(-2), Rational::from(-4)]);

        assert_eq!(c, &a + &b);
        assert_eq!(z, &x + &y);
        assert_eq!(a, a);
        assert_ne!(a, b);
        assert_ne!(a, x); //same coordinates but different space
        assert_ne!(x.scalar_mul(&Rational::from(-2)), z);
        assert_eq!(x.scalar_mul(&Rational::from(-2)), w);
    }
}
