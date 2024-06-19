use std::borrow::Borrow;

use crate::rings::linear::matrix::Matrix;

use super::*;

// #[derive(Debug, Clone)]
// pub struct Point<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> {
//     abmient_space: &'a Space<FS>,
//     coordinates: Vec<FS::Set>, //length equal to abmient_space.dimension()
// }

// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> PartialEq for Point<FS, SP> {
//     fn eq(&self, other: &Self) -> bool {
//         let space = common_space(&self.abmient_space, &other.abmient_space);
//         let n = space.dimension();
//         (0..n).all(|i| {
//             space
//                 .ordered_field()
//                 .equal(self.coordinate(i), other.coordinate(i))
//         })
//     }
// }

// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> Eq for Point<FS, SP> {}

#[derive(Debug, Clone)]
pub struct Vector<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> {
    abmient_space: SP,
    coordinates: Vec<FS::Set>, //length equal to abmient_space.dimension()
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> PartialEq
    for Vector<FS, SP>
{
    fn eq(&self, other: &Self) -> bool {
        match common_space(self.abmient_space.borrow(), other.abmient_space.borrow()) {
            Some(space) => {
                let n = space.dimension();
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

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> Eq for Vector<FS, SP> {}

// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> Point<FS, SP> {
//     pub fn new(abmient_space: &'a Space<FS>, coordinates: Vec<FS::Set>) -> Self {
//         assert_eq!(abmient_space.dimension(), coordinates.len());
//         Self {
//             abmient_space,
//             coordinates,
//         }
//     }

//     fn into_vector(self) -> Vector<FS, SP> {
//         Vector {
//             abmient_space: self.abmient_space,
//             coordinates: self.coordinates,
//         }
//     }

//     pub fn ordered_field(&self) -> Rc<FS> {
//         self.abmient_space.ordered_field()
//     }

//     pub fn coordinate(&self, i: usize) -> &FS::Set {
//         self.coordinates.get(i).unwrap()
//     }
// }

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>>> Vector<FS, SP> {
    pub fn new(abmient_space: SP, coordinates: Vec<FS::Set>) -> Self {
        assert_eq!(abmient_space.borrow().dimension(), coordinates.len());
        Self {
            abmient_space,
            coordinates,
        }
    }

    // fn into_point(self) -> Point<FS, SP> {
    //     Point {
    //         abmient_space: self.abmient_space,
    //         coordinates: self.coordinates,
    //     }
    // }

    pub fn abmient_space(&self) -> &SP {
        &self.abmient_space
    }

    pub fn ordered_field(&self) -> Rc<FS> {
        self.abmient_space.borrow().ordered_field()
    }

    pub fn dimension(&self) -> usize {
        self.abmient_space.borrow().dimension()
    }

    pub fn coordinate(&self, i: usize) -> &FS::Set {
        self.coordinates.get(i).unwrap()
    }

    pub fn into_row(&self) -> Matrix<FS::Set> {
        Matrix::construct(self.dimension(), 1, |r, _c| self.coordinate(r).clone())
    }

    pub fn into_col(&self) -> Matrix<FS::Set> {
        Matrix::construct(1, self.dimension(), |_r, c| self.coordinate(c).clone())
    }
}

// -vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>> + Clone> std::ops::Neg
    for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn neg(self) -> Self::Output {
        Vector {
            abmient_space: self.abmient_space.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ordered_field().neg(x))
                .collect(),
        }
    }
}

// vector + vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>> + Clone>
    std::ops::Add<&Vector<FS, SP>> for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn add(self, other: &Vector<FS, SP>) -> Self::Output {
        match common_space(self.abmient_space.clone(), other.abmient_space.clone()) {
            Some(space) => {
                let n = space.borrow().dimension();
                let coordinates = (0..n)
                    .map(|i| {
                        space
                            .borrow()
                            .ordered_field()
                            .add(self.coordinate(i), other.coordinate(i))
                    })
                    .collect();
                Vector {
                    abmient_space: space,
                    coordinates,
                }
            }
            None => panic!("Can't add vectors belonging to different spaces"),
        }
    }
}

// vector - vector
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>> + Clone>
    std::ops::Sub<&Vector<FS, SP>> for &Vector<FS, SP>
{
    type Output = Vector<FS, SP>;

    fn sub(self, other: &Vector<FS, SP>) -> Self::Output {
        self + &(-other)
    }
}

// vector * scalar
impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<LinearSpace<FS>> + Clone>
    Vector<FS, SP>
{
    fn scalar_mul(&self, other: &FS::Set) -> Vector<FS, SP> {
        Vector {
            abmient_space: self.abmient_space.clone(),
            coordinates: self
                .coordinates
                .iter()
                .map(|x| self.ordered_field().mul(x, other))
                .collect(),
        }
    }
}

// //point - vector -> point
// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> std::ops::Sub<&Vector<FS, SP>>
//     for &Point<FS, SP>
// {
//     type Output = Point<FS, SP>;

//     fn sub(self, other: &Vector<FS, SP>) -> Self::Output {
//         (&self.clone().into_vector() - other).into_point()
//     }
// }

// //point + vector -> point
// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> std::ops::Add<&Vector<FS, SP>>
//     for &Point<FS, SP>
// {
//     type Output = Point<FS, SP>;

//     fn add(self, other: &Vector<FS, SP>) -> Self::Output {
//         (&self.clone().into_vector() + other).into_point()
//     }
// }

// //point - point -> vector
// impl<FS: OrderedRingStructure + FieldStructure, SP : Borrow<LinearSpace<FS>>> std::ops::Sub<&Point<FS, SP>>
//     for &Point<FS, SP>
// {
//     type Output = Vector<FS, SP>;

//     fn sub(self, other: &Point<FS, SP>) -> Self::Output {
//         &self.clone().into_vector() - &other.clone().into_vector()
//     }
// }

#[cfg(test)]
mod tests {
    use malachite_q::Rational;

    use crate::rings::structure::StructuredType;

    use super::*;

    #[test]
    fn test_abgroup() {
        let space_ab = LinearSpace::new(Rational::structure(), 2);
        let a = Vector::new(&space_ab, vec![Rational::from(1), Rational::from(2)]);
        let b = Vector::new(&space_ab, vec![Rational::from(6), Rational::from(3)]);
        let c = Vector::new(&space_ab, vec![Rational::from(7), Rational::from(5)]);

        let space_xy = LinearSpace::new(Rational::structure(), 2);
        let x = Vector::new(&space_xy, vec![Rational::from(1), Rational::from(2)]);
        let y = Vector::new(&space_xy, vec![Rational::from(6), Rational::from(3)]);
        let z = Vector::new(&space_xy, vec![Rational::from(7), Rational::from(5)]);
        let w = Vector::new(&space_xy, vec![Rational::from(-2), Rational::from(-4)]);

        assert_eq!(c, &a + &b);
        assert_eq!(z, &x + &y);
        assert_eq!(a, a);
        assert_ne!(a, b);
        assert!(std::panic::catch_unwind(|| a == x).is_err()); //same coordinates but different space
        assert_ne!(x.scalar_mul(&Rational::from(-2)), z);
        assert_eq!(x.scalar_mul(&Rational::from(-2)), w);
    }
}
