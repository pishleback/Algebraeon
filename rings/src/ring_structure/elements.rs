use std::{ops::{Add, Div, Mul, Neg, Sub}, rc::Rc};

use algebraeon_structure::*;
use malachite_nz::integer::Integer;
use super::structure::*;

pub trait IntoRingElem: MetaType {
    fn into_ring(self) -> StructuredElement<Self::Structure> {
        StructuredElement::new(Self::structure(), self)
    }
}
impl<T: MetaType> IntoRingElem for T {}

#[derive(Debug, Clone)]
pub struct StructuredElement<S: Structure> {
    structure: Rc<S>,
    elem: S::Set,
}

impl<S: Structure> StructuredElement<S> {
    pub fn new(structure: Rc<S>, elem: S::Set) -> Self {
        Self { structure, elem }
    }

    pub fn structure(&self) -> Rc<S> {
        self.structure.clone()
    }

    pub fn ref_set(&self) -> &S::Set {
        &self.elem
    }

    pub fn into_set(self) -> S::Set {
        self.elem
    }
}

impl<S: Structure> std::fmt::Display for StructuredElement<S>
where
    S::Set: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.elem, f)
    }
}


impl<RS: RingStructure> PartialEq for StructuredElement<RS> {
    fn eq(&self, other: &Self) -> bool {
        let structure = common_structure(self.structure(), other.structure());
        structure.equal(&self.ref_set(), &other.ref_set())
    }
}

impl<RS: RingStructure> Eq for StructuredElement<RS> {}

impl<RS: RingStructure> Neg for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn neg(self) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().neg(&self.ref_set()),
        )
    }
}

impl<RS: RingStructure> Neg for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn neg(self) -> Self::Output {
        (&self).neg()
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().add(&self.ref_set(), &rhs.ref_set()),
        )
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        self + &rhs
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self + rhs
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self + &rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .add(&self.ref_set(), &self.structure().neg(&rhs.ref_set())),
        )
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        self - &rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self - rhs
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self - &rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().mul(&self.ref_set(), &rhs.ref_set()),
        )
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        self * &rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self * rhs
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self * &rhs
    }
}

impl<RS: IntegralDomainStructure> Div<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .div(&self.ref_set(), &rhs.ref_set())
                .unwrap(),
        )
    }
}

impl<RS: IntegralDomainStructure> Div<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        self / &rhs
    }
}

impl<RS: IntegralDomainStructure> Div<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self / rhs
    }
}

impl<RS: IntegralDomainStructure> Div<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self / &rhs
    }
}

impl<RS: IntegralDomainStructure> StructuredElement<RS> {
    pub fn pow(&self, n: i32) -> StructuredElement<RS> {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .int_pow(&self.ref_set(), &Integer::from(n))
                .unwrap(),
        )
    }
}

macro_rules! impl_int_ops {
    ($I : ty) => {
        //adding $I

        impl<RS: RingStructure> Add<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Add<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Add<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Add<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &$I) -> Self::Output {
                let ring = self.structure();
                self + StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Add<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) + rhs
            }
        }

        impl<RS: RingStructure> Add<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) + rhs
            }
        }

        impl<RS: RingStructure> Add<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    + rhs
            }
        }

        impl<RS: RingStructure> Add<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    + rhs
            }
        }

        //subbing $I
        impl<RS: RingStructure> Sub<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Sub<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Sub<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Sub<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Sub<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) - rhs
            }
        }

        impl<RS: RingStructure> Sub<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) - rhs
            }
        }

        impl<RS: RingStructure> Sub<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    - rhs
            }
        }

        impl<RS: RingStructure> Sub<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    - rhs
            }
        }

        //multiplying $I
        impl<RS: RingStructure> Mul<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Mul<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Mul<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
            }
        }

        impl<RS: RingStructure> Mul<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(
                    ring.clone(),
                    ring.from_int(&Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingStructure> Mul<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) * rhs
            }
        }

        impl<RS: RingStructure> Mul<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) * rhs
            }
        }

        impl<RS: RingStructure> Mul<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    * rhs
            }
        }

        impl<RS: RingStructure> Mul<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self.clone())))
                    * rhs
            }
        }
    };
}

// impl_int_ops!(u8);
// impl_int_ops!(u16);
// impl_int_ops!(u32);
// impl_int_ops!(u64);
// impl_int_ops!(u128);

// impl_int_ops!(i8);
// impl_int_ops!(i16);
impl_int_ops!(i32);
// impl_int_ops!(i64);
// impl_int_ops!(i128);

impl_int_ops!(Integer);

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::super::super::polynomial::polynomial::*;

    use super::*;

    #[test]
    fn test_poly_elem_opps() {
        let x = &Polynomial::<Integer>::var().into_ring();

        let f = 4 * x.pow(2) - 1;
        let g = 2 * x + 1;

        // test operations by value and reference
        println!("{}", -&f.clone());
        println!("{}", -f.clone());
        println!("{}", f.clone() + g.clone());
        println!("{}", &f + g.clone());
        println!("{}", f.clone() + &g);
        println!("{}", &f + &g);
        println!("{}", f.clone() - g.clone());
        println!("{}", &f - g.clone());
        println!("{}", f.clone() - &g);
        println!("{}", &f - &g);
        println!("{}", f.clone() * g.clone());
        println!("{}", &f * g.clone());
        println!("{}", f.clone() * &g);
        println!("{}", &f * &g);
        println!("{}", f.clone() / g.clone());
        println!("{}", &f / g.clone());
        println!("{}", f.clone() / &g);
        println!("{}", &f / &g);
    }
}
