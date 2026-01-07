use std::{
    borrow::Borrow,
    ops::{Add, Div, Mul, Neg, Sub},
};

use super::*;
use algebraeon_sets::structure::*;

use algebraeon_nzq::*;

pub trait IntoErgonomicSignature: SetSignature {
    fn into_ergonomic(&self, elem: Self::Set) -> StructuredElement<Self> {
        StructuredElement::new(self.clone(), elem)
    }
}
impl<S: SetSignature> IntoErgonomicSignature for S {}

pub trait IntoErgonomic: MetaType {
    fn into_ergonomic(self) -> StructuredElement<Self::Signature> {
        StructuredElement::new(Self::structure(), self)
    }
}
impl<T: MetaType> IntoErgonomic for T {}

fn common_structure<S: Signature>(structure1: impl Borrow<S>, structure2: impl Borrow<S>) -> S {
    if structure1.borrow() == structure2.borrow() {
        structure1.borrow().clone()
    } else {
        panic!("Unequal ring structures")
    }
}

#[derive(Debug, Clone)]
pub struct StructuredElement<S: SetSignature> {
    structure: S,
    elem: S::Set,
}

impl<S: SetSignature> StructuredElement<S> {
    pub fn new(structure: S, elem: S::Set) -> Self {
        Self { structure, elem }
    }

    pub fn structure(&self) -> S {
        self.structure.clone()
    }

    pub fn ref_set(&self) -> &S::Set {
        &self.elem
    }

    pub fn into_verbose(self) -> S::Set {
        self.elem
    }
}

impl<S: SetSignature> std::fmt::Display for StructuredElement<S>
where
    S::Set: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.elem, f)
    }
}

impl<RS: RingEqSignature> PartialEq for StructuredElement<RS> {
    fn eq(&self, other: &Self) -> bool {
        let structure = common_structure::<RS>(self.structure(), other.structure());
        structure.equal(self.ref_set(), other.ref_set())
    }
}

impl<RS: RingEqSignature> Eq for StructuredElement<RS> {}

impl<RS: RingSignature> Neg for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn neg(self) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().neg(self.ref_set()),
        )
    }
}

impl<RS: RingSignature> Neg for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn neg(self) -> Self::Output {
        (&self).neg()
    }
}

impl<RS: RingSignature> Add<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().add(self.ref_set(), rhs.ref_set()),
        )
    }
}

impl<RS: RingSignature> Add<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        self + &rhs
    }
}

impl<RS: RingSignature> Add<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self + rhs
    }
}

impl<RS: RingSignature> Add<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self + &rhs
    }
}

impl<RS: RingSignature> Sub<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .add(self.ref_set(), &self.structure().neg(rhs.ref_set())),
        )
    }
}

impl<RS: RingSignature> Sub<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        self - &rhs
    }
}

impl<RS: RingSignature> Sub<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self - rhs
    }
}

impl<RS: RingSignature> Sub<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self - &rhs
    }
}

impl<RS: RingSignature> Mul<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().mul(self.ref_set(), rhs.ref_set()),
        )
    }
}

impl<RS: RingSignature> Mul<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        self * &rhs
    }
}

impl<RS: RingSignature> Mul<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self * rhs
    }
}

impl<RS: RingSignature> Mul<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self * &rhs
    }
}

impl<RS: IntegralDomainSignature> Div<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .try_div(self.ref_set(), rhs.ref_set())
                .unwrap(),
        )
    }
}

impl<RS: IntegralDomainSignature> Div<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        self / &rhs
    }
}

impl<RS: IntegralDomainSignature> Div<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self / rhs
    }
}

impl<RS: IntegralDomainSignature> Div<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self / &rhs
    }
}

impl<RS: IntegralDomainSignature> StructuredElement<RS> {
    pub fn pow(&self, n: i32) -> StructuredElement<RS> {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .try_int_pow(self.ref_set(), &Integer::from(n))
                .unwrap(),
        )
    }
}

macro_rules! impl_int_ops {
    ($I : ty) => {
        //adding $I

        impl<RS: RingSignature> Add<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Add<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Add<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self + StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Add<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &$I) -> Self::Output {
                let ring = self.structure();
                self + StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Add<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) + rhs
            }
        }

        impl<RS: RingSignature> Add<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) + rhs
            }
        }

        impl<RS: RingSignature> Add<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
                    + rhs
            }
        }

        impl<RS: RingSignature> Add<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
                    + rhs
            }
        }

        //subbing $I
        impl<RS: RingSignature> Sub<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Sub<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Sub<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Sub<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self - StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Sub<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) - rhs
            }
        }

        impl<RS: RingSignature> Sub<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) - rhs
            }
        }

        impl<RS: RingSignature> Sub<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
                    - rhs
            }
        }

        impl<RS: RingSignature> Sub<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
                    - rhs
            }
        }

        //multiplying $I
        impl<RS: RingSignature> Mul<$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Mul<&$I> for StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Mul<$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: $I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(ring.clone(), ring.from_int(Integer::from(rhs)))
            }
        }

        impl<RS: RingSignature> Mul<&$I> for &StructuredElement<RS> {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &$I) -> Self::Output {
                let ring = self.structure().clone();
                self * StructuredElement::new(
                    ring.clone(),
                    ring.from_int(Integer::from(rhs.clone())),
                )
            }
        }

        impl<RS: RingSignature> Mul<StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) * rhs
            }
        }

        impl<RS: RingSignature> Mul<&StructuredElement<RS>> for $I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self))) * rhs
            }
        }

        impl<RS: RingSignature> Mul<StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
                    * rhs
            }
        }

        impl<RS: RingSignature> Mul<&StructuredElement<RS>> for &$I {
            type Output = StructuredElement<RS>;

            fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
                let ring = rhs.structure().clone();
                StructuredElement::new(ring.clone(), ring.from_int(Integer::from(self.clone())))
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

    use crate::polynomial::*;

    use super::*;

    #[test]
    fn test_poly_elem_operations() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

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
