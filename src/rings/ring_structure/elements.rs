use std::{
    fmt::{Debug, Display},
    ops::{Add, Div, Mul, Neg, Sub},
    rc::Rc,
};

use malachite_nz::integer::Integer;

use super::super::structure::*;
use super::structure::*;

impl<RS: RingStructure> PartialEq for StructuredElement<RS> {
    fn eq(&self, other: &Self) -> bool {
        let structure = common_structure(self.structure(), other.structure());
        structure.equal(&self.elem(), &other.elem())
    }
}

impl<RS: RingStructure> Eq for StructuredElement<RS> {}

impl<RS: RingStructure> Neg for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn neg(self) -> Self::Output {
        StructuredElement::new(self.structure().clone(), self.structure().neg(&self.elem()))
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
            self.structure().add(&self.elem(), &rhs.elem()),
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
                .add(&self.elem(), &self.structure().neg(&rhs.elem())),
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
            self.structure().mul(&self.elem(), &rhs.elem()),
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

impl<RS: RingStructure> Div<&StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        StructuredElement::new(
            self.structure().clone(),
            self.structure().div(&self.elem(), &rhs.elem()).unwrap(),
        )
    }
}

impl<RS: RingStructure> Div<StructuredElement<RS>> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        self / &rhs
    }
}

impl<RS: RingStructure> Div<&StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: &StructuredElement<RS>) -> Self::Output {
        &self / rhs
    }
}

impl<RS: RingStructure> Div<StructuredElement<RS>> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn div(self, rhs: StructuredElement<RS>) -> Self::Output {
        &self / &rhs
    }
}

impl<RS: RingStructure> StructuredElement<RS> {
    pub fn pow(&self, n: i32) -> StructuredElement<RS> {
        StructuredElement::new(
            self.structure().clone(),
            self.structure()
                .int_pow(&self.elem(), &Integer::from(n))
                .unwrap(),
        )
    }
}

//adding i32
impl<RS: RingStructure> Add<i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Add<&i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Add<i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self + StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(rhs)),
        )
    }
}

impl<RS: RingStructure> Add<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self + StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(*rhs)),
        )
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) + rhs
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) + rhs
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) + rhs
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) + rhs
    }
}

//subbing i32
impl<RS: RingStructure> Sub<i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Sub<&i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Sub<i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self - StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(rhs)),
        )
    }
}

impl<RS: RingStructure> Sub<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self - StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(*rhs)),
        )
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) - rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) - rhs
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) - rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) - rhs
    }
}

//multiplying i32
impl<RS: RingStructure> Mul<i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Mul<&i32> for StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Mul<i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: i32) -> Self::Output {
        let ring = self.structure().clone();
        self * StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(rhs)),
        )
    }
}

impl<RS: RingStructure> Mul<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self * StructuredElement::new(
            self.structure().clone(),
            self.structure().from_int(&Integer::from(*rhs)),
        )
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) * rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(self)),
        ) * rhs
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) * rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(
            rhs.structure().clone(),
            rhs.structure().from_int(&Integer::from(*self)),
        ) * rhs
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::super::super::number::integer::*;
    use super::super::super::polynomial::polynomial::*;
    use super::super::cannonical::*;

    use super::*;

    #[test]
    fn test_poly_elem_opps() {
        let x = &Polynomial::<Integer>::var().as_elem();

        let f = 4 * x.pow(2) - 1;
        let g = 2 * x + 1;

        //test operations by value and reference
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
