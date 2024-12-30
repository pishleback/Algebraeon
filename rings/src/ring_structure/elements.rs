use std::ops::{Add, Div, Mul, Neg, Sub};

use algebraeon_structure::*;
use malachite_nz::integer::Integer;

use super::super::elements::*;
use super::structure::*;

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
        self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Add<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &i32) -> Self::Output {
        let ring = self.structure();
        self + StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) + rhs
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) + rhs
    }
}

impl<RS: RingStructure> Add<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) + rhs
    }
}

impl<RS: RingStructure> Add<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn add(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) + rhs
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
        self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Sub<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self - StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) - rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) - rhs
    }
}

impl<RS: RingStructure> Sub<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) - rhs
    }
}

impl<RS: RingStructure> Sub<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn sub(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) - rhs
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
        self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(rhs)))
    }
}

impl<RS: RingStructure> Mul<&i32> for &StructuredElement<RS> {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &i32) -> Self::Output {
        let ring = self.structure().clone();
        self * StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*rhs)))
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) * rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(self))) * rhs
    }
}

impl<RS: RingStructure> Mul<StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) * rhs
    }
}

impl<RS: RingStructure> Mul<&StructuredElement<RS>> for &i32 {
    type Output = StructuredElement<RS>;

    fn mul(self, rhs: &StructuredElement<RS>) -> Self::Output {
        let ring = rhs.structure().clone();
        StructuredElement::new(ring.clone(), ring.from_int(&Integer::from(*self))) * rhs
    }
}

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
