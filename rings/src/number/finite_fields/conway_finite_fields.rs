use super::conway_polynomials::conway_polynomial;
use crate::{
    polynomial::{FieldExtensionByPolynomialQuotientStructure, Polynomial, PolynomialStructure},
    structure::{
        FieldStructure, FiniteFieldStructure, FiniteUnitsStructure, IntegralDomainStructure,
        QuotientStructure, RingStructure, SemiRingStructure, UnitsStructure,
    },
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Natural, traits::DivMod};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConwayFiniteFieldStructure {
    p: Natural,
    n: Natural,
    structure: FieldExtensionByPolynomialQuotientStructure<
        QuotientStructure<IntegerCanonicalStructure, true>,
    >,
}

impl ConwayFiniteFieldStructure {
    pub fn new(p: Natural, n: Natural) -> Result<Self, ()> {
        let f = conway_polynomial(&p, &n)?;
        Ok(Self {
            structure: FieldExtensionByPolynomialQuotientStructure::new_field_unchecked(
                PolynomialStructure::new(QuotientStructure::new_field_unchecked(
                    Integer::structure(),
                    Integer::from(&p),
                )),
                f.clone(),
            ),
            p,
            n,
        })
    }

    pub fn reduce(&self, f: Polynomial<Integer>) -> Polynomial<Integer> {
        self.structure.reduce(f)
    }
}

impl Structure for ConwayFiniteFieldStructure {}

impl SetStructure for ConwayFiniteFieldStructure {
    type Set = Polynomial<Integer>;

    fn is_element(&self, _: &Self::Set) -> bool {
        true
    }
}

impl EqStructure for ConwayFiniteFieldStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.structure.equal(a, b)
    }
}

impl SemiRingStructure for ConwayFiniteFieldStructure {
    fn zero(&self) -> Self::Set {
        self.structure.zero()
    }

    fn one(&self) -> Self::Set {
        self.structure.one()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.structure.add(a, b)
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.structure.mul(a, b)
    }
}

impl RingStructure for ConwayFiniteFieldStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.structure.neg(a)
    }
}

impl UnitsStructure for ConwayFiniteFieldStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        self.structure.inv(a)
    }
}

impl IntegralDomainStructure for ConwayFiniteFieldStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        self.structure.div(a, b)
    }
}

impl FieldStructure for ConwayFiniteFieldStructure {}

impl FiniteUnitsStructure for ConwayFiniteFieldStructure {
    fn all_units(&self) -> Vec<Self::Set> {
        self.structure.all_units()
    }
}

impl FiniteFieldStructure for ConwayFiniteFieldStructure {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.p.clone(), self.n.clone())
    }
}

#[derive(Debug, Clone)]
pub struct ConwayFiniteFieldInclusion {
    // a prime
    p: Natural,
    // finite field of order p^m
    domain: ConwayFiniteFieldStructure,
    // finite field of order p^n
    range: ConwayFiniteFieldStructure,
    // n/m
    degree: Natural,
    // (p^n-1)/(p^m-1)
    r: Natural,
}

impl ConwayFiniteFieldInclusion {
    pub fn new(p: Natural, m: Natural, n: Natural) -> Result<Self, ()> {
        let (degree, r) = (&n).div_mod(&m);
        if r == Natural::ZERO {
            Ok(Self {
                r: ((&p).pow(&n) - Natural::ONE) / ((&p).pow(&m) - Natural::ONE),
                domain: ConwayFiniteFieldStructure::new(p.clone(), m).unwrap(),
                range: ConwayFiniteFieldStructure::new(p.clone(), n).unwrap(),
                p,
                degree,
            })
        } else {
            Err(())
        }
    }
}

impl Morphism<ConwayFiniteFieldStructure, ConwayFiniteFieldStructure>
    for ConwayFiniteFieldInclusion
{
    fn domain(&self) -> &ConwayFiniteFieldStructure {
        &self.domain
    }

    fn range(&self) -> &ConwayFiniteFieldStructure {
        &self.range
    }
}

impl Function<ConwayFiniteFieldStructure, ConwayFiniteFieldStructure>
    for ConwayFiniteFieldInclusion
{
    fn image(&self, x: &Polynomial<Integer>) -> Polynomial<Integer> {
        self.range()
            .reduce(x.eval_var_pow(self.r.clone().try_into().unwrap()))
    }
}

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct ConwayFiniteFieldInclusionStructure {
//     p: Natural,
// }
// impl ConwayFiniteFieldInclusionStructure {
//     pub fn new(p: Natural) -> Self {
//         debug_assert!(is_prime(&p));
//         Self { p }
//     }
// }

// impl Structure for ConwayFiniteFieldInclusionStructure {}

// impl SetStructure for ConwayFiniteFieldInclusionStructure {
//     type Set = ConwayFiniteFieldInclusion;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         self.p == x.p
//     }
// }

// impl MorphismsStructure<ConwayFiniteFieldStructure, ConwayFiniteFieldStructure>
//     for ConwayFiniteFieldInclusionStructure
// {
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomicStructure;

    #[test]
    fn conway_finite_fields() {
        let f = ConwayFiniteFieldInclusion::new(2u32.into(), 2u32.into(), 4u32.into()).unwrap();

        let a = f.domain().into_ergonomic(Polynomial::var());
        let b = f.range().into_ergonomic(f.image(&a.clone().into_verbose()));

        assert!(
            f.domain()
                .equal(&a.pow(3).into_verbose(), &f.domain().one())
        );
        assert!(f.range().equal(&b.pow(3).into_verbose(), &f.range().one()));
    }
}
