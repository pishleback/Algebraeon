use super::conway_polynomials::conway_polynomial;
use crate::{
    linear::matrix::{Matrix, MatrixStructure},
    polynomial::*,
    rings::quotient::QuotientStructure,
    structure::*,
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Natural};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConwayFiniteFieldStructure {
    p: usize,
    n: usize,
    structure: FieldExtensionByPolynomialQuotientStructure<
        QuotientStructure<IntegerCanonicalStructure, true>,
    >,
}

impl ConwayFiniteFieldStructure {
    pub fn new(p: usize, n: usize) -> Result<Self, ()> {
        let f = conway_polynomial(p, n)?;
        Ok(Self {
            p,
            n,
            structure: FieldExtensionByPolynomialQuotientStructure::new_field_unchecked(
                PolynomialStructure::new(QuotientStructure::new_field_unchecked(
                    Integer::structure(),
                    Integer::from(p),
                )),
                f.clone(),
            ),
        })
    }

    pub fn reduce(&self, f: Polynomial<Integer>) -> Polynomial<Integer> {
        self.structure.reduce(f)
    }

    pub fn to_row_vector(&self, f: &Polynomial<Integer>) -> Matrix<Integer> {
        self.structure.to_row_vector(f)
    }

    pub fn to_col_vector(&self, f: &Polynomial<Integer>) -> Matrix<Integer> {
        self.structure.to_col_vector(f)
    }

    pub fn to_vector(&self, f: &Polynomial<Integer>) -> Vec<Integer> {
        self.structure.to_vector(f)
    }

    pub fn from_row_vector(&self, v: Matrix<Integer>) -> Polynomial<Integer> {
        self.structure.from_row_vector(v)
    }

    pub fn from_col_vector(&self, v: Matrix<Integer>) -> Polynomial<Integer> {
        self.structure.from_col_vector(v)
    }

    pub fn from_vector(&self, v: Vec<Integer>) -> Polynomial<Integer> {
        self.structure.from_vector(v)
    }
}

impl Signature for ConwayFiniteFieldStructure {}

impl SetSignature for ConwayFiniteFieldStructure {
    type Set = Polynomial<Integer>;

    fn is_element(&self, _: &Self::Set) -> bool {
        true
    }
}

impl EqSignature for ConwayFiniteFieldStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.structure.equal(a, b)
    }
}

impl SemiRingSignature for ConwayFiniteFieldStructure {
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

impl CharacteristicSignature for ConwayFiniteFieldStructure {
    fn characteristic(&self) -> Natural {
        self.characteristic_and_power().0
    }
}

impl RingSignature for ConwayFiniteFieldStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.structure.neg(a)
    }
}

impl UnitsSignature for ConwayFiniteFieldStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        self.structure.inv(a)
    }
}

impl IntegralDomainSignature for ConwayFiniteFieldStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        self.structure.div(a, b)
    }
}

impl FieldSignature for ConwayFiniteFieldStructure {}

impl FiniteUnitsSignature for ConwayFiniteFieldStructure {
    fn all_units(&self) -> Vec<Self::Set> {
        self.structure.all_units()
    }
}

impl FiniteFieldSignature for ConwayFiniteFieldStructure {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.p.into(), self.n.into())
    }
}

#[derive(Debug, Clone)]
pub struct ConwayFiniteFieldInclusion {
    // a prime
    p: usize,
    // finite field of order p^m
    domain: ConwayFiniteFieldStructure,
    // finite field of order p^n
    range: ConwayFiniteFieldStructure,
    // n/m
    degree: usize,
    // Linear map F_{p^m} -> F_{p^n} of column vectors of polynomial coefficients over F_p
    inclusion: Matrix<Integer>,
    // matricies modulo p
    mat_mod_p: MatrixStructure<QuotientStructure<IntegerCanonicalStructure, true>>,
}

impl ConwayFiniteFieldInclusion {
    pub fn new(p: usize, m: usize, n: usize) -> Result<Self, ()> {
        if n % m == 0 {
            let degree = n / m;

            let domain = ConwayFiniteFieldStructure::new(p, m).unwrap();
            let range = ConwayFiniteFieldStructure::new(p, n).unwrap();

            // r = (p^n-1)/(p^m-1)
            let r = (Natural::from(p).pow(&n.into()) - Natural::ONE)
                / (Natural::from(p).pow(&m.into()) - Natural::ONE);

            // x -> x^r
            let inclusion = Matrix::join_cols(
                n,
                (0..m)
                    .map(|i| range.to_col_vector(&range.nat_pow(&Polynomial::var_pow(i), &r)))
                    .collect(),
            );
            assert_eq!(inclusion.rows(), n);
            assert_eq!(inclusion.cols(), m);
            Ok(Self {
                p,
                domain,
                range,
                degree,
                inclusion,
                mat_mod_p: MatrixStructure::new(QuotientStructure::new_field(
                    Integer::structure(),
                    p.into(),
                )),
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
        self.range().from_col_vector(
            self.mat_mod_p
                .mul(&self.inclusion, &self.domain().to_col_vector(x))
                .unwrap(),
        )
    }
}

impl InjectiveFunction<ConwayFiniteFieldStructure, ConwayFiniteFieldStructure>
    for ConwayFiniteFieldInclusion
{
    fn try_preimage(&self, x: &Polynomial<Integer>) -> Option<Polynomial<Integer>> {
        Some(
            self.domain().from_vector(
                self.inclusion
                    .clone()
                    .col_solve(&self.range().to_vector(x))?,
            ),
        )
    }
}

impl RingHomomorphism<ConwayFiniteFieldStructure, ConwayFiniteFieldStructure>
    for ConwayFiniteFieldInclusion
{
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

    #[test]
    fn conway_finite_fields_homomorphism() {
        let f = ConwayFiniteFieldInclusion::new(2, 2, 4).unwrap();

        f.inclusion.pprint();

        let a = Polynomial::var();
        println!("a = {}", a);
        let b = f.image(&a.clone());
        println!("b = {}", b);

        assert!(
            f.domain()
                .equal(&f.domain().nat_pow(&a, &3u32.into()), &f.domain().one())
        );
        assert!(
            f.range()
                .equal(&f.range().nat_pow(&b, &3u32.into()), &f.range().one())
        );

        let a2 = f.try_preimage(&b).unwrap();
        println!("a2 = {}", a2);
        assert!(f.domain().equal(&a, &a2));
    }

    #[test]
    fn conway_finite_fields_compatability() {
        /*
                 x
        f_{3^2} ===> f_{3^4}
           ||           ||
         y ||           || z
           \/           \/
        f_{3^6} ===> f_{3^12}
                 w

        */

        let f_3_2 = ConwayFiniteFieldStructure::new(3, 2).unwrap();
        let f_3_12 = ConwayFiniteFieldStructure::new(3, 12).unwrap();

        let x = ConwayFiniteFieldInclusion::new(3, 2, 4).unwrap();
        let y = ConwayFiniteFieldInclusion::new(3, 2, 6).unwrap();
        let z = ConwayFiniteFieldInclusion::new(3, 4, 12).unwrap();
        let w = ConwayFiniteFieldInclusion::new(3, 6, 12).unwrap();

        for a in f_3_2.all_elements() {
            assert!(f_3_12.equal(&z.image(&x.image(&a)), &w.image(&y.image(&a))));
        }
    }
}
