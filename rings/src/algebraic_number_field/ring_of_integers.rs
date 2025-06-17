use super::polynomial_quotient_number_field::AlgebraicNumberFieldPolynomialQuotientStructure;
use crate::{
    algebraic_number_field::number_field_traits::AlgebraicNumberFieldSignature,
    matrix::Matrix,
    module::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    polynomial::Polynomial,
    structure::*,
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Natural, Rational};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct RingOfIntegersWithIntegralBasisStructure {
    algebraic_number_field: AlgebraicNumberFieldPolynomialQuotientStructure,
    integral_basis: Vec<Polynomial<Rational>>,
    discriminant: Integer,
    // The below just aid in the efficiency of calculations
    one: Option<Vec<Integer>>, // store 1
    mul_crossterms: Option<Vec<Vec<Vec<Integer>>>>,
}

impl ToStringSignature for RingOfIntegersWithIntegralBasisStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        self.roi_to_anf(elem).to_string()
    }
}

impl PartialEq for RingOfIntegersWithIntegralBasisStructure {
    fn eq(&self, other: &Self) -> bool {
        self.algebraic_number_field == other.algebraic_number_field
            && self.integral_basis == other.integral_basis
            && self.discriminant == other.discriminant
    }
}
impl Eq for RingOfIntegersWithIntegralBasisStructure {}

impl RingOfIntegersWithIntegralBasisStructure {
    pub fn new(
        algebraic_number_field: AlgebraicNumberFieldPolynomialQuotientStructure,
        integral_basis: Vec<Polynomial<Rational>>,
        discriminant: Integer,
    ) -> Self {
        debug_assert_eq!(integral_basis.len(), algebraic_number_field.degree());
        debug_assert_eq!(
            algebraic_number_field
                .rational_extension()
                .discriminant(&integral_basis),
            discriminant
        );
        let (_, true_discriminant) =
            algebraic_number_field.compute_integral_basis_and_discriminant();
        debug_assert_eq!(discriminant, true_discriminant);
        for a in &integral_basis {
            debug_assert!(algebraic_number_field.is_algebraic_integer(a));
        }
        let mut roi = Self {
            algebraic_number_field,
            integral_basis,
            discriminant,
            one: None,
            mul_crossterms: None,
        };
        let n = roi.degree();
        roi.one = Some(
            roi.try_anf_to_roi(&roi.algebraic_number_field.one())
                .unwrap(),
        );
        roi.mul_crossterms = Some(
            (0..n)
                .map(|j| {
                    (0..=j)
                        .map(|i| {
                            roi.try_anf_to_roi(
                                &roi.algebraic_number_field
                                    .mul(&roi.integral_basis[i], &roi.integral_basis[j]),
                            )
                            .unwrap()
                        })
                        .collect()
                })
                .collect(),
        );
        roi
    }

    pub fn degree(&self) -> usize {
        debug_assert_eq!(
            self.integral_basis.len(),
            self.algebraic_number_field.degree()
        );
        self.integral_basis.len()
    }

    pub fn z_module(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, IntegerCanonicalStructure> {
        Integer::structure().into_free_module(self.degree())
    }

    pub fn integral_basis(&self) -> &Vec<Polynomial<Rational>> {
        &self.integral_basis
    }

    pub fn basis_element(&self, i: usize) -> &Polynomial<Rational> {
        assert!(i < self.degree());
        &self.integral_basis[i]
    }

    pub fn anf(&self) -> &AlgebraicNumberFieldPolynomialQuotientStructure {
        &self.algebraic_number_field
    }
}

// impl RingOfIntegersWithIntegralBasisElement {
//     pub fn basis_element(n: usize, i: usize) -> Self {
//         Self {
//             coefficients: (0..n)
//                 .map(|j| if i == j { Integer::ONE } else { Integer::ZERO })
//                 .collect(),
//         }
//     }

//     pub fn into_col(self) -> Matrix<Integer> {
//         Matrix::from_cols(vec![self.coefficients])
//     }

//     pub fn from_col(m: &Matrix<Integer>) -> Self {
//         debug_assert_eq!(m.cols(), 1);
//         let n = m.rows();
//         Self {
//             coefficients: (0..n).map(|i| m.at(i, 0).unwrap().clone()).collect(),
//         }
//     }

//     pub fn into_row(self) -> Matrix<Integer> {
//         Matrix::from_cols(vec![self.coefficients])
//     }

//     pub fn from_row(m: &Matrix<Integer>) -> Self {
//         debug_assert_eq!(m.rows(), 1);
//         let n = m.cols();
//         Self {
//             coefficients: (0..n).map(|i| m.at(0, i).unwrap().clone()).collect(),
//         }
//     }

//     pub fn into_coefficients(self) -> Vec<Integer> {
//         self.coefficients
//     }

//     pub fn coefficients(&self) -> &Vec<Integer> {
//         &self.coefficients
//     }

//     pub fn from_coefficients(coefficients: Vec<Integer>) -> Self {
//         Self {
//             coefficients: coefficients,
//         }
//     }

//     pub fn scalar_mul(self, a: &Integer) -> Self {
//         Self {
//             coefficients: self.coefficients.into_iter().map(|c| c * a).collect(),
//         }
//     }

//     pub fn scalar_mul_ref(&self, a: &Integer) -> Self {
//         Self {
//             coefficients: self.coefficients.iter().map(|c| c * a).collect(),
//         }
//     }

//     pub fn neg(self) -> Self {
//         Self {
//             coefficients: self.coefficients.into_iter().map(|c| -c).collect(),
//         }
//     }

//     pub fn neg_ref(&self) -> Self {
//         Self {
//             coefficients: self.coefficients.iter().map(|c| -c).collect(),
//         }
//     }
// }

impl RingOfIntegersWithIntegralBasisStructure {
    pub fn roi_to_anf(&self, elem: &Vec<Integer>) -> Polynomial<Rational> {
        debug_assert!(self.is_element(elem).is_ok());
        let n = self.degree();
        self.algebraic_number_field.sum(
            (0..n)
                .map(|i| self.integral_basis[i].mul_scalar(&Rational::from(&elem[i])))
                .collect(),
        )
    }

    pub fn try_anf_to_roi(&self, elem: &Polynomial<Rational>) -> Option<Vec<Integer>> {
        let n = self.degree();
        let y = self.algebraic_number_field.to_vec(elem);
        let m = Matrix::join_cols(
            n,
            (0..n)
                .map(|i| self.algebraic_number_field.to_col(&self.integral_basis[i]))
                .collect(),
        );
        if let Some(s) = m.col_solve(&y) {
            debug_assert_eq!(s.len(), n);
            #[allow(clippy::manual_ok_err)]
            if let Ok(coefficients) = (0..n).map(|i| Integer::try_from(s[i].clone())).collect() {
                Some(coefficients)
            } else {
                None
            }
        } else {
            // the integral basis is a basis of the anf as a rational vector space
            unreachable!()
        }
    }
}

impl Signature for RingOfIntegersWithIntegralBasisStructure {}

impl SetSignature for RingOfIntegersWithIntegralBasisStructure {
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        if x.len() != self.degree() {
            return Err("wrong length".to_string());
        }
        Ok(())
    }
}

impl EqSignature for RingOfIntegersWithIntegralBasisStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl AdditiveMonoidSignature for RingOfIntegersWithIntegralBasisStructure {
    fn zero(&self) -> Self::Set {
        self.z_module().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.z_module().add(a, b)
    }
}

impl AdditiveGroupSignature for RingOfIntegersWithIntegralBasisStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.z_module().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.z_module().sub(a, b)
    }
}

impl SemiRingSignature for RingOfIntegersWithIntegralBasisStructure {
    fn one(&self) -> Self::Set {
        match &self.one {
            Some(one) => one.clone(),
            None => self
                .try_anf_to_roi(&self.algebraic_number_field.one())
                .unwrap(),
        }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let n = self.degree();
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match &self.mul_crossterms {
            Some(mul_crossterms) => {
                // Used cached cross-terms
                let mut t = self.zero();
                for i in 0..n {
                    for j in 0..n {
                        let c = &a[i] * &b[j];
                        // Sort i and j for indexing into mul_crossterms which is a triangular array
                        let (i, j) = if j <= i { (i, j) } else { (j, i) };
                        t = self.add(&t, &mul_crossterms[i][j].iter().map(|x| x * &c).collect());
                    }
                }
                debug_assert!(
                    self.equal(
                        &self
                            .try_anf_to_roi(
                                &self
                                    .algebraic_number_field
                                    .mul(&self.roi_to_anf(a), &self.roi_to_anf(b)),
                            )
                            .unwrap(),
                        &t
                    )
                );
                t
            }
            None => {
                // Compute using anf mul
                self.try_anf_to_roi(
                    &self
                        .algebraic_number_field
                        .mul(&self.roi_to_anf(a), &self.roi_to_anf(b)),
                )
                .unwrap()
            }
        }
    }
}

impl RingSignature for RingOfIntegersWithIntegralBasisStructure {}

impl CharacteristicSignature for RingOfIntegersWithIntegralBasisStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl SemiRingUnitsSignature for RingOfIntegersWithIntegralBasisStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        #[allow(clippy::collapsible_else_if)]
        if self.is_zero(a) {
            Err(RingDivisionError::DivideByZero)
        } else {
            if let Some(a_inv) = self.try_anf_to_roi(
                &self
                    .algebraic_number_field
                    .inv(&self.roi_to_anf(a))
                    .unwrap(),
            ) {
                Ok(a_inv)
            } else {
                Err(RingDivisionError::NotDivisible)
            }
        }
    }
}

impl IntegralDomainSignature for RingOfIntegersWithIntegralBasisStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        match self.inv(b) {
            Ok(b_inv) => Ok(self.mul(a, &b_inv)),
            Err(err) => Err(err),
        }
    }
}

impl DedekindDomainSignature for RingOfIntegersWithIntegralBasisStructure {}

impl CharZeroRingSignature for RingOfIntegersWithIntegralBasisStructure {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        self.anf().try_to_int(&self.roi_to_anf(x))
    }
}

#[derive(Debug, Clone)]
pub struct RingOfIntegersToAlgebraicNumberFieldInclusion {
    roi: RingOfIntegersWithIntegralBasisStructure,
    anf: AlgebraicNumberFieldPolynomialQuotientStructure,
}

impl RingOfIntegersToAlgebraicNumberFieldInclusion {
    pub fn from_algebraic_number_field(
        anf: AlgebraicNumberFieldPolynomialQuotientStructure,
    ) -> Self {
        Self {
            roi: anf.compute_ring_of_integers(),
            anf,
        }
    }

    pub fn from_ring_of_integers(roi: RingOfIntegersWithIntegralBasisStructure) -> Self {
        Self {
            anf: roi.anf().clone(),
            roi,
        }
    }
}

impl
    Morphism<
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    > for RingOfIntegersToAlgebraicNumberFieldInclusion
{
    fn domain(&self) -> &RingOfIntegersWithIntegralBasisStructure {
        &self.roi
    }

    fn range(&self) -> &AlgebraicNumberFieldPolynomialQuotientStructure {
        &self.anf
    }
}

impl
    Function<
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    > for RingOfIntegersToAlgebraicNumberFieldInclusion
{
    fn image(
        &self,
        x: &<RingOfIntegersWithIntegralBasisStructure as SetSignature>::Set,
    ) -> <AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set {
        self.roi.roi_to_anf(x)
    }
}

impl
    InjectiveFunction<
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    > for RingOfIntegersToAlgebraicNumberFieldInclusion
{
    fn try_preimage(
        &self,
        x: &<AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set,
    ) -> Option<<RingOfIntegersWithIntegralBasisStructure as SetSignature>::Set> {
        self.roi.try_anf_to_roi(x)
    }
}

impl
    RingHomomorphism<
        RingOfIntegersWithIntegralBasisStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    > for RingOfIntegersToAlgebraicNumberFieldInclusion
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomic;

    #[test]
    fn ring_of_integer_arithmetic() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Take the integral basis (0 + x, 1/2 + 1/2x)
        let a = Polynomial::<Rational>::from_coeffs(vec![Rational::ZERO, Rational::ONE]);
        let b = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE_HALF, Rational::ONE_HALF]);

        let anf = (x.pow(2) + 7)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = RingOfIntegersWithIntegralBasisStructure::new(
            anf.clone(),
            vec![a.clone(), b.clone()],
            Integer::from(-7),
        );

        {
            assert_eq!(
                roi.roi_to_anf(&vec![Integer::from(1), Integer::from(4)]),
                (2 + 3 * &x).into_verbose()
            );
        }

        {
            assert!(
                roi.try_anf_to_roi(&Polynomial::<Rational>::from_coeffs(vec![
                    Rational::ONE_HALF,
                    Rational::ONE,
                ]))
                .is_none()
            );

            let c = roi
                .try_anf_to_roi(&Polynomial::<Rational>::from_coeffs(vec![
                    Rational::from(2),
                    Rational::from(3),
                ]))
                .unwrap();
            assert_eq!(c.len(), 2);
            assert_eq!(c[0], Integer::from(1));
            assert_eq!(c[1], Integer::from(4));
        }

        {
            // 0 = 0 * (0+x) + 0 * (1/2 + 1/2x)
            let zero = roi.zero();
            assert_eq!(zero.len(), 2);
            assert_eq!(zero[0], Integer::ZERO);
            assert_eq!(zero[1], Integer::ZERO);
        }

        {
            // 1 = -1 * (0+x) + 2 * (1/2 + 1/2x)
            let one = roi.one();
            assert_eq!(one.len(), 2);
            assert_eq!(one[0], Integer::from(-1));
            assert_eq!(one[1], Integer::from(2));
        }

        {
            let alpha = roi.try_anf_to_roi(&(2 + 3 * &x).into_verbose()).unwrap();
            let beta = roi.try_anf_to_roi(&(-1 + 2 * &x).into_verbose()).unwrap();

            {
                let gamma = roi.try_anf_to_roi(&(1 + 5 * &x).into_verbose()).unwrap();
                // (2 + 3x) + (-1 + 2x) = 1 + 5x
                assert_eq!(roi.add(&alpha, &beta), gamma);
            }

            {
                let gamma = roi.try_anf_to_roi(&(-44 + &x).into_verbose()).unwrap();
                // x^2 = -7 so
                // (2 + 3x) * (-1 + 2x) = -44 + x
                assert_eq!(roi.mul(&alpha, &beta), gamma);
            }

            {
                let gamma = roi.try_anf_to_roi(&(-2 - 3 * &x).into_verbose()).unwrap();
                // -(2 + 3x) = -2 - 3x
                assert_eq!(roi.neg(&alpha), gamma);
            }

            {
                assert_eq!(roi.inv(&roi.neg(&roi.one())).unwrap(), roi.neg(&roi.one()));
                assert_eq!(roi.inv(&roi.one()).unwrap(), roi.one());
                assert!(roi.inv(&alpha).is_err());
                assert!(roi.inv(&beta).is_err());
            }
        }

        println!("{:?}", roi);
    }

    #[test]
    fn ideal_factoring() {
        // Construct the ring of integers Z[i]
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(3) + x + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.compute_ring_of_integers();
        let roi_ideals = roi.ideals();

        // The ideal (27i - 9) in Z[i]
        let ideal =
            roi_ideals.principal_ideal(&roi.try_anf_to_roi(&(27 * x - 9).into_verbose()).unwrap());

        let (a, b) = roi_ideals.ideal_two_generators(&ideal);
        println!("I = ({}, {})", roi.to_string(&a), roi.to_string(&b));

        // Factor the ideal
        for (prime_ideal, power) in roi_ideals
            .factorizations()
            .into_powers(roi_ideals.factor_ideal(&ideal).unwrap())
        {
            let (a, b) = roi_ideals.ideal_two_generators(prime_ideal.ideal());
            println!(
                "P = ({}, {})    power = {power}",
                roi.to_string(&a),
                roi.to_string(&b)
            );
        }
    }
}
