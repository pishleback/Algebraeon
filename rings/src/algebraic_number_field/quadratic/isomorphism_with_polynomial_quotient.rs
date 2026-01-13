use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldPolynomialQuotientStructure, QuadraticNumberFieldElement,
        QuadraticNumberFieldStructure,
    },
    polynomial::Polynomial,
    structure::{
        AdditiveMonoidSignature, CharZeroFieldSignature, FieldSignature, MetaFactoringMonoid,
        MetaMultiplicativeMonoidSignature, RingHomomorphism, SemiModuleSignature,
    },
};
use algebraeon_nzq::{Integer, Natural, Rational, traits::DivMod};
use algebraeon_sets::structure::{
    BijectiveFunction, BorrowedSet, BorrowedStructure, Function, InjectiveFunction, Morphism,
};

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldIsomorphism<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> {
    anf_polyquo: AlgebraicNumberFieldPolynomialQuotientStructureBorrowed,
    anf_quadratic: QuadraticNumberFieldStructure<Integer>,
    // where does the generator of anf_polyquo get mapped to inside anf_quadratic?
    generator_image: QuadraticNumberFieldElement,
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    #[allow(unused)]
    fn new(
        anf_polyquo: AlgebraicNumberFieldPolynomialQuotientStructureBorrowed,
    ) -> Result<Self, ()> {
        let anf_polyquo_borrowed = anf_polyquo.borrow();
        if anf_polyquo_borrowed.degree() == 2 {
            // let g be the generator of this ANF, so we are QQ[g]
            // g is a root of an integer polynomial
            // ax^2 + bx + c
            let poly = anf_polyquo_borrowed.modulus().primitive_part_fof();
            let poly_coeffs = poly.clone().into_coeffs();
            debug_assert_eq!(poly_coeffs.len(), 3);
            let two_a = Integer::TWO * &poly_coeffs[2];
            let neg_b = -&poly_coeffs[1];

            // find s, d such that s^2d = b^2-4ac and d is squarefree
            let (s, d) = {
                let disc = poly.discriminant().unwrap();
                let (mut d, powers) = disc.factor().into_unit_and_powers().unwrap();
                let mut s = Integer::ONE;
                for (p, k) in powers {
                    let (q, r) = k.div_mod(Natural::TWO);
                    s *= p.nat_pow(&q);
                    if r == Natural::ONE {
                        d *= p;
                    }
                }
                (s, d)
            };

            // now g = (-b ± s sqrt(d)) / 2a
            //       = (-b/2a) + (±s/2a) sqrt(d)
            // by the quadratic formula
            // this determines the desired isomorphism
            Ok(Self {
                anf_polyquo,
                anf_quadratic: QuadraticNumberFieldStructure::new_unchecked(d),
                generator_image: QuadraticNumberFieldElement {
                    rational_part: Rational::from_integers(neg_b, two_a.clone()),
                    algebraic_part: Rational::from_integers(s, two_a),
                },
            })
        } else {
            Err(())
        }
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> Morphism<AlgebraicNumberFieldPolynomialQuotientStructure, QuadraticNumberFieldStructure<Integer>>
    for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn domain(&self) -> &AlgebraicNumberFieldPolynomialQuotientStructure {
        self.anf_polyquo.borrow()
    }

    fn range(&self) -> &QuadraticNumberFieldStructure<Integer> {
        &self.anf_quadratic
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> Function<AlgebraicNumberFieldPolynomialQuotientStructure, QuadraticNumberFieldStructure<Integer>>
    for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn image(&self, x: &Polynomial<Rational>) -> QuadraticNumberFieldElement {
        let x = self.anf_polyquo.borrow().to_vec(x);
        debug_assert!(x.len() == 2);
        self.anf_quadratic.add(
            &self.anf_quadratic.from_rat(&x[0]),
            &self
                .anf_quadratic
                .inbound_principal_rational_map()
                .range_module_structure()
                .scalar_mul(&self.generator_image, &x[1]),
        )
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
>
    InjectiveFunction<
        AlgebraicNumberFieldPolynomialQuotientStructure,
        QuadraticNumberFieldStructure<Integer>,
    > for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn try_preimage(&self, y: &QuadraticNumberFieldElement) -> Option<Polynomial<Rational>> {
        Some(self.preimage(y))
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
>
    BijectiveFunction<
        AlgebraicNumberFieldPolynomialQuotientStructure,
        QuadraticNumberFieldStructure<Integer>,
    > for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn preimage(&self, y: &QuadraticNumberFieldElement) -> Polynomial<Rational> {
        let gen_coeff = &y.algebraic_part / &self.generator_image.algebraic_part;
        let rat_coeff = &y.rational_part - &gen_coeff * &self.generator_image.rational_part;
        Polynomial::from_coeffs(vec![rat_coeff, gen_coeff])
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
>
    RingHomomorphism<
        AlgebraicNumberFieldPolynomialQuotientStructure,
        QuadraticNumberFieldStructure<Integer>,
    > for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
}

impl AlgebraicNumberFieldPolynomialQuotientStructure {
    /// Returns an isomorphism from this degree 2 polynomial quotient number field representation to the quadratic number field implementation .
    /// Returns `Err` if this number field is not quadratic.
    pub fn quadratic_anf_isomorphism(&self) -> Result<QuadraticNumberFieldIsomorphism<&Self>, ()> {
        QuadraticNumberFieldIsomorphism::new(self)
    }

    /// Returns an isomorphism from the quadratic number field implementation to this degree 2 polynomial quotient number field representation.
    /// Returns `Err` if this number field is not quadratic.
    pub fn into_quadratic_anf_isomorphism(
        self,
    ) -> Result<QuadraticNumberFieldIsomorphism<Self>, ()> {
        QuadraticNumberFieldIsomorphism::new(self)
    }
}

impl<D: BorrowedSet<Integer>> QuadraticNumberFieldStructure<D> {}

#[cfg(test)]
mod tests {
    use algebraeon_sets::structure::EqSignature;

    use crate::parsing::parse_rational_polynomial;

    use super::*;

    // ZZ[i]
    #[test]
    fn qanf_neg_1() {
        let poly = parse_rational_polynomial("x^2-6*x+13", "x").unwrap();
        let anf_poly = poly.algebraic_number_field().unwrap();
        let iso = anf_poly.quadratic_anf_isomorphism().unwrap();
        let anf_quad = iso.range();

        assert!(anf_quad.equal(
            &iso.image(&parse_rational_polynomial("1", "x").unwrap()),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(1),
                algebraic_part: Rational::from(0)
            }
        ));
        assert!(anf_quad.equal(
            &iso.image(&parse_rational_polynomial("x", "x").unwrap()),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(3),
                algebraic_part: Rational::from(2)
            }
        ));
        assert!(anf_quad.equal(
            &iso.image(&parse_rational_polynomial("x^2", "x").unwrap()),
            &QuadraticNumberFieldElement {
                rational_part: Rational::from(5),
                algebraic_part: Rational::from(12)
            }
        ));

        assert!(anf_poly.equal(
            &iso.preimage(&QuadraticNumberFieldElement {
                rational_part: Rational::from(1),
                algebraic_part: Rational::from(0)
            }),
            &parse_rational_polynomial("1", "x").unwrap()
        ));

        assert!(anf_poly.equal(
            &iso.preimage(&QuadraticNumberFieldElement {
                rational_part: Rational::from(0),
                algebraic_part: Rational::from(2)
            }),
            &parse_rational_polynomial("x-3", "x").unwrap()
        ));
    }
}
