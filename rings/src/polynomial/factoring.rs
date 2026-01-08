use super::{Polynomial, polynomial_structure::*};
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::{BorrowedStructure, SetSignature};

impl<
    RS: UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + GreatestCommonDivisorSignature
        + CharZeroRingSignature,
    RSB: BorrowedStructure<RS>,
> PolynomialStructure<RS, RSB>
where
    PolynomialStructure<RS, RSB>: SetSignature<Set = Polynomial<RS::Set>>
        + GreatestCommonDivisorSignature
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    /// Reduce a factorization problem for polynomials over a ring of characteristic 0 to a factorization of primitive squarefree polynomials over the ring
    //https://en.wikipedia.org/wiki/Square-free_polynomial#Yun's_algorithm
    pub fn factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(
        &self,
        f: Polynomial<RS::Set>,
        factor_coeff: impl Fn(&RS::Set) -> Factored<RS::Set, Natural>,
        primitive_sqfree_factorize: &impl Fn(
            Polynomial<RS::Set>,
        ) -> Factored<Polynomial<RS::Set>, Natural>,
    ) -> Factored<Polynomial<RS::Set>, Natural> {
        debug_assert!(!self.is_zero(&f));
        //look for a squarefree factorization of the form
        //  f = x * a_1^1 * a_2^2 * a_3^3 * ... * a_k^k
        //where each a_i is a squarefree primitive polynomial and x is an element of R

        let (content, prim) = self.factor_primitive(f.clone()).unwrap();
        let (content_unit, content_factors) =
            factor_coeff(&content).into_unit_and_powers().unwrap();

        let mut factors = self.factorizations().new_unit_and_powers_impl(
            Polynomial::constant(content_unit),
            content_factors
                .into_iter()
                .map(|(factor, power)| (Polynomial::constant(factor), power))
                .collect(),
        );

        let mut f = prim;
        let f_prime = self.derivative(f.clone());
        let mut i: usize = 1;
        let mut a = self.gcd_by_primitive_subresultant(f.clone(), f_prime.clone());
        let mut b = self.try_div(&f, &a).unwrap();
        let mut c = self.try_div(&f_prime, &a).unwrap();
        let mut d = self.add(&self.neg(&self.derivative(b.clone())), &c);

        while self.degree(&f).unwrap() != 0 {
            a = self.gcd(&b, &d);

            //a^i is a power of a squarefree factor of f
            self.factorizations().mul_mut(
                &mut factors,
                &self
                    .factorizations()
                    .pow(&primitive_sqfree_factorize(a.clone()), &Natural::from(i)),
            );
            f = self
                .try_div(&f, &self.nat_pow(&a, &Natural::from(i)))
                .unwrap();

            (b, c) = (self.try_div(&b, &a).unwrap(), self.try_div(&d, &a).unwrap());
            i += 1;
            d = self.add(&self.neg(&self.derivative(b.clone())), &c);
        }
        self.factorizations()
            .mul_mut(&mut factors, &self.factorizations().new_unit_impl(f));
        factors
    }
}

impl<
    RS: FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + GreatestCommonDivisorSignature
        + FiniteUnitsSignature,
    RSB: BorrowedStructure<RS>,
> PolynomialStructure<RS, RSB>
where
    PolynomialStructure<RS, RSB>: SetSignature<Set = Polynomial<RS::Set>>
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    fn factor_primitive_linear_part(
        &self,
        mut f: Polynomial<RS::Set>,
    ) -> (Factored<Polynomial<RS::Set>, Natural>, Polynomial<RS::Set>) {
        debug_assert!(self.is_primitive(f.clone()));

        /*
        If f is a primitive polynomial over R then any
        linear factor `(a+bx)` of `f(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n`
        is such that `a` divides `c0` and `b` divides `cn`
        so just factor c0 and cn and check all divisors
        */

        let mut linear_factors = self.factorizations().one();
        'seek_linear_factor: while self.degree(&f).unwrap() > 0 {
            let c0 = self.coeff(&f, 0);
            #[allow(clippy::redundant_else)]
            if self.coeff_ring().is_zero(c0.as_ref()) {
                //linear factor of x
                f = self.try_div(&f, &self.var()).unwrap();
                linear_factors = self.factorizations().mul(
                    &linear_factors,
                    &self.factorizations().new_irreducible_impl(&self.var()),
                );
                continue 'seek_linear_factor;
            } else {
                //look for linear factors of the form (a+bx)
                let c0fs = self.coeff_ring().factor(self.coeff(&f, 0).as_ref());
                let cnfs = self
                    .coeff_ring()
                    .factor(self.coeff(&f, self.degree(&f).unwrap()).as_ref());
                for a_assoc in self.coeff_ring().factorizations().divisors(&c0fs).unwrap() {
                    for u in self.coeff_ring().all_units() {
                        let a = self.coeff_ring().mul(&u, &a_assoc);
                        for b in self.coeff_ring().factorizations().divisors(&cnfs).unwrap() {
                            //a ranges over all divisors of c0
                            //b ranges over all divisors factors of cn up to associates
                            //try the linear factor (a+bx)
                            let lin = Polynomial::from_coeffs(vec![a.clone(), b]);
                            debug_assert!(!self.is_zero(&lin));
                            if let Some(new_f) = self.try_div(&f, &lin) {
                                f = new_f;
                                linear_factors = self.factorizations().mul(
                                    &linear_factors,
                                    &self.factorizations().new_irreducible_impl(&lin),
                                );
                                continue 'seek_linear_factor;
                            }
                        }
                    }
                }
            }

            break;
        }
        (linear_factors, f)
    }
}

impl<
    RS: UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + IntegralDomainSignature
        + CharZeroRingSignature
        + FiniteUnitsSignature
        + 'static,
    RSB: BorrowedStructure<RS>,
> PolynomialStructure<RS, RSB>
where
    PolynomialStructure<RS, RSB>: SetSignature<Set = Polynomial<RS::Set>>,
{
    fn find_factor_primitive_by_kroneckers_algorithm(
        &self,
        f: &Polynomial<RS::Set>,
        factor_coeff: impl Fn(&RS::Set) -> Factored<RS::Set, Natural>,
    ) -> FindFactorResult<Polynomial<RS::Set>> {
        /*
        Suppose we want to factor f(x) = 2 + x + x^2 + x^4 + x^5
        Assume it has a proper factor g(x). wlog g(x) has degree <= 2
        g(x) is determined by its value at 3 points, say at x=0, x=1, x=-1
        f(0)=2, f(1)=6, f(-1)=2     if one of these was zero, then we would have found a linear factor
        g(0) divides 2, g(1) divides 6, g(-1) divides 2
        there are finitely many possible values of g(0), g(1) and g(-1) which satisfy these
        infact there are 4*8*4=128 possible triples
        however, only 64 need to be checked as the other half are their negatives
        more abstractly, some possibilities can be avoided because we only care about g up to multiplication by a unit
         */
        let f_deg = self.degree(f).unwrap();
        if f_deg == 1 {
            //linear factor is irreducible
            FindFactorResult::Irreducible
        } else {
            let max_factor_degree = f_deg / 2;
            let mut f_points = vec![];
            let mut elem_gen = self.coeff_ring().generate_distinct_elements();
            //take more samples than necessary, then take the subset with the smallest number of divisors
            while f_points.len() < 3 * (max_factor_degree + 1) {
                //loop terminates because polynomial over integral domain has finitely many roots
                let x = elem_gen.next().unwrap();
                let y = self.evaluate(f, &x);
                if !self.coeff_ring().is_zero(&y) {
                    f_points.push((x, factor_coeff(&y)));
                }
            }

            //compute all factors of each y value. choose the y with the most divisors to only factor up to units
            f_points.sort_by_cached_key(|(_x, yf)| {
                self.coeff_ring().factorizations().count_divisors(yf)
            });
            let _ = f_points.split_off(max_factor_degree + 1);
            //possible_g_points is (x, possible_y_values)
            let all_possible_g_points: Vec<(RS::Set, Vec<RS::Set>)> = f_points
                .into_iter()
                .rev()
                .enumerate()
                .map(|(i, (x, yf))| {
                    let mut y_divs = vec![];
                    for d in self.coeff_ring().factorizations().divisors(&yf).unwrap() {
                        if i == 0 {
                            //take divisors up to associates for one, because we only care about g up to associates
                            y_divs.push(d);
                        } else {
                            //take _all_ divisors for the rest
                            for u in self.coeff_ring().all_units() {
                                y_divs.push(self.coeff_ring().mul(&u, &d));
                            }
                        }
                    }
                    (x, y_divs)
                })
                .collect();

            for possible_g_points in itertools::Itertools::multi_cartesian_product(
                all_possible_g_points
                    .into_iter()
                    .map(|(x, y_divs)| y_divs.into_iter().map(move |y_div| (x.clone(), y_div))),
            ) {
                // println!("{:?}", possible_g_points);
                if let Some(g) = self.interpolate_by_lagrange_basis(&possible_g_points)
                    && self.degree(&g).unwrap() >= 1
                {
                    //g is a possible proper divisor of f
                    debug_assert!(!self.is_zero(&g));
                    if let Some(h) = self.try_div(f, &g) {
                        //g really is a proper divisor of f
                        return FindFactorResult::Composite(g, h);
                    }
                }
            }
            //f is irreducible
            FindFactorResult::Irreducible
        }
    }
}

impl<
    RS: UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + GreatestCommonDivisorSignature
        + CharZeroRingSignature
        + FiniteUnitsSignature
        + 'static,
    RSB: BorrowedStructure<RS>,
> PolynomialStructure<RS, RSB>
where
    PolynomialStructure<RS, RSB>: SetSignature<Set = Polynomial<RS::Set>>
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    pub fn factorize_by_kroneckers_method(
        &self,
        f: Polynomial<RS::Set>,
        factor_coeff: impl Fn(&RS::Set) -> Factored<RS::Set, Natural>,
    ) -> Factored<Polynomial<RS::Set>, Natural> {
        if self.is_zero(&f) {
            Factored::Zero
        } else {
            let (g, f_prim) = self.factor_primitive(f).unwrap();
            let g_factored = factor_coeff(&g);
            let (g_unit, g_factors) = g_factored.into_unit_and_powers().unwrap();
            let g_factored = self.factorizations().new_unit_and_powers_impl(
                Polynomial::constant(g_unit),
                g_factors
                    .into_iter()
                    .map(|(g_factor, power)| (Polynomial::constant(g_factor), power))
                    .collect(),
            );
            let mut factored = factorize_by_find_factor(self, f_prim, &|f| {
                self.find_factor_primitive_by_kroneckers_algorithm(&f, |c| factor_coeff(c))
            });
            self.factorizations().mul_mut(&mut factored, &g_factored);
            factored
        }
    }
}

impl<
    RS: UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>
        + GreatestCommonDivisorSignature
        + CharZeroRingSignature
        + FiniteUnitsSignature
        + 'static,
    RSB: BorrowedStructure<RS>,
> PolynomialStructure<RS, RSB>
where
    PolynomialStructure<RS, RSB>: SetSignature<Set = Polynomial<RS::Set>>
        + GreatestCommonDivisorSignature
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    pub fn factorize_by_yuns_and_kroneckers_method(
        &self,
        f: &Polynomial<RS::Set>,
        factor_coeff: impl Fn(&RS::Set) -> Factored<RS::Set, Natural>,
    ) -> Factored<Polynomial<RS::Set>, Natural> {
        if self.is_zero(f) {
            Factored::Zero
        } else {
            self.factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(
                f.clone(),
                &factor_coeff,
                &|f| {
                    let big_ff = factorize_by_find_factor(self, f, &|f| {
                        self.find_factor_primitive_by_kroneckers_algorithm(&f, |c| factor_coeff(c))
                    });
                    #[allow(clippy::let_and_return)]
                    big_ff
                },
            )
        }
    }
}

// impl<R: MetaType> Polynomial<R>
// where
//     R::Signature: UniqueFactorizationDomainSignature
//         + GreatestCommonDivisorSignature
//         + CharZeroRingSignature
//         + FiniteUnitsSignature
//         + 'static,
//     PolynomialStructure<R::Signature, R::Signature>: SetSignature<Set = Polynomial<R>>
//         + GreatestCommonDivisorSignature
//         + UniqueFactorizationDomainSignature,
// {
//     pub fn factorize_by_kroneckers_method(
//         &self,
//         factor_coeff: impl Fn(&R) -> Option<FactoredRingElement<R>>,
//     ) -> Option<FactoredRingElement<Polynomial<R>>> {
//         Self::structure().factorize_by_yuns_and_kroneckers_method(self, factor_coeff)
//     }
// }

pub fn factorize_by_factorize_primitive_part<
    Ring: RingSignature,
    Field: FieldSignature,
    FieldB: BorrowedStructure<Field>,
    Fof: FieldOfFractionsInclusion<Ring, Field>,
>(
    fof_inclusion: &Fof,
    poly_ring: &PolynomialStructure<Field, FieldB>,
    f: &Polynomial<Field::Set>,
) -> Factored<Polynomial<Field::Set>, Natural>
where
    PolynomialStructure<Field, FieldB>: SetSignature<Set = Polynomial<Field::Set>>
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    for<'a> PolynomialStructure<Ring, &'a Ring>: SetSignature<Set = Polynomial<Ring::Set>>
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    Ring: GreatestCommonDivisorSignature,
{
    let (unit, prim) = factor_primitive_fof(fof_inclusion, f);
    if let Some((prim_unit, prim_factors)) = fof_inclusion
        .domain()
        .polynomials()
        .factor(&prim)
        .into_unit_and_powers()
    {
        let mut fof_unit = prim_unit.apply_map(|c| fof_inclusion.image(c));
        let mut fof_factors = vec![];
        for (factor, power) in prim_factors {
            let fof_factor = factor.apply_map(|c| fof_inclusion.image(c));
            let (fof_factor_unit, fof_factor_prim) = poly_ring.factor_fav_assoc(&fof_factor);
            poly_ring.mul_mut(&mut fof_unit, &fof_factor_unit);
            fof_factors.push((fof_factor_prim, power));
        }

        poly_ring.factorizations().new_unit_and_powers_impl(
            poly_ring.mul(&Polynomial::constant(unit), &fof_unit),
            fof_factors,
        )
    } else {
        Factored::Zero
    }
}

impl<RS: FieldSignature + FiniteUnitsSignature, RSB: BorrowedStructure<RS>>
    PolynomialStructure<RS, RSB>
where
    Self: SetSignature<Set = Polynomial<RS::Set>>
        + UniqueFactorizationMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
{
    fn find_factor_by_trying_all_factors(
        &self,
        f: Polynomial<RS::Set>,
    ) -> FindFactorResult<Polynomial<RS::Set>> {
        let f_deg = self.degree(&f).unwrap();
        let max_factor_degree = f_deg / 2;
        for d in 0..max_factor_degree {
            for mut coeffs in itertools::Itertools::multi_cartesian_product((0..=d).map(|_d| {
                let mut all_elems = vec![self.coeff_ring().zero()];
                all_elems.append(&mut self.coeff_ring().all_units());
                all_elems
            })) {
                coeffs.push(self.coeff_ring().one());
                let g = Polynomial::from_coeffs(coeffs);
                debug_assert!(!self.is_zero(&g));
                if let Some(h) = self.try_div(&f, &g) {
                    return FindFactorResult::Composite(g, h);
                }
            }
        }
        FindFactorResult::Irreducible
    }

    pub fn factorize_by_trying_all_factors(
        &self,
        f: Polynomial<RS::Set>,
    ) -> Factored<Polynomial<RS::Set>, Natural> {
        if self.is_zero(&f) {
            Factored::Zero
        } else {
            factorize_by_find_factor(self, f, &|f| self.find_factor_by_trying_all_factors(f))
        }
    }
}

// impl<F: MetaType> Polynomial<F>
// where
//     F::Signature: FieldSignature + FiniteUnitsSignature,
//     PolynomialStructure<F::Signature, F::Signature>:
//         SetSignature<Set = Polynomial<F>> + UniqueFactorizationDomainSignature,
// {
//     pub fn factorize_by_trying_all_factors(&self) -> Option<FactoredRingElement<Polynomial<F>>> {
//         Self::structure().factorize_by_trying_all_factors(self.clone())
//     }
// }

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomic;
    use algebraeon_sets::structure::{EqSignature, MetaType};

    #[test]
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let int_poly_fs = Integer::structure()
            .into_polynomials()
            .into_factorizations();

        //primitive cases
        let f = ((1 + x).pow(2)).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![((1 + x).into_verbose(), Natural::from(2u8))]
                )
            )
        );

        let f = (-1 - 2 * x).into_verbose();
        let fs1 = Integer::structure()
            .polynomials()
            .factorize_by_kroneckers_method(f, Integer::factor);
        let fs2 = &int_poly_fs.new_unit_and_powers_unchecked(
            Polynomial::neg(&Polynomial::one()),
            vec![((1 + 2 * x).into_verbose(), Natural::from(1u8))],
        );
        println!("fs1={:?} fs2={:?}", fs1, fs2);
        assert!(int_poly_fs.equal(&fs1, fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![
                        ((1 + x + x.pow(2)).into_verbose(), Natural::from(1u8)),
                        ((2 - x + x.pow(3)).into_verbose(), Natural::from(1u8))
                    ]
                )
            )
        );

        let f = (1 + x + x.pow(2)).pow(2).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![((1 + x + x.pow(2)).into_verbose(), Natural::from(2u8))]
                )
            )
        );

        //non-primitive cases
        let f = (2 + 2 * x).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![
                        (Polynomial::from_int(Integer::from(2)), Natural::from(1u8)),
                        ((1 + x).into_verbose(), Natural::from(1u8))
                    ]
                )
            )
        );

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![
                        (Polynomial::from_int(Integer::from(2)), Natural::from(2u8)),
                        (Polynomial::from_int(Integer::from(3)), Natural::from(1u8)),
                        ((2 + 3 * x).into_verbose(), Natural::from(1u8)),
                        ((x - 1).into_verbose(), Natural::from(2u8))
                    ]
                )
            )
        );

        let f = Polynomial::<Integer>::one();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(Polynomial::one(), vec![])
            )
        );

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_verbose();
        assert!(
            int_poly_fs.equal(
                &Integer::structure()
                    .polynomials()
                    .factorize_by_kroneckers_method(f, Integer::factor),
                &int_poly_fs.new_unit_and_powers_unchecked(
                    Polynomial::one(),
                    vec![
                        ((x.pow(4) + x + 1).into_verbose(), Natural::from(1u8)),
                        ((x.pow(3) + x + 1).into_verbose(), Natural::from(1u8))
                    ]
                )
            )
        );
    }

    #[test]
    fn test_fof_factor_over_rationals() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let f = (6 * (x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_verbose();
        let fs = f.factor();
        let rat_poly_fs = Rational::structure()
            .into_polynomials()
            .into_factorizations();

        println!("fs = {}", fs);

        assert!(rat_poly_fs.equal(
            &f.factor(),
            &rat_poly_fs.new_unit_and_powers_unchecked(
                Polynomial::constant(Rational::from(6)),
                vec![
                    ((x.pow(4) + x + 1).into_verbose(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).into_verbose(), Natural::from(1u8))
                ]
            )
        ));
    }
}
