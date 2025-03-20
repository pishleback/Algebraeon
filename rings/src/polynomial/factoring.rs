use super::super::structure::factorization::*;
use super::super::structure::structure::*;
use super::polynomial::*;
use algebraeon_sets::structure::*;

use algebraeon_nzq::natural::*;

impl<RS: UniqueFactorizationStructure + GreatestCommonDivisorStructure + CharZeroStructure>
    PolynomialStructure<RS>
where
    PolynomialStructure<RS>: Structure<Set = Polynomial<RS::Set>>
        + GreatestCommonDivisorStructure
        + UniqueFactorizationStructure,
{
    /// Reduce a factorization problem for polynomials over a ring of characteristic 0 to a factorization of primitive squarefree polynomials over the ring
    //https://en.wikipedia.org/wiki/Square-free_polynomial#Yun's_algorithm
    pub fn factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(
        &self,
        f: Polynomial<RS::Set>,
        primitive_sqfree_factorize: &impl Fn(Polynomial<RS::Set>) -> Factored<PolynomialStructure<RS>>,
    ) -> Factored<PolynomialStructure<RS>> {
        debug_assert!(!self.is_zero(&f));
        //look for a squarefree factorization of the form
        //  f = x * a_1^1 * a_2^2 * a_3^3 * ... * a_k^k
        //where each a_i is a squarefree primitive polynomial and x is an element of R

        let (content, prim) = self.factor_primitive(f.clone()).unwrap();

        let (content_unit, content_factors) = self
            .coeff_ring()
            .factor(&content)
            .unwrap()
            .unit_and_factors();

        let mut factors = Factored::new_unchecked(
            self.clone().into(),
            Polynomial::constant(content_unit),
            content_factors
                .into_iter()
                .map(|(factor, power)| (Polynomial::constant(factor), power))
                .collect(),
        );

        let mut f = prim;
        let f_prime = self.derivative(f.clone());
        let mut i: usize = 1;
        let a = self.gcd_by_primitive_subresultant(f.clone(), f_prime.clone());
        let mut b = self.div(&f, &a).unwrap();
        let mut c = self.div(&f_prime, &a).unwrap();
        let mut d = self.add(&self.neg(&self.derivative(b.clone())), &c);

        while self.degree(&f).unwrap() != 0 {
            let a = self.gcd(&b, &d);
            //a^i is a power of a squarefree factor of f
            factors.mul_mut(primitive_sqfree_factorize(a.clone()).pow(&Natural::from(i)));
            f = self.div(&f, &self.nat_pow(&a, &Natural::from(i))).unwrap();

            (b, c) = (self.div(&b, &a).unwrap(), self.div(&d, &a).unwrap());
            i += 1;
            d = self.add(&self.neg(&self.derivative(b.clone())), &c);
        }

        factors.mul_mut(Factored::factored_unit_unchecked(self.clone().into(), f));
        factors
    }
}

impl<RS: UniqueFactorizationStructure + GreatestCommonDivisorStructure + FiniteUnitsStructure>
    PolynomialStructure<RS>
where
    PolynomialStructure<RS>: Structure<Set = Polynomial<RS::Set>>,
{
    fn factor_primitive_linear_part(
        &self,
        mut f: Polynomial<RS::Set>,
    ) -> (Factored<PolynomialStructure<RS>>, Polynomial<RS::Set>)
    where
        PolynomialStructure<RS>: UniqueFactorizationStructure,
    {
        debug_assert!(self.is_primitive(f.clone()));

        /*
        If f is a primitive polynomial over R then any
        linear factor `(a+bx)` of `f(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n`
        is such that `a` divides `c0` and `b` divides `cn`
        so just factor c0 and cn and check all divisors
        */

        let mut linear_factors = Factored::factored_one(self.clone().into());
        'seek_linear_factor: while self.degree(&f).unwrap() > 0 {
            let c0 = self.coeff(&f, 0);
            if self.coeff_ring().is_zero(c0) {
                //linear factor of x
                f = self.div(&f, &self.var()).unwrap();
                linear_factors = Factored::mul(
                    linear_factors,
                    Factored::factored_irreducible_unchecked(self.clone().into(), self.var()),
                );
                continue 'seek_linear_factor;
            } else {
                //look for linear factors of the form (a+bx)
                let c0fs = self.coeff_ring().factor(self.coeff(&f, 0)).unwrap();
                let cnfs = self
                    .coeff_ring()
                    .factor(self.coeff(&f, self.degree(&f).unwrap()))
                    .unwrap();
                for a_assoc in c0fs.divisors() {
                    for u in self.coeff_ring().all_units() {
                        let a = self.coeff_ring().mul(&u, &a_assoc);
                        for b in cnfs.divisors() {
                            //a ranges over all divisors of c0
                            //b ranges over all divisors factors of cn up to associates
                            //try the linear factor (a+bx)
                            let lin = Polynomial::from_coeffs(vec![a.clone(), b]);
                            match self.div(&f, &lin) {
                                Ok(new_f) => {
                                    f = new_f;
                                    linear_factors = Factored::mul(
                                        linear_factors,
                                        Factored::factored_irreducible_unchecked(
                                            self.clone().into(),
                                            lin,
                                        ),
                                    );
                                    continue 'seek_linear_factor;
                                }
                                Err(RingDivisionError::NotDivisible) => {}
                                Err(RingDivisionError::DivideByZero) => panic!(),
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
    RS: UniqueFactorizationStructure
        + GreatestCommonDivisorStructure
        + CharZeroStructure
        + FiniteUnitsStructure
        + 'static,
> PolynomialStructure<RS>
where
    PolynomialStructure<RS>: Structure<Set = Polynomial<RS::Set>> + GreatestCommonDivisorStructure,
{
    fn find_factor_primitive_by_kroneckers_algorithm(
        &self,
        f: &Polynomial<RS::Set>,
    ) -> FindFactorResult<PolynomialStructure<RS>> {
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
                let y = self.evaluate(&f, &x);
                if !self.coeff_ring().is_zero(&y) {
                    f_points.push((x, self.coeff_ring().factor(&y).unwrap()));
                }
            }

            //compute all factors of each y value. choose the y with the most divisors to only factor up to units
            f_points.sort_by_cached_key(|(_x, yf)| yf.count_divisors());
            let _ = f_points.split_off(max_factor_degree + 1);
            //possible_g_points is (x, possible_y_values)
            let all_possible_g_points: Vec<(RS::Set, Vec<RS::Set>)> = f_points
                .into_iter()
                .rev()
                .enumerate()
                .map(|(i, (x, yf))| {
                    let mut y_divs = vec![];
                    for d in yf.divisors() {
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
                match self.interpolate_by_lagrange_basis(&possible_g_points) {
                    Some(g) => {
                        if self.degree(&g).unwrap() >= 1 {
                            //g is a possible proper divisor of f
                            match self.div(&f, &g) {
                                Ok(h) => {
                                    //g really is a proper divisor of f
                                    return FindFactorResult::Composite(g, h);
                                }
                                Err(RingDivisionError::NotDivisible) => {}
                                Err(RingDivisionError::DivideByZero) => panic!(),
                            }
                        }
                    }
                    None => {}
                }
            }
            //f is irreducible
            FindFactorResult::Irreducible
        }
    }

    pub fn factorize_by_kroneckers_method(&self, f: &Polynomial<RS::Set>) -> Option<Factored<Self>>
    where
        Self: UniqueFactorizationStructure,
    {
        if self.is_zero(f) {
            None
        } else {
            Some(
                self.factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(
                    f.clone(),
                    &|f| {
                        factorize_by_find_factor(self, f, &|f| {
                            self.find_factor_primitive_by_kroneckers_algorithm(&f)
                        })
                    },
                ),
            )
        }
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: UniqueFactorizationStructure
        + GreatestCommonDivisorStructure
        + CharZeroStructure
        + FiniteUnitsStructure
        + 'static,
    PolynomialStructure<R::Structure>: Structure<Set = Polynomial<R>>
        + GreatestCommonDivisorStructure
        + UniqueFactorizationStructure,
{
    pub fn factorize_by_kroneckers_method(
        &self,
    ) -> Option<Factored<PolynomialStructure<R::Structure>>> {
        Self::structure().factorize_by_kroneckers_method(self)
    }
}

impl<Fof: FieldOfFractionsStructure> PolynomialStructure<Fof>
where
    Self: Structure<Set = Polynomial<Fof::Set>> + UniqueFactorizationStructure,
    PolynomialStructure<Fof::RS>:
        Structure<Set = Polynomial<<Fof::RS as Structure>::Set>> + UniqueFactorizationStructure,
    Fof::RS: GreatestCommonDivisorStructure,
{
    pub fn factorize_by_factorize_primitive_part(
        &self,
        f: &Polynomial<Fof::Set>,
    ) -> Option<Factored<PolynomialStructure<Fof>>> {
        let (unit, prim) = self.factor_primitive_fof(f);
        let (prim_unit, prim_factors) =
            PolynomialStructure::new(self.coeff_ring().base_ring_structure())
                .factor(&prim)?
                .unit_and_factors();
        let mut fof_unit = prim_unit.apply_map(|c| self.coeff_ring().from_base_ring(c.clone()));
        let mut fof_factors = vec![];
        for (factor, power) in prim_factors.into_iter() {
            let fof_factor = factor.apply_map(|c| self.coeff_ring().from_base_ring(c.clone()));
            let (fof_factor_unit, fof_factor_prim) = self.factor_fav_assoc(&fof_factor);
            self.mul_mut(&mut fof_unit, &fof_factor_unit);
            fof_factors.push((fof_factor_prim, power));
        }
        let factors = Factored::new_unchecked(
            self.clone().into(),
            self.mul(&Polynomial::constant(unit), &fof_unit),
            fof_factors,
        );
        Some(factors)
    }
}

impl<RS: FieldStructure + FiniteUnitsStructure> PolynomialStructure<RS>
where
    Self: Structure<Set = Polynomial<RS::Set>>,
{
    fn find_factor_by_trying_all_factors(
        &self,
        f: Polynomial<RS::Set>,
    ) -> FindFactorResult<PolynomialStructure<RS>> {
        let f_deg = self.degree(&f).unwrap();
        let max_factor_degree = f_deg / 2;
        for d in 0..max_factor_degree {
            for mut coeffs in
                itertools::Itertools::multi_cartesian_product((0..d + 1).into_iter().map(|_d| {
                    let mut all_elems = vec![self.coeff_ring().zero()];
                    all_elems.append(&mut self.coeff_ring().all_units());
                    all_elems
                }))
            {
                coeffs.push(self.coeff_ring().one());
                let g = Polynomial::from_coeffs(coeffs);
                match self.div(&f, &g) {
                    Ok(h) => {
                        return FindFactorResult::Composite(g, h);
                    }
                    Err(RingDivisionError::NotDivisible) => {}
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
        }
        FindFactorResult::Irreducible
    }

    pub fn factorize_by_trying_all_factors(
        &self,
        f: Polynomial<RS::Set>,
    ) -> Option<Factored<PolynomialStructure<RS>>>
    where
        Self: UniqueFactorizationStructure,
    {
        if self.is_zero(&f) {
            None
        } else {
            Some(factorize_by_find_factor(self, f, &|f| {
                self.find_factor_by_trying_all_factors(f)
            }))
        }
    }
}

impl<F: MetaType> Polynomial<F>
where
    F::Structure: FieldStructure + FiniteUnitsStructure,
    PolynomialStructure<F::Structure>:
        Structure<Set = Polynomial<F>> + UniqueFactorizationStructure,
{
    pub fn factorize_by_trying_all_factors(
        &self,
    ) -> Option<Factored<PolynomialStructure<F::Structure>>> {
        Self::structure().factorize_by_trying_all_factors(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::elements::IntoErgonomic;

    use algebraeon_nzq::integer::*;
    use algebraeon_nzq::rational::*;

    use super::*;

    #[test]
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        //primitive cases
        let f = ((1 + x).pow(2)).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![((1 + x).into_verbose(), Natural::from(2u8))]
            )
        ));

        let f = (-1 - 2 * x).into_verbose();
        let fs1 = f.factorize_by_kroneckers_method().unwrap();
        let fs2 = &Factored::new_unchecked(
            Polynomial::<Integer>::structure().into(),
            Polynomial::neg(&Polynomial::one()),
            vec![((1 + 2 * x).into_verbose(), Natural::from(1u8))],
        );
        println!("fs1={} fs2={}", fs1, fs2);
        assert!(Factored::equal(&fs1, &fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    ((1 + x + x.pow(2)).into_verbose(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).into_verbose(), Natural::from(1u8))
                ]
            )
        ));

        let f = (1 + x + x.pow(2)).pow(2).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![((1 + x + x.pow(2)).into_verbose(), Natural::from(2u8))]
            )
        ));

        //non-primitive cases
        let f = (2 + 2 * x).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).into_verbose(), Natural::from(1u8))
                ]
            )
        ));

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(Integer::from(2)), Natural::from(2u8)),
                    (Polynomial::from_int(Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).into_verbose(), Natural::from(1u8)),
                    ((x - 1).into_verbose(), Natural::from(2u8))
                ]
            )
        ));

        let f = Polynomial::<Integer>::one();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![]
            )
        ));

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_verbose();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    ((x.pow(4) + x + 1).into_verbose(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).into_verbose(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_fof_factor_over_rationals() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let f = (6 * (x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_verbose();
        let fs = f.factor().unwrap();

        println!("fs = {}", fs);

        assert!(Factored::equal(
            &f.factor().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Rational>::structure().into(),
                Polynomial::constant(Rational::from(6)),
                vec![
                    ((x.pow(4) + x + 1).into_verbose(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).into_verbose(), Natural::from(1u8))
                ]
            )
        ));
    }
}
