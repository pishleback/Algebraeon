use itertools::Itertools;
use malachite_nz::natural::Natural;

use super::super::ring_structure::factorization::*;
use super::super::ring_structure::structure::*;
use super::polynomial::*;
use crate::linear::matrix::*;
use crate::linear::subspace::*;
use algebraeon_structure::*;

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// Reduce a factorization problem for polynomials over a finite field to factorizations of squarefree polynomials over the finite field
    // https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
    pub fn factorize_by_sqfree_factorize_over_finite_field(
        &self,
        f: Polynomial<FS::Set>,
        sqfree_factorize: &impl Fn(Polynomial<FS::Set>) -> Factored<PolynomialStructure<FS>>,
    ) -> Factored<PolynomialStructure<FS>>
    where
        PolynomialStructure<FS>: EuclideanDivisionStructure
            + GreatestCommonDivisorStructure
            + UniqueFactorizationStructure,
    {
        let (unit, f) = self.factor_fav_assoc(&f);
        // println!("f = {:?}", f);
        let mut factors = Factored::new_unchecked(self.clone().into(), unit, vec![]);
        //let w be the squarefree product of all factors of f of multiplicity not divisible by p
        let mut c = self.gcd(&f, &self.derivative(f.clone()));
        let mut w = self.div(&f, &c).unwrap();
        // println!("c = {:?}", c);
        // println!("w = {:?}", w);
        // println!("{:?}", self.divisible(&f, &c));
        // println!("{:?}", self.divisible(&f, &w));
        // println!("f' = {:?}", self.derivative(f.clone()));

        //find all the factors in w
        let mut i = Natural::from(1u8);
        while !self.equal(&w, &self.one()) {
            let y = self.gcd(&w, &c);
            factors.mul_mut(sqfree_factorize(self.div(&w, &y).unwrap()).pow(&i));
            let c_over_y = self.div(&c, &y).unwrap();
            (w, c) = (y, c_over_y);
            i += Natural::from(1u8);
        }

        // c is now the product, with multiplicity, of the remaining factors of f
        if !self.equal(&c, &self.one()) {
            //c = c^{1/p}
            let p = self.coeff_ring().characteristic_and_power().0;
            // println!("p = {:?} c = {:?}", p, c);
            let mut reduced_c_coeffs = vec![];
            for (k, coeff) in c.coeffs().into_iter().enumerate() {
                if Natural::from(k) % &p == 0 {
                    reduced_c_coeffs.push(coeff.clone());
                } else {
                    debug_assert!(self.coeff_ring().is_zero(coeff));
                }
            }
            let reduced_c = Polynomial::from_coeffs(reduced_c_coeffs);
            // println!("reduced_c = {}", reduced_c);
            factors.mul_mut(
                self.factorize_by_sqfree_factorize_over_finite_field(reduced_c, sqfree_factorize)
                    .pow(&p),
            );
        }
        factors
    }
}

impl<RS: UniqueFactorizationStructure + GreatestCommonDivisorStructure + CharZeroStructure>
    PolynomialStructure<RS>
where
    PolynomialStructure<RS>: Structure<Set = Polynomial<RS::Set>>
        + GreatestCommonDivisorStructure
        + UniqueFactorizationStructure,
{
    /// Reduce a factorization problem for polynomials over a ring of characteristic 0 to a factorization of primitive squarefree polynomials over the ring
    //https://en.wikipedia.org/wiki/Square-free_polynomial#Yun's_algorithm
    pub fn factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(
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
        let a = self.gcd(&f, &f_prime);
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
    ) -> Option<(Polynomial<RS::Set>, Polynomial<RS::Set>)> {
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
            None
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
                                    return Some((g, h));
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
            None
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
                self.factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(f.clone(), &|f| {
                    factorize_by_find_factor(self, f, &|f| {
                        self.find_factor_primitive_by_kroneckers_algorithm(&f)
                    })
                }),
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
    ) -> Option<(Polynomial<RS::Set>, Polynomial<RS::Set>)> {
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
                        return Some((g, h));
                    }
                    Err(RingDivisionError::NotDivisible) => {}
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
        }
        None
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

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    fn find_factor_by_berlekamps_algorithm(
        &self,
        f: Polynomial<FS::Set>,
    ) -> Option<(Polynomial<FS::Set>, Polynomial<FS::Set>)> {
        // println!("FACTOR {}", f);
        //f is squarefree
        let f_deg = self.degree(&f).unwrap();
        let all_elems = self.coeff_ring().all_elements();
        let q = all_elems.len();

        let row_polys = (0..f_deg)
            .map(|i| {
                self.rem(
                    &self.add(&self.var_pow(i * q), &self.neg(&self.var_pow(i))),
                    &f,
                )
            })
            .collect_vec();
        let mat = Matrix::construct(f_deg, f_deg, |i, j| self.coeff(&row_polys[i], j).clone());
        // mat.pprint();
        //the column kernel gives a basis of berlekamp subalgebra - all polynomials g such that g^q=g
        let mat_struct = MatrixStructure::<FS>::new(self.coeff_ring().into());
        let linlat_struct = LinearLatticeStructure::<FS>::new(self.coeff_ring().into());
        let ker = mat_struct.row_kernel(mat);
        // ker.pprint();
        let ker_rank = linlat_struct.rank(&ker);
        let ker_basis = linlat_struct
            .basis_matrices(&ker)
            .into_iter()
            .map(|col| {
                Polynomial::from_coeffs((0..f_deg).map(|c| col.at(0, c).unwrap().clone()).collect())
            })
            .collect_vec();

        for berlekamp_subspace_coeffs in (0..ker_rank)
            .map(|_i| all_elems.clone().into_iter())
            .multi_cartesian_product()
        {
            let h = self.sum(
                (0..ker_rank)
                    .map(|i| {
                        self.mul(
                            &Polynomial::constant(berlekamp_subspace_coeffs[i].clone()),
                            &ker_basis[i],
                        )
                    })
                    .collect(),
            );
            //g is a possible non-trivial factor
            let g = self.gcd(&h, &f);
            // println!("g = {}", g);
            let g_deg = self.degree(&g).unwrap();
            if g_deg != 0 && g_deg != f_deg {
                match self.div(&f, &g) {
                    Ok(g_prime) => {
                        return Some((g, g_prime));
                    }
                    Err(_) => panic!(),
                }
            }
        }
        None
    }

    pub fn factorize_by_berlekamps_algorithm(
        &self,
        f: Polynomial<FS::Set>,
    ) -> Option<Factored<PolynomialStructure<FS>>>
    where
        PolynomialStructure<FS>: UniqueFactorizationStructure,
    {
        if self.is_zero(&f) {
            None
        } else {
            Some(
                self.factorize_by_sqfree_factorize_over_finite_field(f, &|f| {
                    factorize_by_find_factor(self, f, &|f| {
                        self.find_factor_by_berlekamps_algorithm(f)
                    })
                }),
            )
        }
    }
}

impl<F: MetaType> Polynomial<F>
where
    F::Structure: FiniteFieldStructure,
    PolynomialStructure<F::Structure>:
        Structure<Set = Polynomial<F>> + UniqueFactorizationStructure,
{
    pub fn factorize_by_berlekamps_algorithm(
        &self,
    ) -> Option<Factored<PolynomialStructure<F::Structure>>> {
        Self::structure().factorize_by_berlekamps_algorithm(self.clone())
    }
}

#[cfg(test)]
mod tests {

    use crate::elements::*;
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use crate::number::finite_fields::{modulo::Modulo, quaternary_field::QuaternaryField};

    use super::*;

    #[test]
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Polynomial::<Integer>::var().into_ring();

        //primitive cases
        let f = ((1 + x).pow(2)).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![((1 + x).into_set(), Natural::from(2u8))]
            )
        ));

        let f = (-1 - 2 * x).into_set();
        let fs1 = f.factorize_by_kroneckers_method().unwrap();
        let fs2 = &Factored::new_unchecked(
            Polynomial::<Integer>::structure().into(),
            Polynomial::neg(&Polynomial::one()),
            vec![((1 + 2 * x).into_set(), Natural::from(1u8))],
        );
        println!("fs1={} fs2={}", fs1, fs2);
        assert!(Factored::equal(&fs1, &fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    ((1 + x + x.pow(2)).into_set(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).into_set(), Natural::from(1u8))
                ]
            )
        ));

        let f = (1 + x + x.pow(2)).pow(2).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![((1 + x + x.pow(2)).into_set(), Natural::from(2u8))]
            )
        ));

        //non-primitive cases
        let f = (2 + 2 * x).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).into_set(), Natural::from(1u8))
                ]
            )
        ));

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(2u8)),
                    (Polynomial::from_int(&Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).into_set(), Natural::from(1u8)),
                    ((x - 1).into_set(), Natural::from(2u8))
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

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_set();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Integer>::structure().into(),
                Polynomial::one(),
                vec![
                    ((x.pow(4) + x + 1).into_set(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).into_set(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_fof_factor_over_rationals() {
        let x = &Polynomial::<Rational>::var().into_ring();
        let f = (6 * (x.pow(4) + x + 1) * (x.pow(3) + x + 1)).into_set();

        assert!(Factored::equal(
            &f.factor().unwrap(),
            &Factored::new_unchecked(
                Polynomial::<Rational>::structure().into(),
                Polynomial::constant(Rational::from(6)),
                vec![
                    ((x.pow(4) + x + 1).into_set(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).into_set(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_factorize_over_f2_example1() {
        let x = &Polynomial::<Modulo<2>>::var().into_ring();
        let p = (1 + x.pow(4) + x.pow(5)).pow(12).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<2>>::structure().into(),
            Polynomial::one(),
            vec![
                ((1 + x + x.pow(2)).into_set(), Natural::from(12u32)),
                ((1 + x + x.pow(3)).into_set(), Natural::from(12u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f2_example2() {
        let x = &Polynomial::<Modulo<2>>::var().into_ring();
        let p = ((1 + x.pow(4) + x.pow(7)).pow(6) * (1 + x.pow(6) + x.pow(7)).pow(4)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<2>>::structure().into(),
            Polynomial::one(),
            vec![
                ((1 + x.pow(4) + x.pow(7)).into_set(), Natural::from(6u32)),
                ((1 + x.pow(6) + x.pow(7)).into_set(), Natural::from(4u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example1() {
        let x = &Polynomial::<Modulo<5>>::var().into_ring();
        let p = (1 + x.pow(4)).pow(5).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<5>>::structure().into(),
            Polynomial::one(),
            vec![
                ((2 + x.pow(2)).into_set(), Natural::from(5u8)),
                ((3 + x.pow(2)).into_set(), Natural::from(5u8)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example2() {
        let x = &Polynomial::<Modulo<5>>::var().into_ring();
        let p = (3 + 2 * x.pow(2) + x.pow(4) + x.pow(6)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<5>>::structure().into(),
            Polynomial::one(),
            vec![((2 + x.pow(2)).into_set(), Natural::from(3u8))],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f31_example1() {
        let x = &Polynomial::<Modulo<31>>::var().into_ring();
        let p = (1 + x.pow(27) + 8 * x.pow(30)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<31>>::structure().into(),
            Polynomial::constant(Modulo::from_int(&Integer::from(8))),
            vec![
                ((12 + x.pow(3)).into_set(), Natural::from(1u32)),
                (
                    (25 + 27 * x + 3 * x.pow(2) + 3 * x.pow(3) + 29 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (1 + 24 * x + 3 * x.pow(2) + 15 * x.pow(3) + 12 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (21 + 12 * x.pow(3) + 22 * x.pow(6) + 4 * x.pow(9) + x.pow(12))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (5 + 11 * x + 3 * x.pow(2) + 13 * x.pow(3) + 21 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
            ],
        );

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f4_example1() {
        let x = &Polynomial::<QuaternaryField>::var().into_ring();

        let a = x - Polynomial::constant(QuaternaryField::One).into_ring();
        let b = x - Polynomial::constant(QuaternaryField::Alpha).into_ring();
        let c = x - Polynomial::constant(QuaternaryField::Beta).into_ring();
        let p = (1 - x.pow(3)).pow(48).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<QuaternaryField>::structure().into(),
            Polynomial::one(),
            vec![
                (a.into_set(), Natural::from(48u32)),
                (b.into_set(), Natural::from(48u32)),
                (c.into_set(), Natural::from(48u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }
}
