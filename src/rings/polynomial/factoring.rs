use std::thread::sleep;
use std::time::Duration;

use itertools::Itertools;
use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use crate::rings::linear::matrix::Matrix;

use super::super::numbers::nzq::*;
use super::super::polynomial::poly::*;
use super::super::ring::*;

//https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
pub fn factorize_by_sqfree_factorize_over_finite_field<FS: FiniteFieldStructure>(
    fs: &FS,
    mut f: Polynomial<FS::ElementT>,
    sqfree_factorize: &impl Fn(Polynomial<FS::ElementT>) -> Factored<Polynomial<FS::ElementT>>,
) -> Factored<Polynomial<FS::ElementT>>
where
    Polynomial<FS::ElementT>:
        EuclideanDomain + GreatestCommonDivisorDomain + UniqueFactorizationDomain,
{
    // println!("f = {:?}", f);
    let (unit, f) = f.factor_fav_assoc();
    let mut factors = Factored::new_unchecked(unit, vec![]);
    //let w be the squarefree product of all factors of f of multiplicity not divisible by p
    let mut c = Polynomial::gcd(f.clone(), f.clone().derivative());
    let mut w = Polynomial::div_refs(&f, &c).unwrap();
    // println!("c = {:?}", c);
    // println!("w = {:?}", w);

    //find all the factors in w
    let mut i = Natural::from(1u8);
    while w != Polynomial::one() {
        let y = Polynomial::gcd(w.clone(), c.clone());
        factors.mul_mut(sqfree_factorize(Polynomial::div_rref(w, &y).unwrap()).pow(&i));
        let c_over_y = Polynomial::div_rref(c, &y).unwrap();
        (w, c) = (y, c_over_y);
        i += Natural::from(1u8);
    }

    // c is now the product, with multiplicity, of the remaining factors of f
    if c != Polynomial::one() {
        // println!("c = {}", c);
        //c = c^{1/p}
        let p = fs.characteristic_and_power().0;
        let mut reduced_c_coeffs = vec![];
        for (k, coeff) in c.coeffs().into_iter().enumerate() {
            if Natural::from(k) % &p == 0 {
                reduced_c_coeffs.push(coeff);
            } else {
                debug_assert_eq!(coeff, FS::ElementT::zero());
            }
        }
        let reduced_c = Polynomial::from_coeffs(reduced_c_coeffs);
        // println!("reduced_c = {}", reduced_c);
        factors.mul_mut(factorize_by_sqfree_factorize_over_finite_field(fs, reduced_c, sqfree_factorize).pow(&p));
    }
    factors
}

//https://en.wikipedia.org/wiki/Square-free_polynomial#Yun's_algorithm
pub fn factorize_by_primitive_sqfree_factorize_by_yuns_algorithm<
    R: UniqueFactorizationDomain + GreatestCommonDivisorDomain + CharacteristicZero,
>(
    mut f: Polynomial<R>,
    primitive_sqfree_factorize: &impl Fn(Polynomial<R>) -> Factored<Polynomial<R>>,
) -> Factored<Polynomial<R>>
where
    Polynomial<R>: GreatestCommonDivisorDomain + UniqueFactorizationDomain,
{
    debug_assert_ne!(f, Polynomial::zero());
    //look for a squarefree factorization of the form
    //  f = x * a_1^1 * a_2^2 * a_3^3 * ... * a_k^k
    //where each a_i is a squarefree primitive polynomial and x is an element of R

    let (content, prim) = f.clone().factor_primitive().unwrap();

    let (content_unit, content_factors) = content.factor().unwrap().unit_and_factors();

    let mut factors = Factored::new_unchecked(
        Polynomial::constant(content_unit),
        content_factors
            .into_iter()
            .map(|(factor, power)| (Polynomial::constant(factor), power))
            .collect(),
    );

    let f_prime = f.clone().derivative();
    let mut i: usize = 1;
    let mut a = Polynomial::gcd(f.clone(), f_prime.clone());
    let mut b = Polynomial::div_refs(&f, &a).unwrap();
    let mut c = Polynomial::div_refs(&f_prime, &a).unwrap();
    let mut d = Polynomial::add_ref(b.clone().derivative().neg(), &c);
    while b != Polynomial::one() {
        let a = Polynomial::gcd(b.clone(), d.clone());
        //a^i is a power of a squarefree factor of f
        factors.mul_mut(primitive_sqfree_factorize(a.clone()).pow(&Natural::from(i)));
        (b, c) = (
            Polynomial::div_refs(&b, &a).unwrap(),
            Polynomial::div_refs(&d, &a).unwrap(),
        );
        i += 1;
        d = Polynomial::add_ref(b.clone().derivative().neg(), &c);
    }
    factors
}

impl<Ring: UniqueFactorizationDomain + GreatestCommonDivisorDomain + FiniteUnits> Polynomial<Ring> {
    fn factor_primitive_linear_part(
        &self,
        mut f: Polynomial<Ring>,
    ) -> (Factored<Polynomial<Ring>>, Polynomial<Ring>)
    where
        Polynomial<Ring>: UniqueFactorizationDomain,
    {
        //f should be a primitive polynomial over R
        //linear factor (a+bx) of c0 + c1*x + c2*x^2 + ... + cn*x^n
        //must be such that a divides c0 and b divides cn
        //so just factor c0 and cn and check all divisors

        let mut linear_factors = Factored::factored_one();
        'seek_linear_factor: while f.degree().unwrap() > 0 {
            let c0 = Self::coeff(&f, 0);
            if c0 == Ring::zero() {
                //linear factor of x
                f = Self::div(f, Self::var()).unwrap();
                linear_factors = Factored::mul(
                    linear_factors,
                    Factored::factored_irreducible_unchecked(Self::var()),
                );
                continue 'seek_linear_factor;
            } else {
                //look for linear factors of the form (a+bx)
                let c0fs = Ring::factor(&Self::coeff(&f, 0)).unwrap();
                let cnfs = Ring::factor(&f.coeff(f.degree().unwrap())).unwrap();
                for a_assoc in c0fs.divisors() {
                    for u in Ring::all_units() {
                        let a = Ring::mul_ref(u, &a_assoc);
                        for b in cnfs.divisors() {
                            //a ranges over all divisors of c0
                            //b ranges over all divisors factors of cn up to associates
                            //try the linear factor (a+bx)
                            let lin = Self::from_coeffs(vec![a.clone(), b]);
                            match Self::div_refs(&f, &lin) {
                                Ok(new_f) => {
                                    f = new_f;
                                    linear_factors = Factored::mul(
                                        linear_factors,
                                        Factored::factored_irreducible_unchecked(lin),
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
        Ring: UniqueFactorizationDomain + GreatestCommonDivisorDomain + CharacteristicZero + FiniteUnits,
    > Polynomial<Ring>
where
    Polynomial<Ring>: GreatestCommonDivisorDomain,
{
    fn find_factor_primitive_by_kroneckers_algorithm(
        self,
    ) -> Option<(Polynomial<Ring>, Polynomial<Ring>)> {
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
        let f = self;
        let f_deg = f.degree().unwrap();
        if f_deg == 1 {
            //linear factor is irreducible
            None
        } else {
            let max_factor_degree = f_deg / 2;
            let mut f_points = vec![];
            let mut elem_gen = Ring::generate_distinct_elements();
            //take more samples than necessary, then take the subset with the smallest number of divisors
            while f_points.len() < 3 * (max_factor_degree + 1) {
                //loop terminates because polynomial over integral domain has finitely many roots
                let x = elem_gen.next().unwrap();
                let y = Polynomial::evaluate(&f, &x);
                if y != Ring::zero() {
                    f_points.push((x, Ring::factor(&y).unwrap()));
                }
            }

            //compute all factors of each y value. choose the y with the most divisors to only factor up to units
            f_points.sort_by_cached_key(|(_x, yf)| yf.count_divisors());
            let _ = f_points.split_off(max_factor_degree + 1);
            //possible_g_points is (x, possible_y_values)
            let all_possible_g_points: Vec<(Ring, Vec<Ring>)> = f_points
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
                            for u in Ring::all_units() {
                                y_divs.push(Ring::mul_ref(u, &d));
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
                match Polynomial::interpolate_by_lagrange_basis(&possible_g_points) {
                    Some(g) => {
                        if g.degree().unwrap() >= 1 {
                            //g is a possible proper divisor of f
                            match Polynomial::div_refs(&f, &g) {
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

    pub fn factorize_by_kroneckers_method(self) -> Option<Factored<Self>>
    where
        Self: UniqueFactorizationDomain,
    {
        if self == Self::zero() {
            None
        } else {
            Some(factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(
                self,
                &|f| {
                    factorize_by_find_factor(f, &|f| {
                        f.find_factor_primitive_by_kroneckers_algorithm()
                    })
                },
            ))
        }
    }
}

impl Polynomial<Integer> {
    fn find_factor_primitive_sqfree_by_zassenhaus_algorithm(
        self,
    ) -> Option<(Polynomial<Integer>, Polynomial<Integer>)> {
        let f = self;
        let f_deg = f.degree().unwrap();
        debug_assert_ne!(f_deg, 0);
        println!("zassenhaus: {}", f);
        if f_deg == 1 {
            None
        } else {
            let prime_gen = NaturalPrimeGenerator::new();
            for p in prime_gen.take(20) {
                println!("{:?}", p);
                let f_mod_p = f.apply_map_ref(|c| {
                    UniversalEuclideanQuotient::<true, _>::new(c.clone(), Integer::from(&p))
                });

                println!("{}", f_mod_p);

                sleep(Duration::from_millis(100));
            }

            todo!()
        }
    }

    pub fn factorize_by_zassenhaus_algorithm(self) -> Option<Factored<Self>> {
        if self == Self::zero() {
            None
        } else {
            Some(factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(
                self,
                &|f| {
                    factorize_by_find_factor(f, &|f| {
                        f.find_factor_primitive_sqfree_by_zassenhaus_algorithm()
                    })
                },
            ))
        }
    }
}

impl<Fof: FieldOfFractions> Polynomial<Fof>
where
    Self: UniqueFactorizationDomain,
    Polynomial<Fof::R>: UniqueFactorizationDomain,
    Fof::R: GreatestCommonDivisorDomain,
{
    pub fn factorize_by_factorize_primitive_part(&self) -> Option<Factored<Polynomial<Fof>>> {
        let (unit, prim) = self.factor_primitive_fof();
        let (prim_unit, prim_factors) = prim.factor()?.unit_and_factors();
        let mut fof_unit = prim_unit.apply_map(|c| Fof::from_base_ring(c));
        let mut fof_factors = vec![];
        for (factor, power) in prim_factors.into_iter() {
            let fof_factor = factor.apply_map(Fof::from_base_ring);
            let (fof_factor_unit, fof_factor_prim) = fof_factor.factor_fav_assoc();
            fof_unit.mul_mut(&fof_factor_unit);
            fof_factors.push((fof_factor_prim, power));
        }
        let factors = Factored::new_unchecked(
            Polynomial::mul(Polynomial::constant(unit), fof_unit),
            fof_factors,
        );
        Some(factors)
    }
}

impl<Ring: Field + FiniteUnits> Polynomial<Ring> {
    pub fn factorize_by_trying_all_factors(self) -> Option<Factored<Polynomial<Ring>>>
    where
        Self: UniqueFactorizationDomain,
    {
        let f = self;
        if f == Self::zero() {
            None
        } else {
            fn partial_factor<Ring: Field + FiniteUnits>(
                f: Polynomial<Ring>,
            ) -> Option<(Polynomial<Ring>, Polynomial<Ring>)> {
                let f_deg = f.degree().unwrap();
                let max_factor_degree = f_deg / 2;
                for d in 0..max_factor_degree {
                    for mut coeffs in itertools::Itertools::multi_cartesian_product(
                        (0..d + 1).into_iter().map(|_d| {
                            let mut all_elems = vec![Ring::zero()];
                            all_elems.append(&mut Ring::all_units());
                            all_elems
                        }),
                    ) {
                        coeffs.push(Ring::one());
                        let g = Polynomial::from_coeffs(coeffs);
                        match Polynomial::div_refs(&f, &g) {
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
            Some(factorize_by_find_factor(f, &partial_factor::<Ring>))
        }
    }
}

pub fn factorize_by_berlekamps_algorithm<FS: FiniteFieldStructure>(
    fs: &FS,
    f: Polynomial<FS::ElementT>,
) -> Option<Factored<Polynomial<FS::ElementT>>>
where
    Polynomial<FS::ElementT>: UniqueFactorizationDomain,
{
    if f == Polynomial::zero() {
        None
    } else {
        fn partial_factor<FS: FiniteFieldStructure>(
            fs: &FS,
            f: Polynomial<FS::ElementT>,
        ) -> Option<(Polynomial<FS::ElementT>, Polynomial<FS::ElementT>)> {
            // println!("FACTOR {}", f);
            //f is squarefree
            let f_deg = f.degree().unwrap();
            let all_elems = fs.all_elements();
            let q = all_elems.len();

            let mut row_polys = (0..f_deg)
                .map(|i| {
                    Polynomial::rem_rref(
                        Polynomial::add(Polynomial::var_pow(i * q), Polynomial::var_pow(i).neg()),
                        &f,
                    )
                })
                .collect_vec();
            let mat = Matrix::construct(f_deg, f_deg, |i, j| row_polys[i].coeff(j));
            // mat.pprint();
            //the column kernel gives a basis of berlekamp subalgebra - all polynomials g such that g^q=g
            let ker = mat.row_kernel();
            // ker.pprint();
            let ker_rank = ker.rank();
            let ker_basis = ker
                .basis_matrices()
                .into_iter()
                .map(|col| {
                    Polynomial::from_coeffs(
                        (0..f_deg).map(|c| col.at(0, c).unwrap().clone()).collect(),
                    )
                })
                .collect_vec();

            for berlekamp_subspace_coeffs in (0..ker_rank)
                .map(|i| all_elems.clone().into_iter())
                .multi_cartesian_product()
            {
                let h = Polynomial::sum(
                    (0..ker_rank)
                        .map(|i| {
                            Polynomial::mul_ref(
                                Polynomial::constant(berlekamp_subspace_coeffs[i].clone()),
                                &ker_basis[i],
                            )
                        })
                        .collect(),
                );
                //g is a possible non-trivial factor
                let g = Polynomial::gcd(h.clone(), f.clone());
                // println!("g = {}", g);
                let g_deg = g.degree().unwrap();
                if g_deg != 0 && g_deg != f_deg {
                    match Polynomial::div_rref(f, &g) {
                        Ok(g_prime) => {
                            return Some((g, g_prime));
                        }
                        Err(_) => panic!(),
                    }
                }
            }
            None
        }
        Some(factorize_by_sqfree_factorize_over_finite_field(fs, f, &|f| {
            factorize_by_find_factor(f, &|f| partial_factor(fs, f))
        }))
    }
}

impl<F: FiniteField> Polynomial<F> {
    pub fn factorize_by_berlekamps_algorithm(self: Polynomial<F>) -> Option<Factored<Polynomial<F>>>
    where
        Polynomial<F>: UniqueFactorizationDomain,
    {
        factorize_by_berlekamps_algorithm(&CannonicalFiniteFieldStructure::new(), self)
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use crate::rings::numbers::small_fields::QuaternaryField;

    use super::super::super::ergonomic::*;
    use super::*;

    #[test]
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());

        //primitive cases
        let f = ((1 + x).pow(2)).elem();
        assert!(Factored::equal(
            &Polynomial::factorize_by_kroneckers_method(f).unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![((1 + x).elem(), Natural::from(2u8))]
            )
        ));

        let f = (-1 - 2 * x).elem();
        let fs1 = f.factorize_by_kroneckers_method().unwrap();
        let fs2 = &Factored::new_unchecked(
            Polynomial::neg(Polynomial::one()),
            vec![((1 + 2 * x).elem(), Natural::from(1u8))],
        );
        println!("fs1={:?} fs2={:?}", fs1, fs2);
        assert!(Factored::equal(&fs1, &fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    ((1 + x + x.pow(2)).elem(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (1 + x + x.pow(2)).pow(2).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![((1 + x + x.pow(2)).elem(), Natural::from(2u8))]
            )
        ));

        //non-primitive cases
        let f = (2 + 2 * x).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(2u8)),
                    (Polynomial::from_int(&Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).elem(), Natural::from(1u8)),
                    ((x - 1).elem(), Natural::from(2u8))
                ]
            )
        ));

        let f = Polynomial::<Integer>::one();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(Polynomial::one(), vec![])
        ));

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    ((x.pow(4) + x + 1).elem(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).elem(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_factorize_over_f2_example1() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<Modulo<2>>::var());
        let p = (1 + x.pow(4) + x.pow(5)).pow(12).elem();
        let ans = Factored::new_unchecked(
            Polynomial::one(),
            vec![
                ((1 + x + x.pow(2)).elem(), Natural::from(12u32)),
                ((1 + x + x.pow(3)).elem(), Natural::from(12u32)),
            ],
        );

        let f = p.clone().factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f2_example2() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<Modulo<2>>::var());
        let p = ((1 + x.pow(4) + x.pow(7)).pow(6) * (1 + x.pow(6) + x.pow(7)).pow(4)).elem();
        let ans = Factored::new_unchecked(
            Polynomial::one(),
            vec![
                ((1 + x.pow(4) + x.pow(7)).elem(), Natural::from(6u32)),
                ((1 + x.pow(6) + x.pow(7)).elem(), Natural::from(4u32)),
            ],
        );

        let f = p.clone().factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example1() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<Modulo<5>>::var());
        let p = (1 + x.pow(4)).pow(5).elem();
        let ans = Factored::new_unchecked(
            Polynomial::one(),
            vec![
                ((2 + x.pow(2)).elem(), Natural::from(5u8)),
                ((3 + x.pow(2)).elem(), Natural::from(5u8)),
            ],
        );

        let f = p.clone().factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example2() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<Modulo<5>>::var());
        let p = (3 + 2 * x.pow(2) + x.pow(4) + x.pow(6)).elem();
        let ans = Factored::new_unchecked(
            Polynomial::one(),
            vec![((2 + x.pow(2)).elem(), Natural::from(3u8))],
        );

        let f = p.clone().factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f31_example1() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<Modulo<31>>::var());
        let p = (1 + x.pow(27) + 8 * x.pow(30)).elem();
        let ans = Factored::new_unchecked(
            Polynomial::constant(Modulo::from_int(&Integer::from(8))),
            vec![
                ((12 + x.pow(3)).elem(), Natural::from(1u32)),
                (
                    (25 + 27 * x + 3 * x.pow(2) + 3 * x.pow(3) + 29 * x.pow(4) + x.pow(5)).elem(),
                    Natural::from(1u32),
                ),
                (
                    (1 + 24 * x + 3 * x.pow(2) + 15 * x.pow(3) + 12 * x.pow(4) + x.pow(5)).elem(),
                    Natural::from(1u32),
                ),
                (
                    (21 + 12 * x.pow(3) + 22 * x.pow(6) + 4 * x.pow(9) + x.pow(12)).elem(),
                    Natural::from(1u32),
                ),
                (
                    (5 + 11 * x + 3 * x.pow(2) + 13 * x.pow(3) + 21 * x.pow(4) + x.pow(5)).elem(),
                    Natural::from(1u32),
                ),
            ],
        );

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f4_example1() {
        use super::super::super::numbers::small_modulo::*;

        let x = &Ergonomic::new(Polynomial::<QuaternaryField>::var());

        let a = x - Ergonomic::new(Polynomial::constant(QuaternaryField::One));
        let b = x - Ergonomic::new(Polynomial::constant(QuaternaryField::Alpha));
        let c = x - Ergonomic::new(Polynomial::constant(QuaternaryField::Beta));
        let p = (1 - x.pow(3)).pow(48).elem();
        let ans = Factored::new_unchecked(
            Polynomial::one(),
            vec![
                (a.elem(), Natural::from(48u32)),
                (b.elem(), Natural::from(48u32)),
                (c.elem(), Natural::from(48u32)),
            ],
        );

        let f = p.clone().factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.clone().factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }
}
