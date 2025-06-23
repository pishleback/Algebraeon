use super::polynomial_quotient_number_field::AlgebraicNumberFieldPolynomialQuotientStructure;
use crate::{matrix::Matrix, polynomial::*, structure::*};
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

//convert between nested polynomials with bounded degree and rational vectors
//e.g. (a + bx) + (c + dx)x <-> (a, b, c, d)
fn double_poly_to_row(
    outer_poly_len: usize,
    inner_poly_len: usize,
    a: Polynomial<Polynomial<Rational>>,
) -> Matrix<Rational> {
    // let n = outer_poly_len * inner_poly_len;
    let rat_poly_poly = Rational::structure()
        .into_polynomial_ring()
        .into_polynomial_ring();

    debug_assert!(rat_poly_poly.num_coeffs(&a) <= outer_poly_len);
    for c in a.coeffs() {
        debug_assert!(c.num_coeffs() <= inner_poly_len);
    }
    Matrix::from_rows(vec![
        (0..outer_poly_len)
            .map(|i| rat_poly_poly.coeff(&a, i).into_owned())
            .flat_map(|c| {
                (0..inner_poly_len)
                    .map(|j| c.coeff(j).into_owned())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>(),
    ])
}

fn row_to_double_poly(
    outer_poly_len: usize,
    inner_poly_len: usize,
    a: Matrix<Rational>,
) -> Polynomial<Polynomial<Rational>> {
    let n = outer_poly_len * inner_poly_len;
    debug_assert_eq!(a.cols(), n);
    debug_assert_eq!(a.rows(), 1);
    let mut a_vec = (0..a.cols())
        .rev()
        .map(|c| a.at(0, c).unwrap())
        .collect::<Vec<_>>();
    debug_assert_eq!(a_vec.len(), n);
    let mut coeffs = vec![];
    for _i in 0..outer_poly_len {
        let mut coeff_coeffs = vec![];
        for _j in 0..inner_poly_len {
            coeff_coeffs.push(a_vec.pop().unwrap().clone());
        }
        coeffs.push(Polynomial::from_coeffs(coeff_coeffs));
    }
    Polynomial::from_coeffs(coeffs)
}

impl<B: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>>
    PolynomialStructure<AlgebraicNumberFieldPolynomialQuotientStructure, B>
{
    /*
        input:  A polynomial f(x) over an algebraic number field K, return
        output: The polynomial \prod_{i=1}^n \sigma_i(f) over \mathbb{Q}
                where \sigma_1,... ,\sigma_n are the n=deg(K) complex embeddings of K
                and \sigma_i(f) is the application of \sigma_i to the coefficients of f

        example:
                K = Q[θ] where θ has minimal polynomial x + yθ + zθ² = 0
                Let σ₀ σ₁ σ₂ be the 3 complex embeddings of K into the complex numbers
                f(x) = a + bθ + cθ²
                a = a(θ) for some degree 3 polynomial a
                b = b(θ) for some degree 3 polynomial b
                c = c(θ) for some degree 3 polynomial c
                N(f(x)) = (σ₀(a) + σ₀(b)x + σ₀(c)x²) (σ₁(a) + σ₁(b)x + σ₁(c)x²) (σ₂(a) + σ₂(b)x + σ₂(c)x²)

                        = ( σ₀(a)σ₁(a)σ₂(a) )
                        + ( 2σ₀(a)σ₁(a)σ₂(b) + 2σ₀(a)σ₁(b)σ₂(a) + 2σ₀(b)σ₁(a)σ₂(a) ) x
                        + ( 2σ₀(a)σ₁(a)σ₂(c) + 2σ₀(a)σ₁(c)σ₂(a) + 2σ₀(c)σ₁(a)σ₂(a)
                          + 2σ₀(a)σ₁(b)σ₂(b) + 2σ₀(b)σ₁(a)σ₂(b) + 2σ₀(b)σ₁(b)σ₂(a) ) x²
                        + ...

                        = ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) )
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x²
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x³
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x⁴
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x⁵
                        + ( symmetric polynomial in σ₀(θ) σ₁(θ) σ₂(θ) ) x⁶
                The elementary symmetric polynomials in σ₀(θ) σ₁(θ) σ₂(θ) are (up to sign flips) the (rational) coefficients of the minimal polynomial of θ
    */
    pub fn polynomial_norm(&self, f: &Polynomial<Polynomial<Rational>>) -> Polynomial<Rational> {
        // println!("f = {}", f);
        // panic!();

        // println!("start");

        let n = self.coeff_ring().degree();
        // println!("f = {:?}", f);

        let embedding_vars = (0..n)
            .map(|i| Variable::new(format!("σ{}", ss_num(i))))
            .collect::<Vec<_>>();
        // println!("embeddings = {:?}", embedding_vars);

        let rational_poly_multipoly_structure = Rational::structure()
            .into_multivariable_polynomial_ring()
            .into_polynomial_ring();

        let norm_f_sym = rational_poly_multipoly_structure.product(
            embedding_vars
                .iter()
                .map(|sigma| {
                    f.apply_map(|c| {
                        rational_poly_multipoly_structure.evaluate(
                            &c.apply_map(|x| MultiPolynomial::constant(x.clone())),
                            &MultiPolynomial::var(sigma.clone()),
                        )
                    })
                })
                .collect(),
        );

        // for (i, coeff) in norm_f_sym.coeffs().iter().enumerate() {
        //     println!("norm_f_sym {} = {}", i, coeff);
        // }

        let e_vals = {
            let mut min_poly_coeffs = self
                .coeff_ring()
                .modulus()
                .coeffs()
                .into_iter()
                .cloned()
                .collect::<Vec<_>>();

            let lc = min_poly_coeffs.pop().unwrap();
            debug_assert_eq!(lc, Rational::ONE);

            min_poly_coeffs
                .into_iter()
                .rev()
                .enumerate()
                .map(|(i, c)| if i % 2 == 0 { -c } else { c })
                .collect::<Vec<_>>()
        };
        debug_assert_eq!(e_vals.len(), n);
        debug_assert_eq!(embedding_vars.len(), n);

        // println!("e_vals = {:?}", e_vals);

        let norm_f: Polynomial<Rational> = norm_f_sym.apply_map(|sym_sigma| {
            let (elem_sym_vars, sym_elem) = sym_sigma
                .as_elementary_symmetric_polynomials_unchecked(embedding_vars.iter().collect());

            // println!("sym_elem = {}", sym_elem);

            debug_assert_eq!(elem_sym_vars.len(), n);
            sym_elem.evaluate(
                (0..n)
                    .map(|i| (elem_sym_vars[i].clone(), &e_vals[i]))
                    .collect(),
            )
        });

        // println!("end");

        norm_f
    }

    pub fn factor_primitive_sqfree_by_symmetric_root_polynomials(
        &self,
        p: &<Self as SetSignature>::Set,
    ) -> crate::structure::FactoredRingElement<<Self as SetSignature>::Set> {
        //https://www.cse.iitk.ac.in/users/nitin/courses/scribed2-WS2011-12.pdf

        let rat_poly_poly_poly = Rational::structure()
            .into_polynomial_ring()
            .into_polynomial_ring()
            .into_polynomial_ring();

        debug_assert!(!self.is_zero(p));
        //Let K = Q[θ] be the number field over which we are factoring
        let anf = self.coeff_ring();
        let theta = anf.generator();

        // println!("sqfree factor");
        // println!("p = {}", p);

        let mut k: usize = 0;
        let mut t;
        loop {
            // println!("k = {}", k);
            //q(x) = p(x - kθ)
            let q = rat_poly_poly_poly.evaluate(
                &p.apply_map(|c| Polynomial::constant(c.clone())),
                &Polynomial::from_coeffs(vec![
                    anf.mul(&anf.from_int(-Integer::from(k)), &theta),
                    anf.one(),
                ]),
            );
            // println!("q = {}", q);
            t = self.polynomial_norm(&q);
            // println!("t = {}", t);

            if !Polynomial::resultant(&t, &t.clone().derivative()).is_zero() {
                break;
            }

            k += 1;
        }
        // println!("t = {} k = {}", t, k);
        // println!("t = {}", t.factor().unwrap());
        // println!("{:?}", Polynomial::euclidean_gcd(t.clone(), t.clone().derivative()));

        //found k such that t(x) = N(q(x)) = N(p(x - kθ)) is now squarefree
        //the factors of p are gcd(p(x), t_i(x + kθ)) for each squarefree factor t_i of t

        let mut p_factors = vec![];
        for (ti, ti_pow) in Polynomial::<Rational>::structure()
            .factorizations()
            .to_powers(&t.factor().unwrap())
        {
            // println!("ti = {}", ti);
            debug_assert_eq!(ti_pow, &Natural::ONE);
            p_factors.push(self.euclidean_gcd(
                p.clone(),
                rat_poly_poly_poly.evaluate(
                    &ti.apply_map(|c| Polynomial::constant(Polynomial::constant(c.clone()))),
                    &Polynomial::from_coeffs(vec![
                        anf.mul(&anf.from_int(Integer::from(k)), &theta),
                        anf.one(),
                    ]),
                ),
            ));
        }

        // println!("p_factors = {:?}", p_factors);

        self.factorizations().from_unit_and_factor_powers(
            self.one(),
            p_factors
                .into_iter()
                .map(|p_factor| (p_factor, Natural::ONE))
                .collect(),
        )
    }

    pub fn factor_primitive_sqfree_by_reduced_ring(
        &self,
        p: &<Self as SetSignature>::Set,
    ) -> FactoredRingElement<<Self as SetSignature>::Set> {
        debug_assert!(!self.is_zero(p));

        /*
            input: A squarefree polynomial p(x) over K which we want to factor over K = Q[θ]
            output: A factorization of p(x) over K = Q[θ]
            idea:
                Let L = K[x] / p(x)
                Note that L is a reduced ring
                So L is product of distinct fields corresponding to the distinct factors of p
                Find with good probability an element α in L such that 1, α, α², ..., αⁿ⁻¹ span L over Q
                Find degree n-1 q(y) in Q[y] such that αⁿ = q(α)
                Factor q(y) in Q[y]
                L is product of distinct fields corresponding to the distinct factors of q since L = K[y]/q(y)
                Compute the image of x in each factor field of L
                The factors of p(x) are the minimal polynomials of x in each factor field of L as a vector space over K
        */
        // println!("p = {}", p);

        let l_reduced_ring = self.quotient_ring(p.clone());
        //n = degree over L over Q
        let k_deg = self.coeff_ring().degree();
        if k_deg == 1 {
            let (unit, factors) = Polynomial::<Rational>::from_coeffs(
                p.coeffs()
                    .into_iter()
                    .map(|c| self.coeff_ring().reduce(c).as_constant().unwrap())
                    .collect(),
            )
            .factor()
            .unwrap()
            .into_unit_and_powers();
            self.factorizations().from_unit_and_factor_powers(
                Polynomial::constant(unit),
                factors
                    .into_iter()
                    .map(|(f, f_pow)| {
                        (
                            Polynomial::<Polynomial<Rational>>::from_coeffs(
                                f.coeffs()
                                    .into_iter()
                                    .map(|c| Polynomial::constant(c.clone()))
                                    .collect(),
                            ),
                            f_pow,
                        )
                    })
                    .collect(),
            )
        } else {
            debug_assert!(k_deg >= 2);
            let p_deg = l_reduced_ring.degree();
            let n = k_deg * p_deg;
            // println!("k_deg = {}", k_deg);
            // println!("p_deg = {}", p_deg);
            // println!("n = {}", n);
            //turn elements of L into length n vectors over Q
            let l_to_vec = |x: Polynomial<Polynomial<Rational>>| {
                double_poly_to_row(p_deg, k_deg, l_reduced_ring.reduce(&x))
            };
            //turn length n vectors over Q into elements of L
            let vec_to_l = |x: Matrix<Rational>| -> Polynomial<Polynomial<Rational>> {
                row_to_double_poly(p_deg, k_deg, x)
            };

            // println!("l_reduced_ring = {:?}", l_reduced_ring);

            //find alpha in L which generates L
            let mut alpha: Polynomial<Polynomial<Rational>>;
            let mut alpha_pow_mat;
            'alpha_search: {
                //generate random alpha in L until we find one which generates L
                //each random alpha has a good probability of generating L by the primitive element theorem
                let mut rng = <rand::rngs::StdRng as rand::SeedableRng>::seed_from_u64(0);
                let mut rat_pool = vec![];
                let mut rat_gen = Rational::exhaustive_rationals();
                rat_pool.push(rat_gen.next().unwrap()); //0
                rat_pool.push(rat_gen.next().unwrap()); //1
                for rat in rat_gen {
                    // println!("rat_pool = {:?}", rat_pool);
                    for _ in 0..1 {
                        alpha = Polynomial::from_coeffs(
                            (0..p_deg)
                                .map(|_i| {
                                    Polynomial::from_coeffs(
                                        (0..k_deg)
                                            .map(|_j| {
                                                if rand::Rng::random_range(&mut rng, 0..=(n / 3))
                                                    != 0
                                                //try to keep the choice of alpha simple by choosing lots of zeros
                                                {
                                                    Rational::ZERO
                                                } else {
                                                    rand::seq::IteratorRandom::choose(
                                                        rat_pool.iter(),
                                                        &mut rng,
                                                    )
                                                    .unwrap()
                                                    .clone()
                                                }
                                            })
                                            .collect(),
                                    )
                                })
                                .collect(),
                        );
                        // println!("{}", n);
                        // println!("possible alpha = {}", alpha);

                        assert!(l_reduced_ring.equal(&alpha, &vec_to_l(l_to_vec(alpha.clone()))));

                        alpha_pow_mat = Matrix::<Rational>::join_rows(n, {
                            let mut alpha_pow = l_reduced_ring.one();
                            let mut rows = vec![];
                            for k in 0..n {
                                if k != 0 {
                                    l_reduced_ring.mul_mut(&mut alpha_pow, &alpha);
                                }
                                rows.push(l_to_vec(alpha_pow.clone()));
                            }
                            rows
                        });
                        // alpha_pow_mat.pprint();
                        // println!("{}", alpha_pow_mat.clone().primitive_part().unwrap().rank());

                        // println!("alpha_pow_mat rank = {:?}", alpha_pow_mat.rank());

                        if alpha_pow_mat.rank() == n {
                            break 'alpha_search;
                        }
                    }
                    rat_pool.push(rat);
                }

                unreachable!();
            }
            debug_assert_eq!(alpha_pow_mat.rank(), n);
            // let alpha_vec = l_to_vec(alpha.clone());
            // println!("alpha = {}", alpha);
            // println!("alpha_vec = {:?}", alpha_vec);
            // alpha_pow_mat.pprint();
            // println!("αⁿ = {}", l_reduced_ring.nat_pow(&alpha, &Natural::from(n)));
            // println!(
            //     "αⁿ = {:?}",
            //     l_to_vec(l_reduced_ring.nat_pow(&alpha, &Natural::from(n)))
            // );

            //compute q_prime(y) in Q[y] such that αⁿ = q_prime(α)
            //then q(y) = yⁿ - q_prime(y) is such that q(α) = 0
            let q_prime_row = alpha_pow_mat
                .clone()
                .row_solve(&l_to_vec(l_reduced_ring.nat_pow(&alpha, &Natural::from(n))).get_row(0))
                .unwrap();
            let q_prime = Polynomial::<Rational>::from_coeffs(
                (0..n).map(|c| q_prime_row[c].clone()).collect(),
            );
            let q = Polynomial::add(&Polynomial::var_pow(n), &q_prime.neg());
            //assert q(α) = 0
            debug_assert!(
                l_reduced_ring.is_zero(&l_reduced_ring.polynomial_ring().evaluate(
                    &q.apply_map::<Polynomial<Polynomial<Rational>>>(|c| {
                        Polynomial::from_coeffs(vec![Polynomial::from_coeffs(vec![c.clone()])])
                    }),
                    &alpha
                ))
            );

            //Let la be L represented by rational polynomials in α modulo q(x) and set up isomorphisms between l and la
            #[cfg(any())]
            let la_reduced_ring = QuotientStructure::new_ring(
                PolynomialStructure::new(Rational::structure()).into(),
                q.clone(),
            );
            let l_to_la = |x_in_l: Polynomial<Polynomial<Rational>>| -> Polynomial<Rational> {
                let x_in_q = l_to_vec(x_in_l);
                let x_in_la_vec = alpha_pow_mat.clone().row_solve(&x_in_q.get_row(0)).unwrap();
                Polynomial::from_coeffs((0..n).map(|c| x_in_la_vec[c].clone()).collect())
            };
            #[cfg(any())]
            let la_to_l = |x_in_la: Polynomial<Rational>| -> Polynomial<Polynomial<Rational>> {
                PolynomialStructure::new(l_reduced_ring.clone().into()).evaluate(
                    &x_in_la.apply_map::<Polynomial<Polynomial<Rational>>>(|c| {
                        Polynomial::from_coeffs(vec![Polynomial::from_coeffs(vec![c.clone()])])
                    }),
                    &alpha,
                )
            };

            // l_reduced_ring and la_reduced_ring are isomorphic
            // so the unique decomposition of l_reduced_ring corresponding to the distinct factors pi(x) of p(x)
            // is the same as the unique decomposition of la_reduced_ring corresponding to the distinct factors qi(y) of q(y)
            // Q[y]/qi = K[x]/pi
            // so to compute pi we embed K and x into Q[y]/qi and compute the minimal polynomial of the embedded x over the embedded K

            // println!("q = {}", q);
            // println!("q = {}", q.factor().unwrap());

            let generator = Polynomial::constant(self.coeff_ring().generator());
            let x = self.var();

            let gen_in_la = l_to_la(generator);
            let x_in_la = l_to_la(x);
            // println!("gen in la = {}", gen_in_la);
            // println!("x in la {}", &x_in_la);

            let p_factors = Polynomial::<Rational>::structure()
                .factorizations()
                .to_powers(&q.factor().unwrap())
                .into_iter()
                .map(|(qi, pow)| {
                    // println!("");

                    debug_assert_eq!(pow, &Natural::ONE);
                    let _ = pow;

                    // Q[y]/qi(y)
                    let lai_reduced_ring = Rational::structure()
                        .into_polynomial_ring()
                        .into_quotient_field_unchecked(qi.clone());

                    //pi(x) can now be computed as the degree d minimal polynomial of x in K[x]/pi(x) = Q[y]/qi(y)
                    //this is done by writing x^d as a linear combination of smaller powers of x and powers of the generator t of K
                    //a basis of lai_reduced_ring is given by
                    // 1, t, ..., t^{deg(K)-1}
                    // x, tx, ..., t^{deg(K)-1}x
                    //        ...
                    // x^{d-1}, tx^{d-1}, ..., t^{deg(K)-1}x^{d-1}
                    //lets say that lai_reduced_ring with respect to this basis is called li_reduced_ring

                    //compute the degree of the corresponding factor pi(x) of p(x)
                    let qi_deg = lai_reduced_ring.degree();
                    debug_assert_eq!(qi_deg % k_deg, 0);
                    let pi_deg = qi_deg / k_deg;

                    //the basis (1, t, ..., t^{deg(K)-1}, x, tx, ..., t^{deg(K)-1}x, ..., x^{d-1}, tx^{d-1}, ..., t^{deg(K)-1}x^{d-1}) in that order
                    let lai_basis = (0..pi_deg)
                        .flat_map(|i| {
                            (0..k_deg)
                                .map(|j| {
                                    lai_reduced_ring.mul(
                                        &lai_reduced_ring.nat_pow(&x_in_la, &Natural::from(i)),
                                        &lai_reduced_ring.nat_pow(&gen_in_la, &Natural::from(j)),
                                    )
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>();
                    // for b in lai_basis.iter() {
                    //     println!("b = {}", b);
                    //     lai_reduced_ring.to_row_vector(b).pprint();
                    // }
                    let lai_basis_mat = Matrix::join_rows(
                        qi_deg,
                        lai_basis
                            .iter()
                            .map(|b| lai_reduced_ring.to_row(b))
                            .collect::<Vec<_>>(),
                    );
                    // lai_basis_mat.pprint();
                    debug_assert_eq!(lai_basis_mat.rank(), qi_deg);

                    // println!("pi_deg = {}", pi_deg);
                    // println!("qi = {}", qi);
                    // println!("gen = {}", Polynomial::rem(&gen_in_la, qi));
                    // println!("x = {}", Polynomial::rem(&x_in_la, qi));
                    // println!(
                    //     "x^d = {}",
                    //     lai_reduced_ring.reduce(
                    //         &lai_reduced_ring.nat_pow(&x_in_la, &Natural::from(pi_deg))
                    //     )
                    // );

                    let x_wrapping_pow_vec = lai_reduced_ring
                        .to_vec(&lai_reduced_ring.nat_pow(&x_in_la, &Natural::from(pi_deg)));
                    // x_wrapping_pow_vec.pprint();
                    //this is a vector containing the coefficients of (the coefficients of elements of K of) the polynomial pi_prime(x) in K[x] such that
                    //x^n = pi_prime(x)
                    //so pi(x) = x^n - pi_prime(x) is such that pi(x) = 0 and so is the irreducible factor of p we seek such that K[x]/pi(x) = Q[y]/qi(y)
                    let x_wrapping_pow_vec_coeffs =
                        lai_basis_mat.row_solve(&x_wrapping_pow_vec).unwrap();
                    // println!("x_wrapping_pow_vec_coeffs");
                    // x_wrapping_pow_vec_coeffs.pprint();
                    let pi_prime = row_to_double_poly(
                        pi_deg,
                        k_deg,
                        Matrix::from_rows(vec![x_wrapping_pow_vec_coeffs]),
                    );
                    // println!("pi_prime = {}", pi_prime);
                    let pi = self.sub(&self.var_pow(pi_deg), &pi_prime);
                    // println!("pi = {}", pi);
                    #[allow(clippy::let_and_return)]
                    pi
                })
                .collect::<Vec<_>>();

            // println!("p factors");
            // for pi in p_factors.iter() {
            //     println!("pi = {}", pi);
            // }

            self.factorizations().from_unit_and_factor_powers(
                self.one(),
                p_factors.into_iter().map(|pi| (pi, Natural::ONE)).collect(),
            )
        }
    }

    //factor over the rationals first, then factor each irreducible rational factor over the anf
    pub fn factorize_rational_factorize_first(
        &self,
        f: &<Self as SetSignature>::Set,
        factorize: &impl Fn(
            &<Self as SetSignature>::Set,
        ) -> FactoredRingElement<<Self as SetSignature>::Set>,
    ) -> FactoredRingElement<<Self as SetSignature>::Set> {
        debug_assert!(!self.is_zero(f));
        // println!("f = {}", f);

        let rat_f = {
            let mut rat_coeffs = vec![];
            for c in f.coeffs() {
                match c.as_constant() {
                    Some(rat) => rat_coeffs.push(rat),
                    None => {
                        return factorize(f);
                    }
                }
            }
            Polynomial::<Rational>::from_coeffs(rat_coeffs)
        };

        let (rat_unit, rat_factors) = rat_f.factor().unwrap().into_unit_and_powers();
        let mut factored = self
            .factorizations()
            .from_unit(Polynomial::constant(rat_unit));
        for (rat_factor, _rat_pow) in rat_factors {
            let anf_unfactor = Polynomial::<Polynomial<Rational>>::from_coeffs(
                rat_factor
                    .coeffs()
                    .into_iter()
                    .map(|c| Polynomial::constant(c.clone()))
                    .collect(),
            );
            self.factorizations()
                .mul_mut(&mut factored, factorize(&anf_unfactor));
        }
        factored
    }
}

// #[derive(Debug, Clone, PartialEq, Eq)]
// struct IrreducibleAlgebraicNumberFieldPolynomialsStructure {
//     anf: AlgebraicNumberFieldStructure,
// }

// impl Signature for IrreducibleAlgebraicNumberFieldPolynomialsStructure {}

// impl SetSignature for IrreducibleAlgebraicNumberFieldPolynomialsStructure {
//     type Set = Polynomial<<AlgebraicNumberFieldStructure as SetSignature>::Set>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         todo!()
//     }
// }

// impl EqSignature for IrreducibleAlgebraicNumberFieldPolynomialsStructure {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         self.anf.polynomial_ring().equal(a, b)
//     }
// }

// impl OrdSignature for IrreducibleAlgebraicNumberFieldPolynomialsStructure {
//     fn cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
//         todo!()
//     }
// }

// impl<B: BorrowedStructure<AlgebraicNumberFieldStructure>> UniqueFactorizationSignature
//     for PolynomialStructure<AlgebraicNumberFieldStructure, B>
// {
//     type Irreducibles = IrreducibleAlgebraicNumberFieldPolynomialsStructure;

//     type Factorizations<SelfB: BorrowedStructure<Self>> = FactoredRingElementStructure<Self, SelfB>;

//     fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn into_factorizations(self) -> Self::Factorizations<Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn irreducibles(&self) -> impl std::borrow::Borrow<Self::Irreducibles> {
//         IrreducibleAlgebraicNumberFieldPolynomialsStructure {
//             anf: self.coeff_ring().clone(),
//         }
//     }
// }

impl<B: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>> FactorableSignature
    for PolynomialStructure<AlgebraicNumberFieldPolynomialQuotientStructure, B>
{
    fn factor(&self, a: &Self::Set) -> Option<crate::structure::FactoredRingElement<Self::Set>> {
        if self.is_zero(a) {
            None
        } else {
            Some(
                self.factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(
                    a.clone(),
                    |c| self.coeff_ring().factor(c),
                    &|a| {
                        self.factorize_rational_factorize_first(&a, &|a| {
                            // Unsure which is faster. One might be better in different cases.
                            self.factor_primitive_sqfree_by_reduced_ring(a)
                            // OR
                            // self.factor_primitive_sqfree_by_symmetric_root_polynomials(a)
                        })
                    },
                ),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomic;

    #[test]
    fn test_anf_poly_factor_count() {
        let y = &Polynomial::<Rational>::var().into_ergonomic();
        let k = (y.pow(2) - 3)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let k_poly = k.polynomial_ring();
        let x = k_poly.into_ergonomic(k_poly.var());
        debug_assert_eq!(
            k_poly
                .factorizations()
                .to_powers(&k_poly.factor(&(x.pow(2) - 12).into_verbose()).unwrap())
                .len(),
            2
        );

        let y = &Polynomial::<Rational>::var().into_ergonomic();
        let k = (y.pow(4) - y.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let k_poly = k.polynomial_ring();
        let x = k_poly.into_ergonomic(k_poly.var());
        debug_assert_eq!(
            k_poly
                .factorizations()
                .to_powers(
                    &k_poly
                        .factor(&(x.pow(4) - x.pow(2) + 1).into_verbose())
                        .unwrap()
                )
                .len(),
            4
        );

        let k = (y.pow(3) - y + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let k_poly = k.polynomial_ring();
        debug_assert_eq!(
            k_poly
                .factorizations()
                .to_powers(&k_poly.factor(&(x.pow(3) - x + 1).into_verbose()).unwrap())
                .len(),
            2
        );

        let k = (y.pow(4) - y.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let k_poly = k.polynomial_ring();
        let x = k_poly.into_ergonomic(k_poly.var());
        debug_assert_eq!(
            k_poly
                .factorizations()
                .to_powers(&k_poly.factor(&(x.pow(12) - 1).into_verbose()).unwrap())
                .len(),
            12
        );
    }
}
