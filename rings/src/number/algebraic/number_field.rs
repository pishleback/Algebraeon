
use itertools::Itertools;
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::{exhaustive::exhaustive_rationals, Rational};

use crate::{
    linear::matrix::Matrix,
    number::{
        algebraic::isolated_roots::anf_multi_primitive_element_theorem, natural::nat_to_usize,
    },
    polynomial::{multipoly::*, polynomial::*, symmetric::*},
    ring_structure::{cannonical::*, factorization::Factored, quotient::*, structure::*},
    structure::*,
};

pub type ANFStructure = QuotientStructure<PolynomialStructure<CannonicalStructure<Rational>>, true>;

pub fn new_anf(f: Polynomial<Rational>) -> ANFStructure {
    ANFStructure::new(PolynomialStructure::new(Rational::structure()).into(), f)
}

//return the splitting field and the roots of f in the splitting field
pub fn splitting_field_anf(f: &Polynomial<Rational>) -> (ANFStructure, Vec<Polynomial<Rational>>) {
    let roots = f.primitive_part_fof().all_complex_roots();
    let (g, roots_rel_g) = anf_multi_primitive_element_theorem(roots.iter().collect());
    (new_anf(g.min_poly()), roots_rel_g)
}

impl ANFStructure {
    pub fn trace_form_matrix(&self, elems: &Vec<Polynomial<Rational>>) -> Matrix<Rational> {
        let n = self.degree();
        assert_eq!(n, elems.len());
        Matrix::construct(n, n, |r, c| {
            self.trace(&Polynomial::mul(&elems[r], &elems[c]))
        })
    }

    pub fn discriminant(&self, elems: &Vec<Polynomial<Rational>>) -> Rational {
        self.trace_form_matrix(elems).det().unwrap()
    }

    pub fn compute_integral_basis(&self) -> Vec<Polynomial<Rational>> {
        //https://www.ucl.ac.uk/~ucahmki/intbasis.pdf
        // println!("compute_basis_ring_of_integers");
        let n = self.degree();
        let mut guess = (0..n)
            .map(|i| self.integral_multiple(&Polynomial::<Rational>::var_pow(i)))
            .collect_vec();

        'search: loop {
            for algint in &guess {
                debug_assert!(algint.num_coeffs() <= n); //lets keep our basis alg ints reduced
            }

            let disc = self.discriminant(&guess);
            debug_assert_eq!(disc.denominator_ref(), &Natural::ONE); //discriminant of algebraic integers is an integer
            let disc = Rational::numerator(&disc);
            debug_assert_ne!(disc, Integer::ZERO); //discriminant of a basis is non-zero
                                                   //    println!("{}", disc);
            let (_sign, mut disc_factors) = disc.factor().unwrap().unit_and_factors();
            // If p is a prime such that p^2 divides Disc
            // then can find an alg int of the form
            // 1/p (x_1a_1 + ... + x_na_n)
            // 0 <= x_i <= p - 1 and x_i in Z
            // where {a_i} is the current guess at an integral basis
            // If no algebraic integers of this form exist then we have an actual integral basis
            // If one does exist then we can add it to the integral basis & reduce to get a new guess at a basis

            // println!("guess = {:?}", guess);
            // println!("disc = {:?}", disc);
            // println!("disc_factors = {:?}", disc_factors);
            disc_factors.sort_by_key(|(p, _k)| p.clone()); //try small primes first

            for (p, k) in disc_factors {
                debug_assert!(p >= Integer::ZERO);
                let p = nat_to_usize(p.unsigned_abs_ref()).unwrap(); //if p is too big for usize then this algorithm was doomed to take longer than my lifespan anyway

                if k >= 2 {
                    // println!("p = {}", p);

                    for coeffs in (0..n).map(|_i| 0..p).multi_cartesian_product() {
                        let alpha = Polynomial::from_coeffs(
                            Polynomial::sum(
                                (0..n)
                                    .map(|i| {
                                        Polynomial::mul(
                                            &Polynomial::constant(Rational::from(coeffs[i])),
                                            &guess[i],
                                        )
                                    })
                                    .collect(),
                            )
                            .into_coeffs()
                            .into_iter()
                            .map(|c| c / Rational::from(p))
                            .collect(),
                        );

                        // println!("coeffs = {:?}  alpha = {:?}  min_poly = {}", coeffs, alpha, self.min_poly(&alpha));

                        if !self.is_zero(&alpha) && self.is_algebraic_integer(&alpha) {
                            // println!("alpha = {:?} {}", alpha, self.min_poly(&alpha));

                            guess.push(alpha);
                            let guess_mat = Matrix::construct(n + 1, n, |r, c| guess[r].coeff(c));
                            let (mul, guess_mat_prim) = guess_mat.factor_primitive_fof();
                            let guess_mat_prim_hnf = guess_mat_prim
                                .flip_cols()
                                .row_reduced_hermite_normal_form()
                                .flip_cols();

                            // guess_mat.pprint();
                            // guess_mat_prim_hnf.pprint();

                            // println!("{:?}", mul);

                            guess = (0..n)
                                .rev()
                                .map(|i| {
                                    self.from_row_vector(
                                        guess_mat_prim_hnf
                                            .get_row(i)
                                            .apply_map(|v| Rational::from(v) * &mul),
                                    )
                                })
                                .collect();

                            // println!("new_guess = {:?}", guess);
                            continue 'search;
                        }
                    }
                }
            }
            return guess;
        }
    }

    pub fn is_algebraic_integer(&self, a: &Polynomial<Rational>) -> bool {
        if self.trace(a).denominator_ref() != &Natural::ONE {
            return false;
        }
        if self.norm(a).denominator_ref() != &Natural::ONE {
            return false;
        }
        self.min_poly(a)
            .coeffs()
            .into_iter()
            .all(|c| c.denominator_ref() == &Natural::ONE)
    }

    //return a scalar multiple of a which is an algebraic integer
    fn integral_multiple(&self, a: &Polynomial<Rational>) -> Polynomial<Rational> {
        let m = Integer::lcm_list(
            self.min_poly(a)
                .coeffs()
                .into_iter()
                .map(|c| Integer::from(c.denominator_ref()))
                .collect(),
        );
        let b = Polynomial::mul(&Polynomial::constant(Rational::from(m)), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

impl CharZeroStructure for ANFStructure {}

struct RingOfIntegers {
    anf: ANFStructure,
    basis: Vec<Polynomial<Rational>>,
}

//convert between nested polynomials with bounded degree and rational vectors
//e.g. (a + bx) + (c + dc)x <-> (a, b, c, d)
fn double_poly_to_row(
    outer_poly_len: usize,
    inner_poly_len: usize,
    a: Polynomial<Polynomial<Rational>>,
) -> Matrix<Rational> {
    // let n = outer_poly_len * inner_poly_len;
    debug_assert!(a.num_coeffs() <= outer_poly_len);
    for c in a.coeffs() {
        debug_assert!(c.num_coeffs() <= inner_poly_len);
    }
    Matrix::from_rows(vec![(0..outer_poly_len)
        .map(|i| a.coeff(i))
        .map(|c| {
            (0..inner_poly_len)
                .map(|j| c.coeff(j).clone())
                .collect_vec()
        })
        .flatten()
        .collect_vec()])
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
        .collect_vec();
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

impl PolynomialStructure<ANFStructure> {
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
            .collect_vec();
        // println!("embeddings = {:?}", embedding_vars);

        let rational_poly_multipoly_structure =
            PolynomialStructure::new(MultiPolynomialStructure::new(Rational::structure()).into());

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
                .map(|c| c.clone())
                .collect_vec();

            let lc = min_poly_coeffs.pop().unwrap();
            debug_assert_eq!(lc, Rational::ONE);

            min_poly_coeffs
                .into_iter()
                .rev()
                .enumerate()
                .map(|(i, c)| if i % 2 == 0 { -c } else { c })
                .collect_vec()
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
                &(0..n)
                    .map(|i| (elem_sym_vars[i].clone(), &e_vals[i]))
                    .collect(),
            )
        });

        // println!("end");

        norm_f
    }

    pub fn factor_primitive_sqfree_by_symmetric_root_polynomials(
        &self,
        p: &<Self as Structure>::Set,
    ) -> crate::ring_structure::factorization::Factored<Self> {
        //https://www.cse.iitk.ac.in/users/nitin/courses/scribed2-WS2011-12.pdf

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
            let q = p.apply_map(|c| Polynomial::constant(c.clone())).evaluate(
                &Polynomial::from_coeffs(vec![
                    anf.mul(&anf.from_int(&-Integer::from(k)), &theta),
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
        for (ti, ti_pow) in t.factor().unwrap().factors() {
            // println!("ti = {}", ti);
            debug_assert_eq!(ti_pow, &Natural::ONE);
            p_factors.push(
                self.euclidean_gcd(
                    p.clone(),
                    ti.apply_map(|c| Polynomial::constant(Polynomial::constant(c.clone())))
                        .evaluate(&Polynomial::from_coeffs(vec![
                            anf.mul(&anf.from_int(&Integer::from(k)), &theta),
                            anf.one(),
                        ])),
                ),
            );
        }

        // println!("p_factors = {:?}", p_factors);

        let factored = Factored::new_unchecked(
            self.clone().into(),
            self.one(),
            p_factors
                .into_iter()
                .map(|p_factor| (p_factor, Natural::ONE))
                .collect(),
        );
        // println!("factored = {}", factored);
        factored
    }

    pub fn factor_primitive_sqfree_by_reduced_ring(
        &self,
        p: &<Self as Structure>::Set,
    ) -> Factored<Self> {
        debug_assert!(!self.is_zero(p));

        /*
            input: A squarefree polynomial p(x) over K which we want to factor over K = Q[θ]
            output: A factorisation of p(x) over K = Q[θ]
            idea:
                Let L = K[x] / p(x)
                Note that L is a reduced ring
                So L is product of distinct fields coresponding to the distinct factors of p
                Find with good probability an element α in L such that 1, α, α², ..., αⁿ⁻¹ span L over Q
                Find degree n-1 q(y) in Q[y] such that αⁿ = q(α)
                Factor q(y) in Q[y]
                L is product of distinct fields coresponding to the distinct factors of q since L = K[y]/q(y)
                Compute the image of x in each factor field of L
                The factors of p(x) are the minimal polynomials of x in each factor field of L as a vector space over K
        */
        // println!("p = {}", p);

        let l_reduced_ring = QuotientStructure::new_ring(self.clone().into(), p.clone());
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
            .unit_and_factors();
            Factored::new_unchecked(
                self.clone().into(),
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
                let mut rat_gen = exhaustive_rationals();
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
                                                if rand::Rng::gen_range(&mut rng, 0..(1 + n / 3))
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
                .row_solve(l_to_vec(l_reduced_ring.nat_pow(&alpha, &Natural::from(n))))
                .unwrap();
            let q_prime = Polynomial::<Rational>::from_coeffs(
                (0..n)
                    .map(|c| q_prime_row.at(0, c).unwrap().clone())
                    .collect(),
            );
            let q = Polynomial::add(&Polynomial::var_pow(n), &q_prime.neg());
            //assert q(α) = 0
            debug_assert!(l_reduced_ring.is_zero(
                &PolynomialStructure::new(l_reduced_ring.clone().into()).evaluate(
                    &q.apply_map::<Polynomial<Polynomial<Rational>>>(|c| {
                        Polynomial::from_coeffs(vec![Polynomial::from_coeffs(vec![c.clone()])])
                    }),
                    &alpha
                )
            ));

            //Let la be L represented by rational polynomials in α modulo q(x) and set up isomorphisms between l and la
            #[cfg(any())]
            let la_reduced_ring = QuotientStructure::new_ring(
                PolynomialStructure::new(Rational::structure()).into(),
                q.clone(),
            );
            let l_to_la = |x_in_l: Polynomial<Polynomial<Rational>>| -> Polynomial<Rational> {
                let x_in_q = l_to_vec(x_in_l);
                let x_in_la_vec = alpha_pow_mat.row_solve(x_in_q).unwrap();
                let x_in_la = Polynomial::from_coeffs(
                    (0..n)
                        .map(|c| x_in_la_vec.at(0, c).unwrap().clone())
                        .collect(),
                );
                x_in_la
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
            // so the unique decomposition of l_reduced_ring coresponding to the distinct factors pi(x) of p(x)
            // is the same as the unique decomposition of la_reduced_ring coresponding to the distinct factors qi(y) of q(y)
            // Q[y]/qi = K[x]/pi
            // so to compute pi we embed K and x into Q[y]/qi and compute the minimal polynomial of the embedded x over the embedded K

            // println!("q = {}", q);
            // println!("q = {}", q.factor().unwrap());

            let gen = Polynomial::constant(self.coeff_ring().generator());
            let x = self.var();

            let gen_in_la = l_to_la(gen);
            let x_in_la = l_to_la(x);
            // println!("gen in la = {}", gen_in_la);
            // println!("x in la {}", &x_in_la);

            let p_factors = q
                .factor()
                .unwrap()
                .factors()
                .into_iter()
                .map(|(qi, pow)| {
                    // println!("");

                    debug_assert_eq!(pow, &Natural::ONE);
                    let _ = pow;

                    // Q[y]/qi(y)
                    let lai_reduced_ring = QuotientStructure::new_field(
                        PolynomialStructure::new(Rational::structure()).into(),
                        qi.clone(),
                    );

                    //pi(x) can now be computed as the degree d minimal polynomial of x in K[x]/pi(x) = Q[y]/qi(y)
                    //this is done by writing x^d as a linear combination of smaller powers of x and powers of the generator t of K
                    //a basis of lai_reduced_ring is given by
                    // 1, t, ..., t^{deg(K)-1}
                    // x, tx, ..., t^{deg(K)-1}x
                    //        ...
                    // x^{d-1}, tx^{d-1}, ..., t^{deg(K)-1}x^{d-1}
                    //lets say that lai_reduced_ring with respect to this basis is called li_reduced_ring

                    //compute the degree of the coresponding factor pi(x) of p(x)
                    let qi_deg = lai_reduced_ring.degree();
                    debug_assert_eq!(qi_deg % k_deg, 0);
                    let pi_deg = qi_deg / k_deg;

                    //the basis (1, t, ..., t^{deg(K)-1}, x, tx, ..., t^{deg(K)-1}x, ..., x^{d-1}, tx^{d-1}, ..., t^{deg(K)-1}x^{d-1}) in that order
                    let lai_basis = (0..pi_deg)
                        .map(|i| {
                            (0..k_deg)
                                .map(|j| {
                                    lai_reduced_ring.mul(
                                        &lai_reduced_ring.nat_pow(&x_in_la, &Natural::from(i)),
                                        &lai_reduced_ring.nat_pow(&gen_in_la, &Natural::from(j)),
                                    )
                                })
                                .collect_vec()
                        })
                        .flatten()
                        .collect_vec();
                    // for b in lai_basis.iter() {
                    //     println!("b = {}", b);
                    //     lai_reduced_ring.to_row_vector(b).pprint();
                    // }
                    let lai_basis_mat = Matrix::join_rows(
                        qi_deg,
                        lai_basis
                            .iter()
                            .map(|b| lai_reduced_ring.to_row_vector(b))
                            .collect_vec(),
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
                        .to_row_vector(&lai_reduced_ring.nat_pow(&x_in_la, &Natural::from(pi_deg)));
                    // println!("x_wrapping_pow_vec");
                    // x_wrapping_pow_vec.pprint();
                    //this is a vector containing the coefficients of (the coefficients of elements of K of) the polynomial pi_prime(x) in K[x] such that
                    //x^n = pi_prime(x)
                    //so pi(x) = x^n - pi_prime(x) is such that pi(x) = 0 and so is the irreducible factor of p we seek such that K[x]/pi(x) = Q[y]/qi(y)
                    let x_wrapping_pow_vec_coeffs =
                        lai_basis_mat.row_solve(x_wrapping_pow_vec).unwrap();
                    // println!("x_wrapping_pow_vec_coeffs");
                    // x_wrapping_pow_vec_coeffs.pprint();
                    let pi_prime = row_to_double_poly(pi_deg, k_deg, x_wrapping_pow_vec_coeffs);
                    // println!("pi_prime = {}", pi_prime);
                    let pi = self.add(&self.var_pow(pi_deg), &pi_prime.neg());
                    // println!("pi = {}", pi);
                    pi
                })
                .collect_vec();

            // println!("p factors");
            // for pi in p_factors.iter() {
            //     println!("pi = {}", pi);
            // }

            Factored::new_unchecked(
                self.clone().into(),
                self.one(),
                p_factors.into_iter().map(|pi| (pi, Natural::ONE)).collect(),
            )
        }
    }

    //factor over the rationals first, then factor each irreducible rational factor over the anf
    pub fn factorize_rational_factorize_first(
        &self,
        f: &<Self as Structure>::Set,
        factorize: &impl Fn(&<Self as Structure>::Set) -> Factored<Self>,
    ) -> Factored<Self> {
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

        let (rat_unit, rat_factors) = rat_f.factor().unwrap().unit_and_factors();
        let mut factored =
            Factored::factored_unit_unchecked(self.clone().into(), Polynomial::constant(rat_unit));
        for (rat_factor, _rat_pow) in rat_factors {
            let anf_unfactor = Polynomial::<Polynomial<Rational>>::from_coeffs(
                rat_factor
                    .coeffs()
                    .into_iter()
                    .map(|c| Polynomial::constant(c.clone()))
                    .collect(),
            );
            factored.mul_mut(factorize(&anf_unfactor));
        }
        factored
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<ANFStructure> {
    fn factor(
        &self,
        a: &Self::Set,
    ) -> Option<crate::ring_structure::factorization::Factored<Self>> {
        if self.is_zero(a) {
            None
        } else {
            Some(
                self.factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(a.clone(), &|a| {
                    self.factorize_rational_factorize_first(&a, &|a| {
                        self.factor_primitive_sqfree_by_reduced_ring(a)
                        // self.factor_primitive_sqfree_by_symmetric_root_polynomials(a)
                    })
                }),
            )
        }
    }
}

#[cfg(test)]
mod tests {
    

    use super::*;

    #[test]
    fn test_anf_to_and_from_vector() {
        let x = &Polynomial::<Rational>::var().into_ring();
        let anf = new_anf((x.pow(5) - x + 1).into_set());
        let alpha = (x.pow(9) + 5).into_set();

        println!("{}", alpha);
        println!("{}", anf.reduce(&alpha));
        println!("{}", anf.min_poly(&alpha));
        anf.to_col_vector(&alpha).pprint();

        assert_eq!(
            anf.to_col_vector(&alpha),
            Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])
        );

        assert!(anf.equal(
            &anf.from_col_vector(Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])),
            &alpha
        ));
    }

    #[test]
    fn test_anf_poly_factor_count() {
        let y = &Polynomial::<Rational>::var().into_ring();
        let k = new_anf((y.pow(2) - 3).into_set());
        let k_poly = PolynomialStructure::new(k.clone().into());
        let x = &k_poly.var().into_ring();
        debug_assert_eq!(
            k_poly
                .factor(&(x.pow(2) - 12).into_set())
                .unwrap()
                .factors()
                .len(),
            2
        );

        let y = &Polynomial::<Rational>::var().into_ring();
        let k = new_anf((y.pow(4) - y.pow(2) + 1).into_set());
        let k_poly = PolynomialStructure::new(k.clone().into());
        let x = &k_poly.var().into_ring();
        debug_assert_eq!(
            k_poly
                .factor(&(x.pow(4) - x.pow(2) + 1).into_set())
                .unwrap()
                .factors()
                .len(),
            4
        );

        let k = new_anf((y.pow(3) - y + 1).into_set());
        let k_poly = PolynomialStructure::new(k.clone().into());
        debug_assert_eq!(
            k_poly
                .factor(&(x.pow(3) - x + 1).into_set())
                .unwrap()
                .factors()
                .len(),
            2
        );

        let k = new_anf((y.pow(4) - y.pow(2) + 1).into_set());
        let k_poly = PolynomialStructure::new(k.clone().into());
        let x = &k_poly.var().into_ring();
        debug_assert_eq!(
            k_poly
                .factor(&(x.pow(12) - 1).into_set())
                .unwrap()
                .factors()
                .len(),
            12
        );
    }
}
