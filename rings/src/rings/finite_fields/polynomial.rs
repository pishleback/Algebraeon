use crate::{linear::matrix::*, polynomial::*, rings::quotient::QuotientStructure, structure::*};
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;

// Useful: https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
/*
Factorization over finite fields typically happens in the following steps:
1. Monic factorization: Remove a unit so that the polynomial to be factored is monic.
2. Squarefree factorization: Factor into squarefree factors.
3. Distinct degree factorization: Factor a squarefree polynomial into polynomials such that all irreducible factors of each polynomial have equal degree.
4. Equal degree factorization: Factor a squarefree polynomial each of whose irreducible factors all have the same known degree d.

1 is trivial.
There is a standard algorithm for 2.
Berlekamps algorithm does 3 and 4 at the same time.
There is a standard algorithm for 3.
Cantor–Zassenhaus algorithm does 4.
*/

/// Store a monic factorization
#[derive(Debug, Clone)]
pub struct MonicFactored<FS: FieldSignature, FSB: BorrowedStructure<FS>> {
    poly_ring: PolynomialStructure<FS, FSB>,
    unit: FS::Set,              // a unit
    monic: Polynomial<FS::Set>, // a monic polynomial
}

/// Store a squarefree factorization
#[derive(Debug, Clone)]
pub struct SquarefreeFactored<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> {
    poly_ring: PolynomialStructure<FS, FSB>,
    unit: FS::Set,                                           // a unit
    squarefree_factors: Vec<(Polynomial<FS::Set>, Natural)>, // squarefree monic polynomials and their multiplicities
}

#[derive(Debug, Clone)]
struct DistinctDegreeFactor<FS: FiniteFieldSignature> {
    irreducible_factor_degree: usize, // the degree of the irreducible factors
    polynomial: Polynomial<FS::Set>, // a monic squarefree polynomial equal to a product of irreducibles all of the same degree
}
/// Store a distinct degree factorization
#[derive(Debug, Clone)]
pub struct DistinctDegreeFactored<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> {
    poly_ring: PolynomialStructure<FS, FSB>,
    unit: FS::Set, // a unit
    distinct_degree_factors: Vec<(DistinctDegreeFactor<FS>, Natural)>,
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> MonicFactored<FS, FSB> {
    // pub fn factor(poly_ring: PolynomialStructure<FS>, poly: Polynomial<FS::Set>) -> Self {
    //     let (unit, monic) = poly_ring.factor_fav_assoc(&poly);
    //     let unit = poly_ring.as_constant(&unit).unwrap();
    //     debug_assert!(poly_ring.is_monic(&monic));
    //     Self {
    //         poly_ring,
    //         unit,
    //         monic,
    //     }
    // }

    pub fn new_monic_unchecked(
        poly_ring: PolynomialStructure<FS, FSB>,
        monic: Polynomial<FS::Set>,
    ) -> Self {
        debug_assert!(poly_ring.is_monic(&monic));
        let unit = poly_ring.coeff_ring().one();
        Self {
            poly_ring,
            unit,
            monic,
        }
    }

    pub fn scalar_part(&self) -> &FS::Set {
        &self.unit
    }

    pub fn monic_part(&self) -> &Polynomial<FS::Set> {
        &self.monic
    }

    pub fn into_polynomial(self) -> Polynomial<FS::Set> {
        self.poly_ring
            .mul(&Polynomial::constant(self.unit), &self.monic)
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> SquarefreeFactored<FS, FSB> {
    pub fn unit_unchecked(poly_ring: PolynomialStructure<FS, FSB>, unit: FS::Set) -> Self {
        debug_assert!(poly_ring.coeff_ring().is_unit(&unit));
        Self {
            poly_ring,
            unit,
            squarefree_factors: vec![],
        }
    }

    pub fn new_squarefree_poly_unchecked(
        poly_ring: PolynomialStructure<FS, FSB>,
        poly: Polynomial<FS::Set>,
    ) -> Self {
        debug_assert!(!poly_ring.is_zero(&poly));
        let deg = poly_ring.degree(&poly).unwrap();
        if deg == 0 {
            let unit = poly_ring.as_constant(&poly).unwrap();
            Self {
                poly_ring,
                unit,
                squarefree_factors: vec![],
            }
        } else {
            let unit = poly_ring.coeff_ring().one();
            Self {
                poly_ring,
                unit,
                squarefree_factors: vec![(poly, Natural::ONE)],
            }
        }
    }

    pub fn pow(mut self, power: &Natural) -> Self {
        for (_, poly_power) in &mut self.squarefree_factors {
            *poly_power *= power;
        }
        self
    }

    pub fn into_monic_factored(self) -> MonicFactored<FS, FSB> {
        let monic = self.poly_ring.product(
            self.squarefree_factors
                .iter()
                .map(|(f, k)| self.poly_ring.nat_pow(f, k))
                .collect(),
        );
        MonicFactored {
            poly_ring: self.poly_ring,
            unit: self.unit,
            monic,
        }
    }

    fn mul(&mut self, other: &Self) {
        assert_eq!(self.poly_ring, other.poly_ring);
        let poly_ring = &self.poly_ring;
        self.unit = poly_ring.coeff_ring().mul(&self.unit, &other.unit);
        'OTHER_LOOP: for (poly, power) in &other.squarefree_factors {
            for (self_poly, self_power) in &mut self.squarefree_factors {
                if poly_ring.equal(&poly, self_poly) {
                    *self_power += power;
                    continue 'OTHER_LOOP;
                }
            }
            self.squarefree_factors.push((poly.clone(), power.clone()));
        }
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> DistinctDegreeFactored<FS, FSB> {
    pub fn into_squarefree_factored(self) -> SquarefreeFactored<FS, FSB> {
        let squarefree_factors = self
            .distinct_degree_factors
            .into_iter()
            .map(|(ddf, k)| (ddf.polynomial, k))
            .collect();
        SquarefreeFactored {
            poly_ring: self.poly_ring,
            unit: self.unit,
            squarefree_factors,
        }
    }
}

impl<
    FS: FiniteFieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
> FactoredRingElementStructure<PolynomialStructure<FS, FSB>, FSPB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>> + FactorableSignature,
{
    pub fn into_distinct_degree_factored(
        &self,
        a: FactoredRingElement<Polynomial<FS::Set>>,
    ) -> DistinctDegreeFactored<FS, FSB> {
        let poly_ring = self.ring().clone();
        let (unit, factors) = a.into_unit_and_factor_powers();
        let unit = poly_ring.as_constant(&unit).unwrap();
        let distinct_degree_factors = factors
            .into_iter()
            .map(|(f, k)| {
                (
                    DistinctDegreeFactor {
                        irreducible_factor_degree: poly_ring.degree(&f).unwrap(),
                        polynomial: f,
                    },
                    k,
                )
            })
            .collect();
        DistinctDegreeFactored {
            poly_ring,
            unit,
            distinct_degree_factors,
        }
    }
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> PolynomialStructure<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    /// monic factorization
    pub fn factorize_monic(&self, poly: &Polynomial<FS::Set>) -> Option<MonicFactored<FS, FSB>> {
        if self.is_zero(poly) {
            None
        } else {
            let (unit, monic) = self.factor_fav_assoc(poly);
            let unit = self.as_constant(&unit).unwrap();
            Some(MonicFactored {
                poly_ring: self.clone(),
                unit,
                monic,
            })
        }
    }
}

impl<F: MetaType> Polynomial<F>
where
    F::Signature: FiniteFieldSignature,
    PolynomialStructure<F::Signature, F::Signature>: SetSignature<Set = Polynomial<F>>,
{
    pub fn factorize_monic(&self) -> Option<MonicFactored<F::Signature, F::Signature>> {
        Self::structure().factorize_monic(self)
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> MonicFactored<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    /// squarefree factorization
    pub fn factorize_squarefree(&self) -> SquarefreeFactored<FS, FSB> {
        let mut factors =
            SquarefreeFactored::unit_unchecked(self.poly_ring.clone(), self.unit.clone());

        //let w be the squarefree product of all factors of f of multiplicity not divisible by p
        let mut c = self
            .poly_ring
            .gcd(&self.monic, &self.poly_ring.derivative(self.monic.clone()));
        let mut w = self.poly_ring.div(&self.monic, &c).unwrap();

        //find all the factors in w
        let mut i = Natural::ONE;
        while !self.poly_ring.equal(&w, &self.poly_ring.one()) {
            let y = self.poly_ring.gcd(&w, &c);
            factors.mul(
                &SquarefreeFactored::new_squarefree_poly_unchecked(
                    self.poly_ring.clone(),
                    self.poly_ring.div(&w, &y).unwrap(),
                )
                .pow(&i),
            );
            let c_over_y = self.poly_ring.div(&c, &y).unwrap();
            (w, c) = (y, c_over_y);
            i += Natural::ONE;
        }

        // c is now the product, with multiplicity, of the remaining factors of f
        if !self.poly_ring.equal(&c, &self.poly_ring.one()) {
            //c = c^{1/p}
            let p = self.poly_ring.coeff_ring().characteristic_and_power().0;
            let mut reduced_c_coeffs = vec![];
            for (k, coeff) in c.coeffs().into_iter().enumerate() {
                if Natural::from(k) % &p == Natural::ZERO {
                    reduced_c_coeffs.push(coeff.clone());
                } else {
                    debug_assert!(self.poly_ring.coeff_ring().is_zero(coeff));
                }
            }
            let reduced_c: Polynomial<<FS as SetSignature>::Set> =
                Polynomial::from_coeffs(reduced_c_coeffs);
            factors.mul(
                &MonicFactored::new_monic_unchecked(self.poly_ring.clone(), reduced_c)
                    .factorize_squarefree()
                    .pow(&p),
            );
        }
        factors
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> PolynomialStructure<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    fn find_factor_by_berlekamps_algorithm(
        &self,
        f: Polynomial<FS::Set>,
    ) -> FindFactorResult<Polynomial<FS::Set>> {
        debug_assert!(self.is_squarefree(&f));

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
            .collect::<Vec<_>>();
        let mat = Matrix::construct(f_deg, f_deg, |i, j| self.coeff(&row_polys[i], j).clone());
        // mat.pprint();
        //the column kernel gives a basis of berlekamp subalgebra - all polynomials g such that g^q=g
        let mat_struct = MatrixStructure::<FS, _>::new(self.coeff_ring());
        let ker = mat_struct.row_kernel(mat);
        // ker.pprint();
        let ker_rank = ker.rank();
        let ker_basis = ker
            .basis()
            .into_iter()
            .map(|col| Polynomial::from_coeffs((0..f_deg).map(|c| col[c].clone()).collect()))
            .collect::<Vec<_>>();

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
                #[allow(clippy::single_match_else)]
                match self.div(&f, &g) {
                    Ok(g_prime) => {
                        return FindFactorResult::Composite(g, g_prime);
                    }
                    Err(_) => panic!(),
                }
            }
        }
        FindFactorResult::Irreducible
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> SquarefreeFactored<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    /// use Berlekamps algorithm for a full factorization from a squarefree
    pub fn factorize_berlekamps(&self) -> FactoredRingElement<Polynomial<FS::Set>>
    where
        PolynomialStructure<FS, FSB>: FactorableSignature,
    {
        let mut factors = self
            .poly_ring
            .factorizations()
            .from_unit_and_factor_powers(Polynomial::constant(self.unit.clone()), vec![]);
        for (sqfree_poly, power) in &self.squarefree_factors {
            self.poly_ring.factorizations().mul_mut(
                &mut factors,
                self.poly_ring.factorizations().pow(
                    factorize_by_find_factor(&self.poly_ring, sqfree_poly.clone(), &|f| {
                        self.poly_ring.find_factor_by_berlekamps_algorithm(f)
                    }),
                    power,
                ),
            );
        }
        factors
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> SquarefreeFactored<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    /// distinct degree factorization
    pub fn factorize_distinct_degree(&self) -> DistinctDegreeFactored<FS, FSB> {
        // https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Distinct-degree_factorization
        let (p, k) = self.poly_ring.coeff_ring().characteristic_and_power();
        let q = p.nat_pow(&k);
        let mut distinct_degree_factors = vec![];
        for (poly, sqfree_poly_multiplicity) in &self.squarefree_factors {
            // a key step in this algorithm is the computation of g = gcd(f, x^{q^i} - x)
            // the naive approach isn't great since x^{q^i} gets quite big
            // instead the GCD can be computed after first reducing x^{q^i} modulo poly, since f divides poly
            // to quickly compute x^{q^i} we can compute x^{q^1}, x^{q^2}, ... in sucessive loops, each time raising the previous value to the power of q
            // raising polynomials mod poly to the power of q is a linear map (freshmans dream over finite field)
            // so the qth power can be obtained by pre-computing the qth power matrix for polynomials mod poly

            let n = self.poly_ring.degree(poly).unwrap();
            debug_assert!(n >= 1);

            let mod_poly_ring =
                QuotientStructure::new_ring(self.poly_ring.clone().into(), poly.clone());
            let mat_structure = MatrixStructure::new(self.poly_ring.coeff_ring().clone());
            let xq = mod_poly_ring.nat_pow(&self.poly_ring.var(), &q);
            let qth_power_matrix = Matrix::join_cols(
                n,
                (0..n)
                    .map(|c| {
                        //compute (x^c)^q mod poly as a length n column vector of coefficients
                        mod_poly_ring
                            .to_col_vector(&mod_poly_ring.nat_pow(&self.poly_ring.var_pow(c), &q))
                    })
                    .collect(),
            );

            let mut i = 1;
            let mut xqi = xq.clone(); // x^{q^i} mod poly
            let mut f = poly.clone();
            while self.poly_ring.degree(&f).unwrap() >= 2 * i {
                debug_assert!(self.poly_ring.is_monic(&f));

                let g = self
                    .poly_ring
                    .factorize_monic(
                        &self.poly_ring.gcd_by_primitive_subresultant(
                            f.clone(),
                            self.poly_ring
                                .add(&xqi, &self.poly_ring.neg(&self.poly_ring.var())),
                        ),
                    )
                    .unwrap()
                    .monic;
                debug_assert!(self.poly_ring.is_monic(&g));

                #[cfg(debug_assertions)]
                {
                    debug_assert!(
                        mod_poly_ring.equal(
                            &xqi,
                            &self
                                .poly_ring
                                .var_pow(q.nat_pow(&i.into()).try_into().unwrap())
                        )
                    );
                    let g_naive = self
                        .poly_ring
                        .factorize_monic(
                            &self.poly_ring.gcd_by_primitive_subresultant(
                                f.clone(),
                                self.poly_ring.add(
                                    &self
                                        .poly_ring
                                        .var_pow(q.nat_pow(&i.into()).try_into().unwrap()),
                                    &self.poly_ring.neg(&self.poly_ring.var()),
                                ),
                            ),
                        )
                        .unwrap()
                        .monic;
                    // This GCD was the naive and expensive way to do it. Use it to verify the fast and fancy way.
                    debug_assert!(self.poly_ring.is_monic(&g_naive));
                    debug_assert!(self.poly_ring.equal(&g, &g_naive));
                }

                if !self.poly_ring.equal(&g, &self.poly_ring.one()) {
                    f = self.poly_ring.div(&f, &g).unwrap();
                    distinct_degree_factors.push((
                        DistinctDegreeFactor {
                            irreducible_factor_degree: i,
                            polynomial: g,
                        },
                        sqfree_poly_multiplicity.clone(),
                    ));
                }
                i += 1;
                xqi = mod_poly_ring.from_col_vector(
                    mat_structure
                        .mul(&qth_power_matrix, &mod_poly_ring.to_col_vector(&xqi))
                        .unwrap(),
                );
            }
            if !self.poly_ring.equal(&f, &self.poly_ring.one()) {
                distinct_degree_factors.push((
                    DistinctDegreeFactor {
                        irreducible_factor_degree: self.poly_ring.degree(&f).unwrap(),
                        polynomial: f,
                    },
                    sqfree_poly_multiplicity.clone(),
                ));
            }
        }
        DistinctDegreeFactored {
            poly_ring: self.poly_ring.clone(),
            unit: self.unit.clone(),
            distinct_degree_factors,
        }
    }
}

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> DistinctDegreeFactored<FS, FSB>
where
    PolynomialStructure<FS, FSB>: SetSignature<Set = Polynomial<FS::Set>>,
{
    /// Cantor–Zassenhaus algorithm for equal degree factorization
    pub fn factorize_cantor_zassenhaus(&self) -> FactoredRingElement<Polynomial<FS::Set>>
    where
        PolynomialStructure<FS, FSB>: FactorableSignature,
    {
        let poly_ring = &self.poly_ring;
        let mut fs = poly_ring
            .factorizations()
            .from_unit(Polynomial::constant(self.unit.clone()));
        for (ddf, mult) in &self.distinct_degree_factors {
            let d = ddf.irreducible_factor_degree;
            let n = self.poly_ring.degree(&ddf.polynomial).unwrap();
            debug_assert_eq!(n % d, 0);

            let finite_field = self.poly_ring.coeff_ring();
            let (p, k) = finite_field.characteristic_and_power();
            let q = p.nat_pow(&k);
            let mut prand_elements = finite_field.generate_random_elements(0);

            let mut to_factor = vec![ddf.polynomial.clone()];
            loop {
                // println!(
                //     "{} {:?}",
                //     d,
                //     to_factor
                //         .iter()
                //         .map(|u| self.poly_ring.degree(u).unwrap())
                //         .collect::<Vec<_>>()
                // );
                // Any polynomial in to_factor of degree d is irreducible.
                to_factor = to_factor
                    .into_iter()
                    .filter_map(|f| {
                        if self.poly_ring.degree(&f).unwrap() == d {
                            poly_ring.factorizations().mul_mut(
                                &mut fs,
                                poly_ring
                                    .factorizations()
                                    .pow(poly_ring.factorizations().new_prime(f), mult),
                            );
                            None
                        } else {
                            Some(f)
                        }
                    })
                    .collect();

                // Is there anything left to factor?
                if to_factor.is_empty() {
                    break;
                }

                // Try to factor what's left using Cantor-Zassenhaus
                // let h be a random polynomial of degree <n
                let h = Polynomial::<FS::Set>::from_coeffs(
                    (0..n).map(|_| prand_elements.next().unwrap()).collect(),
                );
                let g = if p == Natural::TWO {
                    // when char = 2 use h + h^2 + h^4 + ... + h^{2^{kd-1}}
                    // https://math.stackexchange.com/questions/1636518/how-do-i-apply-the-cantor-zassenhaus-algorithm-to-mathbbf-2
                    let mut sum = self.poly_ring.zero();
                    let mut square_powers = h.clone();
                    let mut square_pow = 0usize;
                    while Natural::from(square_pow) < &k * Natural::from(d) {
                        self.poly_ring.add_mut(&mut sum, &square_powers);
                        square_powers = self.poly_ring.mul(&square_powers, &square_powers);
                        square_pow += 1;
                    }
                    sum
                } else {
                    // when char != 2 use h^{(q^d-1)/2}-1 mod f
                    let poly_mod_f =
                        QuotientStructure::new_ring(self.poly_ring.clone(), ddf.polynomial.clone());
                    let a = (q.nat_pow(&d.into()) - Natural::ONE) / Natural::TWO;
                    poly_mod_f.add(
                        &poly_mod_f.nat_pow(&h, &a),
                        &poly_mod_f.neg(&poly_mod_f.one()),
                    )
                };
                to_factor = to_factor
                    .into_iter()
                    .flat_map(|u| {
                        let gcd = self
                            .poly_ring
                            .factorize_monic(&self.poly_ring.subresultant_gcd(u.clone(), g.clone()))
                            .unwrap()
                            .monic;
                        let gcd_deg = self.poly_ring.degree(&gcd).unwrap();
                        if gcd_deg == 0 || gcd_deg == self.poly_ring.degree(&u).unwrap() {
                            vec![u]
                        } else {
                            vec![self.poly_ring.div(&u, &gcd).unwrap(), gcd]
                        }
                    })
                    .collect();
            }
        }
        fs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rings::finite_fields::{modulo::*, quaternary_field::QuaternaryField};

    #[test]
    fn test_distinct_degree_and_cantor_zassenhaus_factorization_f2() {
        let x = &Polynomial::<Modulo<2>>::var().into_ergonomic();
        let p = (x.pow(3) + x.pow(2) + 2) * (x.pow(3) + 2 * x.pow(2) + 1);
        let p = p.into_verbose();

        assert!(
            Polynomial::<Modulo<2>>::structure().factorizations().equal(
                &p.factorize_monic()
                    .unwrap()
                    .factorize_squarefree()
                    .factorize_berlekamps(),
                &p.factorize_monic()
                    .unwrap()
                    .factorize_squarefree()
                    .factorize_distinct_degree()
                    .factorize_cantor_zassenhaus()
            )
        );
    }

    #[test]
    fn test_distinct_degree_and_cantor_zassenhaus_factorization_f3() {
        let x = &Polynomial::<Modulo<3>>::var().into_ergonomic();
        let p = (x.pow(3) + x.pow(2) + 2) * (x.pow(3) + 2 * x.pow(2) + 1);
        let p = p.into_verbose();

        assert!(
            Polynomial::<Modulo<3>>::structure().factorizations().equal(
                &p.factorize_monic()
                    .unwrap()
                    .factorize_squarefree()
                    .factorize_berlekamps(),
                &p.factorize_monic()
                    .unwrap()
                    .factorize_squarefree()
                    .factorize_distinct_degree()
                    .factorize_cantor_zassenhaus()
            )
        );
    }

    #[test]
    fn test_factorize_over_f2_example1() {
        let x = &Polynomial::<Modulo<2>>::var().into_ergonomic();
        let p = (1 + x.pow(4) + x.pow(5)).pow(12).into_verbose();
        let ans = Polynomial::<Modulo<2>>::structure()
            .factorizations()
            .from_unit_and_factor_powers(
                Polynomial::one(),
                vec![
                    ((1 + x + x.pow(2)).into_verbose(), Natural::from(12u32)),
                    ((1 + x + x.pow(3)).into_verbose(), Natural::from(12u32)),
                ],
            );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<2>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<2>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );
    }

    #[test]
    fn test_factorize_over_f2_example2() {
        let x = &Polynomial::<Modulo<2>>::var().into_ergonomic();
        let p =
            ((1 + x.pow(4) + x.pow(7)).pow(6) * (1 + x.pow(6) + x.pow(7)).pow(4)).into_verbose();
        let ans = Polynomial::<Modulo<2>>::structure()
            .factorizations()
            .from_unit_and_factor_powers(
                Polynomial::one(),
                vec![
                    (
                        (1 + x.pow(4) + x.pow(7)).into_verbose(),
                        Natural::from(6u32),
                    ),
                    (
                        (1 + x.pow(6) + x.pow(7)).into_verbose(),
                        Natural::from(4u32),
                    ),
                ],
            );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<2>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<2>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );
    }

    #[test]
    fn test_factorize_over_f5_example1() {
        let x = &Polynomial::<Modulo<5>>::var().into_ergonomic();
        let p = (1 + x.pow(4)).pow(5).into_verbose();
        let fs = Polynomial::<Modulo<5>>::structure().into_factorizations();
        let ans = fs.from_unit_and_factor_powers(
            Polynomial::one(),
            vec![
                ((2 + x.pow(2)).into_verbose(), Natural::from(5u8)),
                ((3 + x.pow(2)).into_verbose(), Natural::from(5u8)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(fs.equal(&f, &ans));

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(fs.equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example2() {
        let x = &Polynomial::<Modulo<5>>::var().into_ergonomic();
        let p = (3 + 2 * x.pow(2) + x.pow(4) + x.pow(6)).into_verbose();
        let ans = Polynomial::<Modulo<5>>::structure()
            .factorizations()
            .from_unit_and_factor_powers(
                Polynomial::one(),
                vec![((2 + x.pow(2)).into_verbose(), Natural::from(3u8))],
            );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<5>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<5>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );
    }

    #[test]
    fn test_factorize_over_f31_example1() {
        let x = &Polynomial::<Modulo<31>>::var().into_ergonomic();
        let p = (1 + x.pow(27) + 8 * x.pow(30)).into_verbose();
        let ans = Polynomial::<Modulo<31>>::structure()
            .factorizations()
            .from_unit_and_factor_powers(
                Polynomial::constant(Modulo::from_int(Integer::from(8))),
                vec![
                    ((12 + x.pow(3)).into_verbose(), Natural::from(1u32)),
                    (
                        (25 + 27 * x + 3 * x.pow(2) + 3 * x.pow(3) + 29 * x.pow(4) + x.pow(5))
                            .into_verbose()
                            .clone(),
                        Natural::from(1u32),
                    ),
                    (
                        (1 + 24 * x + 3 * x.pow(2) + 15 * x.pow(3) + 12 * x.pow(4) + x.pow(5))
                            .into_verbose()
                            .clone(),
                        Natural::from(1u32),
                    ),
                    (
                        (21 + 12 * x.pow(3) + 22 * x.pow(6) + 4 * x.pow(9) + x.pow(12))
                            .into_verbose()
                            .clone(),
                        Natural::from(1u32),
                    ),
                    (
                        (5 + 11 * x + 3 * x.pow(2) + 13 * x.pow(3) + 21 * x.pow(4) + x.pow(5))
                            .into_verbose()
                            .clone(),
                        Natural::from(1u32),
                    ),
                ],
            );

        println!("{:?}", p.factorize_monic());

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<Modulo<31>>::structure()
                .factorizations()
                .equal(&f, &ans)
        );
    }

    #[test]
    fn test_factorize_over_f4_example1() {
        let x = &Polynomial::<QuaternaryField>::var().into_ergonomic();

        let a = x - Polynomial::constant(QuaternaryField::One).into_ergonomic();
        let b = x - Polynomial::constant(QuaternaryField::Alpha).into_ergonomic();
        let c = x - Polynomial::constant(QuaternaryField::Beta).into_ergonomic();
        let p = (1 - x.pow(3)).pow(48).into_verbose();
        let ans = Polynomial::<QuaternaryField>::structure()
            .factorizations()
            .from_unit_and_factor_powers(
                Polynomial::one(),
                vec![
                    (a.into_verbose(), Natural::from(48u32)),
                    (b.into_verbose(), Natural::from(48u32)),
                    (c.into_verbose(), Natural::from(48u32)),
                ],
            );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<QuaternaryField>::structure()
                .factorizations()
                .equal(&f, &ans)
        );

        let f = p
            .factorize_monic()
            .unwrap()
            .factorize_squarefree()
            .factorize_berlekamps();
        println!("{} = {}", p, f);
        assert!(
            Polynomial::<QuaternaryField>::structure()
                .factorizations()
                .equal(&f, &ans)
        );
    }
}
