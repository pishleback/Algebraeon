use crate::matrix::{ComplexInnerProduct, RealSymmetricInnerProduct, SymmetricMatrix};
use crate::structure::{
    AdditionSignature, AdditiveGroupSignature, CancellativeMultiplicationSignature,
    EuclideanDivisionSignature, MetaCancellativeMultiplicationSignature, SemiModuleSignature,
};
use crate::{
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature, structure::ZeroEqSignature,
};
use crate::{
    matrix::{
        Matrix, MatrixStructure, RealInnerProduct,
        row_operations::{ElementaryOpp, ElementaryOppType},
    },
    structure::{FieldSignature, OrderedRingSignature, RealRoundingSignature, RealSubsetSignature},
};
use algebraeon_nzq::traits::Fraction;
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Rational};
use algebraeon_sets::structure::EqSignature;
use algebraeon_sets::structure::{BorrowedStructure, MetaType, SetSignature, ToStringSignature};

impl<
    FS: RealSubsetSignature
        + RealRoundingSignature
        + OrderedRingSignature
        + FieldSignature
        + ToStringSignature,
    FSB: BorrowedStructure<FS>,
> MatrixStructure<FS, FSB>
{
    /// Take the rows of M = `mat` as the basis for a lattice.
    ///
    /// Perform an LLL reduction on the lattice, returning matricies (H, H*M) such that H is invertible the rows of H*M are an LLL reduced basis for the lattice.
    ///
    /// `delta` must be in the range (-0.25, 1] with polynomial time complexity only for `delta` in (-0.25, 1). `delta` = 3/4 is the "default" value to use.
    pub fn lll_row_reduction_algorithm(
        &self,
        mut basis: Matrix<FS::Set>,
        inner_product: &impl RealInnerProduct<FS>,
        delta: &Rational,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let n = basis.rows(); // size of the basis
        let m = basis.cols(); // dimension of the ambient space
        assert!(m >= n);
        debug_assert!(self.rank(basis.clone()) == n);
        // 1/4 < delta <= 1
        assert!(Rational::ONE < Rational::from(4) * delta && delta <= &Rational::ONE);
        let vs = self.ring().free_module(m);

        if n == 0 {
            return (self.ident(0), basis);
        }

        #[cfg(debug_assertions)]
        let orig_basis = basis.clone();

        #[derive(Debug)]
        struct CacheEntry<Set> {
            // for i < k the ith entry is m_{k,i} = <bj, bi>/<bi, bi>
            mu: Vec<Set>,
            // the kth gram-schmidt reduced basis vector bk*
            basis_gs: Vec<Set>,
            // <bk, bk>
            basis_gs_sq: Set,
        }
        let mut cache = Vec::<CacheEntry<FS::Set>>::with_capacity(n);

        #[cfg(debug_assertions)]
        let validate_cache = |basis: &Matrix<FS::Set>, cache: &Vec<CacheEntry<FS::Set>>| {
            let true_mat_gs = self.gram_schmidt_row_orthogonalization(basis.clone(), inner_product);
            for (k, cache_entry) in cache.iter().enumerate() {
                // check mu
                assert_eq!(cache_entry.mu.len(), k);
                for i in 0..k {
                    assert!(
                        self.ring().equal(
                            &cache_entry.mu[i],
                            &self
                                .ring()
                                .try_divide(
                                    &inner_product
                                        .inner_product(&basis.get_row(k), &true_mat_gs.get_row(i),),
                                    &inner_product.inner_product(
                                        &true_mat_gs.get_row(i),
                                        &true_mat_gs.get_row(i),
                                    )
                                )
                                .unwrap()
                        )
                    );
                }

                // check squared norm
                assert!(self.ring().equal(
                    &cache_entry.basis_gs_sq,
                    &inner_product.inner_product(&cache_entry.basis_gs, &cache_entry.basis_gs)
                ));
                assert!(!self.ring().is_zero(&cache_entry.basis_gs_sq));

                // check gs-reduction
                assert!(vs.is_element(&cache_entry.basis_gs).is_ok());
                assert!(vs.equal(&cache_entry.basis_gs, &true_mat_gs.get_row(k)));
            }
        };

        cache.push(CacheEntry {
            mu: vec![],
            basis_gs: basis.get_row(0),
            basis_gs_sq: inner_product.inner_product(&basis.get_row(0), &basis.get_row(0)),
        });

        #[cfg(debug_assertions)]
        validate_cache(&basis, &cache);

        let mut h = self.ident(n);

        let reduce = |cache: &mut Vec<CacheEntry<FS::Set>>,
                      h: &mut Matrix<FS::Set>,
                      basis: &mut Matrix<FS::Set>,
                      k: usize,
                      j: usize| {
            let mu_rounded = self.ring().round(&cache[k].mu[j]);
            // b_k -= mu_rounded * b_i
            if mu_rounded != Integer::ZERO {
                let mu_rounded = self.ring().from_int(&mu_rounded);

                // update basis and basis transformation matrix
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul {
                        i: k,
                        j,
                        x: self.ring().neg(&mu_rounded),
                    },
                );
                // println!("{:?}", row_opp);
                row_opp.apply(basis);
                row_opp.apply(h);

                // update cache
                self.ring().sub_mut(&mut cache[k].mu[j], &mu_rounded);
                for i in 0..j {
                    cache[k].mu[i] = self.ring().sub(
                        &cache[k].mu[i],
                        &self.ring().mul(&mu_rounded, &cache[j].mu[i]),
                    );
                }

                #[cfg(debug_assertions)]
                validate_cache(basis, cache);
            }
        };

        let mut k = 1;
        while k < n {
            // extend the cache if necessary
            if k >= cache.len() {
                debug_assert_eq!(k, cache.len());
                let mut mu_k = Vec::with_capacity(k);
                let mut basis_gs = basis.get_row(k);
                #[allow(clippy::needless_range_loop)]
                for i in 0..k {
                    let mu_k_i = self
                        .ring()
                        .try_divide(
                            &inner_product.inner_product(&basis.get_row(k), &cache[i].basis_gs),
                            &cache[i].basis_gs_sq,
                        )
                        .unwrap();

                    vs.sub_mut(&mut basis_gs, &vs.scalar_mul(&cache[i].basis_gs, &mu_k_i));
                    debug_assert!(
                        self.ring()
                            .is_zero(&inner_product.inner_product(&basis_gs, &cache[i].basis_gs))
                    );
                    mu_k.push(mu_k_i);
                }
                let basis_gs_sq = inner_product.inner_product(&basis_gs, &basis_gs);
                cache.push(CacheEntry {
                    mu: mu_k,
                    basis_gs,
                    basis_gs_sq,
                });
                #[cfg(debug_assertions)]
                validate_cache(&basis, &cache);
            }

            // only reduce at k-1 before checking lovasz condition
            reduce(&mut cache, &mut h, &mut basis, k, k - 1);

            // check lovasz condition
            // ||b_k*||^2 >= (delta - mu_{k, k-1}^2) * ||b_{k-1}*||^2
            let mu_k_km1 = &cache[k].mu[k - 1];
            let bk_sq = &cache[k].basis_gs_sq;
            let bkm1_sq = &cache[k - 1].basis_gs_sq;
            let lovasz_condition = self
                .ring()
                .cmp(
                    bk_sq,
                    &self.ring().mul(
                        &self.ring().sub(
                            &self.ring().from_rat(delta),
                            &self.ring().mul(mu_k_km1, mu_k_km1),
                        ),
                        bkm1_sq,
                    ),
                )
                .is_ge();

            if lovasz_condition {
                // reduce all others only if lovasz condition passes
                for j in (0..(k - 1)).rev() {
                    reduce(&mut cache, &mut h, &mut basis, k, j);
                }

                // continue to next basis vector
                k += 1;
            } else {
                // swap b_{k} and b_{k-1}
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::Swap(k, k - 1),
                );
                // println!("{:?}", row_opp);
                row_opp.apply(&mut basis);
                row_opp.apply(&mut h);

                // update cache
                if k > 1 {
                    for j in 0..(k - 1) {
                        (cache[k].mu[j], cache[k - 1].mu[j]) =
                            (cache[k - 1].mu[j].clone(), cache[k].mu[j].clone());
                    }
                }
                let mu = cache[k].mu[k - 1].clone();
                let basis_gs_sq = self.ring().add(
                    &cache[k].basis_gs_sq,
                    &self
                        .ring()
                        .mul(&self.ring().mul(&mu, &mu), &cache[k - 1].basis_gs_sq),
                );
                cache[k].mu[k - 1] = self
                    .ring()
                    .try_divide(
                        &self.ring().mul(&mu, &cache[k - 1].basis_gs_sq),
                        &basis_gs_sq,
                    )
                    .unwrap();
                let b = cache[k - 1].basis_gs.clone();
                cache[k - 1].basis_gs = vs.add(&cache[k].basis_gs, &vs.scalar_mul(&b, &mu));
                cache[k].basis_gs = vs.sub(
                    &vs.scalar_mul(
                        &b,
                        &self
                            .ring()
                            .try_divide(&cache[k].basis_gs_sq, &basis_gs_sq)
                            .unwrap(),
                    ),
                    &vs.scalar_mul(&cache[k].basis_gs, &cache[k].mu[k - 1]),
                );
                cache[k].basis_gs_sq = self
                    .ring()
                    .try_divide(
                        &self
                            .ring()
                            .mul(&cache[k - 1].basis_gs_sq, &cache[k].basis_gs_sq),
                        &basis_gs_sq,
                    )
                    .unwrap();
                cache[k - 1].basis_gs_sq = basis_gs_sq;
                for i in (k + 1)..cache.len() {
                    let t = cache[i].mu[k].clone();
                    cache[i].mu[k] = self
                        .ring()
                        .sub(&cache[i].mu[k - 1], &self.ring().mul(&mu, &t));
                    cache[i].mu[k - 1] = self
                        .ring()
                        .add(&t, &self.ring().mul(&cache[k].mu[k - 1], &cache[i].mu[k]));
                }

                #[cfg(debug_assertions)]
                validate_cache(&basis, &cache);

                // now only up to k-1 are LLL reduced
                k -= 1;
                if k == 0 {
                    // though k=1 is fine since the first vector is always LLL reduced on it's own
                    k += 1;
                };
            }
        }

        #[cfg(debug_assertions)]
        assert!(self.equal(&self.mul(&h, &orig_basis).unwrap(), &basis));

        (h, basis)
    }

    /// Take the cols of M = `mat` as the basis for a lattice.
    ///
    /// Perform an LLL reduction on the lattice, returning a matricies (H, M*H) such that H is invertible and the cols of M*H are an LLL reduced basis for the lattice.
    ///
    /// `delta` must be in the range (-0.25, 1] with polynomial time complexity only for `delta` in (-0.25, 1). `delta` = 3/4 is the "default" value to use.
    pub fn lll_col_reduction_algorithm(
        &self,
        basis: Matrix<FS::Set>,
        inner_product: &impl RealInnerProduct<FS>,
        delta: &Rational,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (h, b) = self.lll_row_reduction_algorithm(basis.transpose(), inner_product, delta);
        (h.transpose(), b.transpose())
    }
}

impl<F: MetaType> Matrix<F>
where
    F::Signature: RealSubsetSignature
        + RealRoundingSignature
        + OrderedRingSignature
        + FieldSignature
        + ToStringSignature,
{
    pub fn lll_row_reduction_algorithm(
        self,
        inner_product: &impl RealInnerProduct<F::Signature>,
        delta: &Rational,
    ) -> (Matrix<F>, Matrix<F>) {
        Self::structure().lll_row_reduction_algorithm(self, inner_product, delta)
    }

    pub fn lll_col_reduction_algorithm(
        self,
        inner_product: &impl RealInnerProduct<F::Signature>,
        delta: &Rational,
    ) -> (Matrix<F>, Matrix<F>) {
        Self::structure().lll_col_reduction_algorithm(self, inner_product, delta)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    MatrixStructure<IntegerCanonicalStructure, B>
{
    /// Take the rows of M = `mat` as the basis for a lattice.
    ///
    /// Perform an LLL reduction on the lattice, returning matricies (H, H*M) such that H is invertible the rows of H*M are an LLL reduced basis for the lattice.
    ///
    /// `delta` must be in the range (-0.25, 1] with polynomial time complexity only for `delta` in (-0.25, 1). `delta` = 3/4 is the "default" value to use.
    pub fn lll_integral_row_reduction_algorithm(
        &self,
        mut basis: Matrix<Integer>,
        inner_product: &impl RealInnerProduct<IntegerCanonicalStructure>,
        delta: &Rational,
    ) -> (Matrix<Integer>, Matrix<Integer>) {
        let n = basis.rows(); // size of the basis
        let m = basis.cols(); // dimension of the ambient space
        assert!(m >= n);
        debug_assert!(self.rank(basis.clone()) == n);
        // 1/4 < delta <= 1
        assert!(Rational::ONE < Rational::from(4) * delta && delta <= &Rational::ONE);
        let vs = self.ring().free_module(m);

        if n == 0 {
            return (self.ident(0), basis);
        }

        let (delta_numerator, delta_denominator) = delta.numerator_and_denominator();
        debug_assert!(delta_numerator > Integer::ZERO);
        let delta_denominator = Integer::from(delta_denominator);

        #[cfg(debug_assertions)]
        let rational_inner_product = RealSymmetricInnerProduct::new(
            Rational::structure(),
            SymmetricMatrix::construct(m, |i, j| {
                Rational::from(
                    inner_product.inner_product(&vs.basis_element(i), &vs.basis_element(j)),
                )
            }),
        );

        #[cfg(debug_assertions)]
        let orig_basis = basis.clone();

        #[derive(Debug)]
        struct CacheEntry {
            // for i < k the ith entry is lambda_{k,i} = d_i mu_{k,i} = d_i <bj, bi>/<bi, bi>
            lambda: Vec<Integer>,
            // the determinant of the top k x k submatrix of the gram matrix of the basis
            d: Integer,
        }
        let mut cache = Vec::<CacheEntry>::with_capacity(n);

        #[cfg(debug_assertions)]
        let validate_cache = |basis: &Matrix<Integer>, cache: &Vec<CacheEntry>| {
            let true_mat_gs = basis
                .apply_map(|x| Rational::from(x))
                .gram_schmidt_row_orthogonalization(&rational_inner_product);
            let true_mu = Matrix::construct(n, n, |i, j| {
                Rational::structure()
                    .try_divide(
                        &rational_inner_product.inner_product(
                            &basis
                                .get_row(i)
                                .iter()
                                .map(Rational::from)
                                .collect::<Vec<_>>(),
                            &true_mat_gs.get_row(j),
                        ),
                        &rational_inner_product
                            .inner_product(&true_mat_gs.get_row(j), &true_mat_gs.get_row(j)),
                    )
                    .unwrap()
            });
            let true_gram_mat = Matrix::construct(n, n, |i, j| {
                inner_product.inner_product(&basis.get_row(i), &basis.get_row(j))
            });
            let true_d = (0..n)
                .map(|i| {
                    true_gram_mat
                        .submatrix((0..(i + 1)).collect(), (0..(i + 1)).collect())
                        .det()
                        .unwrap()
                })
                .collect::<Vec<_>>();
            let true_lambda = Matrix::construct(n, n, |i, j| {
                Integer::try_from(Rational::from(&true_d[j]) * true_mu.at(i, j).unwrap().clone())
                    .unwrap()
            });

            for (k, cache_entry) in cache.iter().enumerate() {
                assert_eq!(cache_entry.lambda.len(), k);
                for j in 0..k {
                    assert_eq!(&cache_entry.lambda[j], true_lambda.at(k, j).unwrap());
                }
                assert_eq!(cache_entry.d, true_d[k]);
            }
        };

        cache.push(CacheEntry {
            lambda: vec![],
            d: inner_product.inner_product(&basis.get_row(0), &basis.get_row(0)),
        });

        #[cfg(debug_assertions)]
        validate_cache(&basis, &cache);

        let mut h = self.ident(n);

        let reduce = |cache: &mut Vec<CacheEntry>,
                      h: &mut Matrix<Integer>,
                      basis: &mut Matrix<Integer>,
                      k: usize,
                      j: usize| {
            // The integer closest to lambda_{k, j}/d_j
            let q = Integer::structure()
                .quo(
                    &(Integer::TWO * &cache[k].lambda[j] + &cache[j].d),
                    &(Integer::TWO * &cache[j].d),
                )
                .unwrap();

            if q != Integer::ZERO {
                let neg_q = -&q;
                // update basis and basis transformation matrix
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::AddRowMul { i: k, j, x: neg_q },
                );
                println!("{:?}", row_opp);
                row_opp.apply(basis);
                row_opp.apply(h);

                // update the cache
                cache[k].lambda[j] = &cache[k].lambda[j] - &q * &cache[j].d;
                for i in 0..j {
                    cache[k].lambda[i] = &cache[k].lambda[i] - &q * &cache[j].lambda[i];
                }

                #[cfg(debug_assertions)]
                validate_cache(basis, cache);
            }
        };

        let mut k = 1;
        while k < n {
            // extend the cache if necessary
            if k >= cache.len() {
                debug_assert_eq!(k, cache.len());

                cache.push(CacheEntry {
                    lambda: vec![],   // to be filled in below
                    d: Integer::ZERO, // to be filled in below
                });
                for j in 0..(k + 1) {
                    let mut u = inner_product.inner_product(&basis.get_row(k), &basis.get_row(j));
                    for i in 0..j {
                        u = &cache[i].d * u - &cache[k].lambda[i] * &cache[j].lambda[i];
                        if i > 0 {
                            debug_assert!(u.divisible(&cache[i - 1].d));
                            u = u / &cache[i - 1].d;
                        }
                    }
                    if j < k {
                        cache[k].lambda.push(u);
                    } else {
                        debug_assert_eq!(j, k);
                        cache[k].d = u;
                    }
                }

                #[cfg(debug_assertions)]
                validate_cache(&basis, &cache);
                debug_assert_eq!(k + 1, cache.len());
            }

            // only reduce at k-1 before checking lovasz condition
            reduce(&mut cache, &mut h, &mut basis, k, k - 1);

            // The Lovasz condition is
            // d_k * d_{k-2} >= delta * d_{k-1}^2 - lambda_{k, k-1}^2
            let lambda_k_km1 = &cache[k].lambda[k - 1];
            let lovasz_condition = {
                if k == 1 {
                    &delta_denominator * &cache[1].d
                } else {
                    &delta_denominator * &cache[k].d * &cache[k - 2].d
                }
            }
            .cmp(
                &(&delta_numerator * &cache[k - 1].d * &cache[k - 1].d
                    - lambda_k_km1 * lambda_k_km1),
            )
            .is_ge();

            if lovasz_condition {
                // reduce all others only if lovasz condition passes
                for j in (0..(k - 1)).rev() {
                    reduce(&mut cache, &mut h, &mut basis, k, j);
                }

                // continue to next basis vector
                k += 1
            } else {
                // swap b_{k} and b_{k-1}
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring().clone(),
                    ElementaryOppType::Swap(k, k - 1),
                );
                println!("{:?}", row_opp);
                row_opp.apply(&mut basis);
                row_opp.apply(&mut h);

                // update the cache
                if k > 1 {
                    for j in 0..(k - 1) {
                        (cache[k].lambda[j], cache[k - 1].lambda[j]) =
                            (cache[k - 1].lambda[j].clone(), cache[k].lambda[j].clone());
                    }
                }
                let lambda = cache[k].lambda[k - 1].clone();
                let d_km1 = {
                    // compute x = (d_{k-2} * d_k + lambda^2) / d_{k-1}
                    let mut x = cache[k].d.clone();
                    if k > 1 {
                        x *= &cache[k - 2].d;
                    }
                    x += &lambda * &lambda;
                    debug_assert!(x.divisible(&cache[k - 1].d));
                    x = x / &cache[k - 1].d;
                    x
                };
                for i in (k + 1)..cache.len() {
                    let t = cache[i].lambda[k].clone();
                    cache[i].lambda[k] = {
                        // compute x = dk
                        let mut x = &cache[k].d * &cache[i].lambda[k - 1] - &lambda * &t;
                        debug_assert!(x.divisible(&cache[k - 1].d));
                        x = x / &cache[k - 1].d;
                        x
                    };
                    cache[i].lambda[k - 1] = {
                        // compute x = (d_km1 * t + lambda * lambda_{i, k}) / d_k
                        let mut x = &d_km1 * &t + &lambda * &cache[i].lambda[k];
                        debug_assert!(x.divisible(&cache[k].d));
                        x = x / &cache[k].d;
                        x
                    };
                }
                cache[k - 1].d = d_km1;

                #[cfg(debug_assertions)]
                validate_cache(&basis, &cache);

                // now only up to k-1 are LLL reduced
                k -= 1;
                if k == 0 {
                    // though k=1 is fine since the first vector is always LLL reduced on it's own
                    k += 1;
                };
            }
        }

        #[cfg(debug_assertions)]
        {
            use crate::matrix::RingMatricesSignature;

            let (rational_lll_h, _) = basis
                .apply_map(|x| Rational::from(x))
                .lll_row_reduction_algorithm(&rational_inner_product, delta);
            assert_eq!(
                rational_lll_h,
                Rational::structure().matrix_structure().ident(n)
            ); // since basis should now be LLL reduced

            assert!(self.equal(&self.mul(&h, &orig_basis).unwrap(), &basis));
        }

        (h, basis)
    }

    /// Take the cols of M = `mat` as the basis for a lattice.
    ///
    /// Perform an LLL reduction on the lattice, returning a matricies (H, M*H) such that H is invertible and the cols of M*H are an LLL reduced basis for the lattice.
    ///
    /// `delta` must be in the range (-0.25, 1] with polynomial time complexity only for `delta` in (-0.25, 1). `delta` = 3/4 is the "default" value to use.
    pub fn lll_integral_col_reduction_algorithm(
        &self,
        basis: Matrix<Integer>,
        inner_product: &impl RealInnerProduct<IntegerCanonicalStructure>,
        delta: &Rational,
    ) -> (Matrix<Integer>, Matrix<Integer>) {
        let (h, b) =
            self.lll_integral_row_reduction_algorithm(basis.transpose(), inner_product, delta);
        (h.transpose(), b.transpose())
    }
}

impl Matrix<Integer> {
    pub fn lll_integral_row_reduction_algorithm(
        self,
        inner_product: &impl RealInnerProduct<IntegerCanonicalStructure>,
        delta: &Rational,
    ) -> (Matrix<Integer>, Matrix<Integer>) {
        Self::structure().lll_integral_row_reduction_algorithm(self, inner_product, delta)
    }

    pub fn lll_integral_col_reduction_algorithm(
        self,
        inner_product: &impl RealInnerProduct<IntegerCanonicalStructure>,
        delta: &Rational,
    ) -> (Matrix<Integer>, Matrix<Integer>) {
        Self::structure().lll_integral_col_reduction_algorithm(self, inner_product, delta)
    }
}

#[cfg(test)]
mod tests {
    // use crate::matrix::{
    //     Matrix, SymmetricMatrix,
    //     lll_reduction::{integral_lll_reduction, is_positive_definite},
    // };
    // use algebraeon_nzq::Integer;

    // #[test]
    // fn test_is_positive_definite() {
    //     let mat = Matrix::from_rows(vec![
    //         vec![Integer::from(2), Integer::from(-1), Integer::from(0)],
    //         vec![Integer::from(-1), Integer::from(2), Integer::from(-1)],
    //         vec![Integer::from(0), Integer::from(-1), Integer::from(0)],
    //     ]);
    //     assert!(is_positive_definite(&mat.try_into().unwrap()));
    // }

    use std::str::FromStr;

    use crate::matrix::{Matrix, StandardInnerProduct};
    use algebraeon_nzq::{Integer, Rational};
    use algebraeon_sets::structure::MetaType;

    #[test]
    fn approx_golden_ratio() {
        let m = Matrix::<Rational>::from_rows(vec![
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("10000").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("16180").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("26180").unwrap(),
            ],
        ]);

        m.pprint();

        assert_eq!(
            m.clone()
                .lll_row_reduction_algorithm(
                    &StandardInnerProduct::new(Rational::structure()),
                    &Rational::from_str("3/4").unwrap(),
                )
                .0
                .get_row(0),
            vec![
                Rational::from_str("-1").unwrap(),
                Rational::from_str("-1").unwrap(),
                Rational::from_str("1").unwrap(),
            ],
        );
    }

    #[test]
    fn approx_degree_5_root() {
        // Find the polynomial x^5-3x^4-2x^3+x^2+7x-3 of which ~1.1614471390 is a root

        let m = Matrix::<Rational>::from_rows(vec![
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("10000000000").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("11614471390").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("13489594567").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("15667451017").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("18196916159").unwrap(),
            ],
            vec![
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("0").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("21134756212").unwrap(),
            ],
        ]);

        m.pprint();

        assert_eq!(
            m.clone()
                .lll_row_reduction_algorithm(
                    &StandardInnerProduct::new(Rational::structure()),
                    &Rational::from_str("3/4").unwrap(),
                )
                .0
                .get_row(0),
            vec![
                Rational::from_str("-3").unwrap(),
                Rational::from_str("7").unwrap(),
                Rational::from_str("1").unwrap(),
                Rational::from_str("-2").unwrap(),
                Rational::from_str("-3").unwrap(),
                Rational::from_str("1").unwrap(),
            ],
        );
    }

    #[test]
    fn approx_golden_ratio_integral() {
        let m = Matrix::<Integer>::from_rows(vec![
            vec![
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("10000").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("16180").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("26180").unwrap(),
            ],
        ]);

        m.pprint();

        assert_eq!(
            m.clone()
                .lll_integral_row_reduction_algorithm(
                    &StandardInnerProduct::new(Integer::structure()),
                    &Rational::from_str("3/4").unwrap(),
                )
                .0
                .get_row(0),
            vec![
                Integer::from_str("-1").unwrap(),
                Integer::from_str("-1").unwrap(),
                Integer::from_str("1").unwrap(),
            ],
        );
    }

    #[test]
    fn approx_degree_5_root_integral() {
        // Find the polynomial x^5-3x^4-2x^3+x^2+7x-3 of which ~1.1614471390 is a root

        let m = Matrix::<Integer>::from_rows(vec![
            vec![
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("10000000000").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("11614471390").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("13489594567").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("15667451017").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("18196916159").unwrap(),
            ],
            vec![
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("0").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("21134756212").unwrap(),
            ],
        ]);

        m.pprint();

        assert_eq!(
            m.clone()
                .lll_integral_row_reduction_algorithm(
                    &StandardInnerProduct::new(Integer::structure()),
                    &Rational::from_str("3/4").unwrap(),
                )
                .0
                .get_row(0),
            vec![
                Integer::from_str("-3").unwrap(),
                Integer::from_str("7").unwrap(),
                Integer::from_str("1").unwrap(),
                Integer::from_str("-2").unwrap(),
                Integer::from_str("-3").unwrap(),
                Integer::from_str("1").unwrap(),
            ],
        );
    }
}
