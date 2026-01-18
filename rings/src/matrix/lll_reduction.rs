use crate::structure::{AdditionSignature, AdditiveGroupSignature, SemiModuleSignature};
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
use algebraeon_nzq::{Integer, Rational};
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
    /// Perform an LLL reduction on the lattice, returning a matrix H such that the rows of H*M are an LLL reduced basis for the lattice.
    ///
    /// `delta` must be in the range (-0.25, 1] with polynomial time complexity only for `delta` in (-0.25, 1). `delta` = 3/4 is the "default" value to use.
    pub fn lll_row_reduction_algorithm(
        &self,
        mut basis: Matrix<FS::Set>,
        inner_product: &impl RealInnerProduct<FS>,
        delta: &Rational,
    ) -> Matrix<FS::Set> {
        let n = basis.rows(); // size of the basis
        let m = basis.cols(); // dimension of the ambient space
        assert!(m >= n);
        debug_assert!(self.rank(basis.clone()) == n);
        // 1/4 < delta <= 1
        assert!(Rational::ONE < Rational::from(4) * delta && delta <= &Rational::ONE);
        let vs = self.ring().free_module(m);

        if n == 0 {
            return self.ident(0);
        }

        // variables to hold cached values
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

            for j in 0..k {
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
                    row_opp.apply(&mut basis);
                    row_opp.apply(&mut h);

                    // update cache
                    self.ring().sub_mut(&mut cache[k].mu[j], &mu_rounded);
                    for i in 0..j {
                        cache[k].mu[i] = self.ring().sub(
                            &cache[k].mu[i],
                            &self.ring().mul(&mu_rounded, &cache[j].mu[i]),
                        );
                    }

                    #[cfg(debug_assertions)]
                    validate_cache(&basis, &cache);
                }
            }

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

        h
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
    ) -> Matrix<F> {
        Self::structure().lll_row_reduction_algorithm(self, inner_product, delta)
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
    use algebraeon_nzq::Rational;
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
            m.clone().lll_row_reduction_algorithm(
                &StandardInnerProduct::new(Rational::structure()),
                &Rational::from_str("3/4").unwrap(),
            ),
            Matrix::<Rational>::from_rows(vec![
                vec![
                    Rational::from_str("-1").unwrap(),
                    Rational::from_str("-1").unwrap(),
                    Rational::from_str("1").unwrap(),
                ],
                vec![
                    Rational::from_str("-48").unwrap(),
                    Rational::from_str("41").unwrap(),
                    Rational::from_str("-7").unwrap(),
                ],
                vec![
                    Rational::from_str("78").unwrap(),
                    Rational::from_str("-66").unwrap(),
                    Rational::from_str("11").unwrap(),
                ],
            ])
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
            m.clone().lll_row_reduction_algorithm(
                &StandardInnerProduct::new(Rational::structure()),
                &Rational::from_str("3/4").unwrap(),
            ),
            Matrix::<Rational>::from_rows(vec![
                vec![
                    Rational::from_str("-3").unwrap(),
                    Rational::from_str("7").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("-2").unwrap(),
                    Rational::from_str("-3").unwrap(),
                    Rational::from_str("1").unwrap(),
                ],
                vec![
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("-14").unwrap(),
                    Rational::from_str("45").unwrap(),
                    Rational::from_str("4").unwrap(),
                    Rational::from_str("-47").unwrap(),
                    Rational::from_str("16").unwrap(),
                ],
                vec![
                    Rational::from_str("35").unwrap(),
                    Rational::from_str("-9").unwrap(),
                    Rational::from_str("9").unwrap(),
                    Rational::from_str("-48").unwrap(),
                    Rational::from_str("-16").unwrap(),
                    Rational::from_str("32").unwrap(),
                ],
                vec![
                    Rational::from_str("-61").unwrap(),
                    Rational::from_str("0").unwrap(),
                    Rational::from_str("21").unwrap(),
                    Rational::from_str("-25").unwrap(),
                    Rational::from_str("72").unwrap(),
                    Rational::from_str("-28").unwrap(),
                ],
                vec![
                    Rational::from_str("-2").unwrap(),
                    Rational::from_str("56").unwrap(),
                    Rational::from_str("-43").unwrap(),
                    Rational::from_str("5").unwrap(),
                    Rational::from_str("51").unwrap(),
                    Rational::from_str("-50").unwrap(),
                ],
                vec![
                    Rational::from_str("64").unwrap(),
                    Rational::from_str("7").unwrap(),
                    Rational::from_str("15").unwrap(),
                    Rational::from_str("1").unwrap(),
                    Rational::from_str("-4").unwrap(),
                    Rational::from_str("-41").unwrap(),
                ],
            ])
        );
    }
}
