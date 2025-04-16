use super::*;

#[derive(Debug, Clone)]
pub struct JordanBlock<FS: AlgebraicClosureStructure>
where
    PolynomialStructure<FS::BFS>:
        FactorableStructure + SetStructure<Set = Polynomial<<FS::BFS as SetStructure>::Set>>,
{
    eigenvalue: FS::Set,
    blocksize: usize,
}

impl<FS: AlgebraicClosureStructure> JordanBlock<FS>
where
    PolynomialStructure<FS::BFS>:
        FactorableStructure + SetStructure<Set = Polynomial<<FS::BFS as SetStructure>::Set>>,
{
    pub fn matrix(&self, field: &FS) -> Matrix<FS::Set> {
        // let base_field = field.base_field();
        Matrix::construct(self.blocksize, self.blocksize, |r, c| {
            if r == c {
                self.eigenvalue.clone()
            } else if r + 1 == c {
                field.one()
            } else {
                field.zero()
            }
        })
    }
}

#[derive(Debug, Clone)]
pub struct JordanNormalForm<FS: AlgebraicClosureStructure>
where
    PolynomialStructure<FS::BFS>:
        FactorableStructure + SetStructure<Set = Polynomial<<FS::BFS as SetStructure>::Set>>,
{
    field: FS,
    blocks: Vec<JordanBlock<FS>>,
}

impl<FS: AlgebraicClosureStructure> JordanNormalForm<FS>
where
    PolynomialStructure<FS::BFS>:
        FactorableStructure + SetStructure<Set = Polynomial<<FS::BFS as SetStructure>::Set>>,
{
    pub fn matrix(&self) -> Matrix<FS::Set> {
        let ac_mat_structure = MatrixStructure::new(self.field.clone());
        ac_mat_structure.join_diag(
            self.blocks
                .iter()
                .map(|block| block.matrix(&self.field))
                .collect(),
        )
    }
}

impl<FS: AlgebraicClosureStructure> MatrixStructure<FS>
where
    PolynomialStructure<FS::BFS>:
        FactorableStructure + SetStructure<Set = Polynomial<<FS::BFS as SetStructure>::Set>>,
{
    pub fn eigenvalues_list(&self, mat: Matrix<<FS::BFS as SetStructure>::Set>) -> Vec<FS::Set> {
        let base_field_mat_structure = MatrixStructure::new(self.ring().base_field().clone());
        self.ring()
            .all_roots_list(
                &base_field_mat_structure
                    .characteristic_polynomial(mat)
                    .unwrap(),
            )
            .unwrap()
    }

    pub fn eigenvalues_unique(&self, mat: Matrix<<FS::BFS as SetStructure>::Set>) -> Vec<FS::Set> {
        let base_field_mat_structure = MatrixStructure::new(self.ring().base_field().clone());
        self.ring()
            .all_roots_unique(
                &base_field_mat_structure
                    .characteristic_polynomial(mat)
                    .unwrap(),
            )
            .unwrap()
    }

    pub fn eigenvalues_powers(
        &self,
        mat: Matrix<<FS::BFS as SetStructure>::Set>,
    ) -> Vec<(FS::Set, usize)> {
        let base_field_mat_structure = MatrixStructure::new(self.ring().base_field().clone());
        self.ring()
            .all_roots_powers(
                &base_field_mat_structure
                    .characteristic_polynomial(mat)
                    .unwrap(),
            )
            .unwrap()
    }

    pub fn generalized_col_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
        eigenvalue: &FS::Set,
        k: usize,
    ) -> LinearSubspace<FS::Set> {
        let n = mat.rows();
        assert_eq!(n, mat.cols());
        //compute ker((M - xI)^k)
        self.col_kernel(
            self.nat_pow(
                &self
                    .add(
                        &mat.apply_map(|x| self.ring().base_field_inclusion(x)),
                        &self.neg(self.mul_scalar(self.ident(n), eigenvalue)),
                    )
                    .unwrap(),
                &Natural::from(k),
            )
            .unwrap(),
        )
    }

    pub fn generalized_row_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
        eigenvalue: &FS::Set,
        k: usize,
    ) -> LinearSubspace<FS::Set> {
        LinearSubspaceStructure::new(self.ring().clone())
            .transpose(&self.generalized_col_eigenspace(&mat.transpose_ref(), eigenvalue, k))
    }

    pub fn col_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
        eigenvalue: &FS::Set,
    ) -> LinearSubspace<FS::Set> {
        self.generalized_col_eigenspace(mat, eigenvalue, 1)
    }

    pub fn row_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
        eigenvalue: &FS::Set,
    ) -> LinearSubspace<FS::Set> {
        self.generalized_row_eigenspace(mat, eigenvalue, 1)
    }

    //return the jordan normal form F of the matrix M and a basis matrix B such that
    // B^-1 M B = J
    pub fn jordan_algorithm(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
    ) -> (JordanNormalForm<FS>, Matrix<FS::Set>) {
        let n = mat.rows();
        assert_eq!(n, mat.cols());

        let ac_field = self.ring();
        let ac_mat_structure = self;
        let ac_linlat_structure = LinearSubspaceStructure::new(ac_field.clone());

        let ac_mat = mat.apply_map(|x| self.ring().base_field_inclusion(x));

        let mut basis = vec![];
        let mut eigenvalues = vec![]; //store (gesp_basis, eigenvalue, multiplicity)
        for (eigenvalue, multiplicity) in self.eigenvalues_powers(mat.clone()) {
            let eigenspace = self.generalized_col_eigenspace(mat, &eigenvalue, multiplicity);
            debug_assert_eq!(ac_linlat_structure.rank(&eigenspace), multiplicity);
            basis.append(&mut ac_linlat_structure.basis_matrices(&eigenspace));
            eigenvalues.push((eigenvalue, multiplicity));
        }

        //b = direct sum of generalized eigenspace
        let gesp_basis = Matrix::join_cols(mat.rows(), basis);
        //b^-1 * mat * b = block diagonal of generalized eigenspaces
        let gesp_blocks_mat = ac_mat_structure
            .mul(
                &ac_mat_structure.inv(gesp_basis.clone()).unwrap(),
                &ac_mat_structure.mul(&ac_mat, &gesp_basis).unwrap(),
            )
            .unwrap();

        let mut idx_to_block = vec![];
        for (b, (_eval, mult)) in eigenvalues.iter().enumerate() {
            for _i in 0..*mult {
                idx_to_block.push(b);
            }
        }

        // println!("{:?}", idx_to_block);
        // println!("{:?}", eigenvalues);
        // ac_mat_structure.pprint(&gesp_blocks_mat);

        //extract the blocks from the block diagonal gesp_blocks_mat
        let mut gesp_blocks = vec![];
        let mut cum_mult = 0;
        for (eval, mult) in eigenvalues {
            gesp_blocks.push((
                eval,
                mult,
                Matrix::construct(mult, mult, |r, c| {
                    gesp_blocks_mat
                        .at(cum_mult + r, cum_mult + c)
                        .unwrap()
                        .clone()
                }),
            ));
            cum_mult += mult;
        }
        debug_assert_eq!(cum_mult, n);
        drop(gesp_blocks_mat);

        // Vec<(eval, multiplicity, Vec<Jordan Block>)>
        let jnf_info = gesp_blocks
            .into_iter()
            .map(|(eval, m, mat_t)| {
                debug_assert_eq!(mat_t.rows(), m);
                debug_assert_eq!(mat_t.cols(), m);
                // println!("eval = {:?} m={}", eval, m);
                // ac_mat_structure.pprint(&mat_t);
                //all eigenvalues of T are eval
                //let S = T - x I so that all eigenvlues of S are zero
                let mat_s = ac_mat_structure
                    .add(
                        &mat_t,
                        &ac_mat_structure
                            .mul_scalar(ac_mat_structure.ident(m), &ac_field.neg(&eval)),
                    )
                    .unwrap();

                let jb_basis = {
                    debug_assert!(m >= 1);

                    let mut mat_s_pows = vec![ac_mat_structure.ident(m), mat_s.clone()];
                    for _i in 0..(m - 1) {
                        mat_s_pows.push(
                            ac_mat_structure
                                .mul(mat_s_pows.last().unwrap(), &mat_s)
                                .unwrap(),
                        );
                    }
                    debug_assert!(
                        ac_mat_structure.equal(&ac_mat_structure.zero(m, m), &mat_s_pows[m])
                    );
                    // for (i, spow) in mat_s_pows.iter().enumerate() {
                    //     println!("s^{}", i);
                    //     ac_mat_structure.pprint(&spow);
                    // }
                    let mat_s_pow_kers = mat_s_pows
                        .into_iter()
                        .map(|s_mat_pow| ac_mat_structure.col_kernel(s_mat_pow))
                        .collect_vec();
                    // ker(S) in ker(S^2) in ker(S^3) in ...
                    // for (i, ker) in mat_s_pow_kers.iter().enumerate() {
                    //     println!("ker(s^{})", i);
                    //     ac_linlat_structure.pprint(&ker);
                    // }

                    let mut accounted = ac_linlat_structure.zero(m, 1);
                    let mut jordan_block_bases = vec![];
                    for k in (0..m).rev() {
                        //extend the basis by stuff in ker(S^{k+1}) but not in ker(S^k) and their images under S, and which are not already accounted for
                        // println!("k = {} {}", k + 1, k);
                        let ker_ext = ac_linlat_structure.from_basis(
                            m,
                            1,
                            ac_linlat_structure
                                .extension_basis(&mat_s_pow_kers[k], &mat_s_pow_kers[k + 1]),
                        );

                        let unaccounted_ker_ext_basis = ac_linlat_structure.extension_basis(
                            &ac_linlat_structure.intersect_pair(m, 1, &accounted, &ker_ext),
                            &ker_ext,
                        );
                        // let unaccounted_ker_ext =
                        //     ac_linlat_structure.from_basis(m, 1, unaccounted_ker_ext_basis.clone());

                        for ukeb in unaccounted_ker_ext_basis {
                            //one new jordan block for each ukeb
                            // println!("ukeb");
                            // ac_mat_structure.pprint(&ukeb);

                            let mut jb_basis = vec![ukeb];
                            for _i in 0..k {
                                let ukeb_img = ac_mat_structure
                                    .mul(&mat_s, &jb_basis.last().unwrap())
                                    .unwrap();
                                // println!("ukeb_img #{}", i);
                                // ac_mat_structure.pprint(&ukeb_img);
                                jb_basis.push(ukeb_img);
                            }

                            accounted = ac_linlat_structure.sum_pair(
                                m,
                                1,
                                &accounted,
                                &ac_linlat_structure.from_basis(m, 1, jb_basis.clone()),
                            );

                            jordan_block_bases.push(jb_basis.into_iter().rev().collect_vec());
                        }
                    }

                    // println!(
                    //     "jb sizes = {:?}",
                    //     jordan_block_bases.iter().map(|v| v.len()).collect_vec()
                    // );

                    jordan_block_bases

                    // Matrix::join_cols(m, jordan_block_bases.into_iter().flatten().collect_vec())
                };

                // println!("gesp_jordan_basis");
                // ac_mat_structure.pprint(&jb_basis);
                // ac_mat_structure.pprint(
                //     &ac_mat_structure
                //         .mul(
                //             &ac_mat_structure.inv(jb_basis.clone()).unwrap(),
                //             &ac_mat_structure.mul(&mat_t, &jb_basis).unwrap(),
                //         )
                //         .unwrap(),
                // );

                (eval, m, jb_basis)
            })
            .collect_vec();

        let mut jordan_blocks = vec![];
        let mut jnf_basis_rel_gesp_basis: Vec<Matrix<FS::Set>> = vec![];
        for (eval, mult, blocks) in jnf_info {
            // println!("eval={:?}, mult={}", eval, mult);
            let mut eigenblock_basis = vec![];
            for mut block in blocks {
                jordan_blocks.push(JordanBlock {
                    eigenvalue: eval.clone(),
                    blocksize: block.len(),
                });
                eigenblock_basis.append(&mut block);
            }
            jnf_basis_rel_gesp_basis.push(Matrix::join_cols(mult, eigenblock_basis));
        }
        let jnf = JordanNormalForm {
            field: self.ring().clone(),
            blocks: jordan_blocks,
        };
        let jordan_blocks_basis = ac_mat_structure.join_diag(jnf_basis_rel_gesp_basis);

        // ac_mat_structure.pprint(&jnf.matrix());

        // println!("jordan_blocks_basis");
        // ac_mat_structure.pprint(&jordan_blocks_basis);

        let jnf_basis = ac_mat_structure
            .mul(&gesp_basis, &jordan_blocks_basis)
            .unwrap();
        // println!("jnf_basis");
        // ac_mat_structure.pprint(&jnf_basis);

        //check that B^-1 M B = JNF
        debug_assert!(
            ac_mat_structure.equal(
                &ac_mat_structure
                    .mul(
                        &ac_mat_structure.inv(jnf_basis.clone()).unwrap(),
                        &ac_mat_structure.mul(&ac_mat, &jnf_basis).unwrap(),
                    )
                    .unwrap(),
                &jnf.matrix()
            )
        );
        // println!("jnf");
        // ac_mat_structure.pprint(&jnf);

        // todo!()

        (jnf, jnf_basis)
    }

    pub fn jordan_normal_form(
        &self,
        mat: &Matrix<<FS::BFS as SetStructure>::Set>,
    ) -> Matrix<FS::Set> {
        self.jordan_algorithm(mat).0.matrix()
    }

    //TODO: find basis which make two matricies similar if one exists
    /*
    def similar_basis(self, other):
    #find a basis in which self looks like other
    #equivelently, find P such that P^-1*self*P == other
    if type(self) == type(other) == Matrix:
        if self.n == other.n:
            if self.jcf_spec() == other.jcf_spec():
                #need to find a jcf basis for self and other such that the jcf matricies are identical (including order (thats the only hard part))
                #NOTE - by the implementation of the algorithm used, each eigen block will be consistently ordered - largest first
                #HOWEVER, the order of the eigen block is still unknown (and an order cant be imposed in the algorhtm becasue in general, arbitrary number things cant be ordered in a consistent way)
                self_jcf_info = self.jcf_info()
                other_jcf_info = other.jcf_info()
                #rewrite these in terms of {e_val : info}
                self_jcf_info = {info["ev"] : info for info in self_jcf_info}
                other_jcf_info = {info["ev"] : info for info in other_jcf_info}
                assert self_jcf_info.keys() == other_jcf_info.keys()
                keys = list(self_jcf_info.keys()) #decide a consistent order here
                #reorder the info
                self_jcf_info = [self_jcf_info[ev] for ev in keys]
                other_jcf_info = [other_jcf_info[ev] for ev in keys]
                #now both info lists have the eigen values in the same order
                #as well as all blocks within each eigen block being in the right order
                self_jcf_basis = Matrix.jcf_info_to_jcf_basis(self_jcf_info)
                other_jcf_basis = Matrix.jcf_info_to_jcf_basis(other_jcf_info)
                return self_jcf_basis * other_jcf_basis ** -1
            else:
                raise Exception("Matricies are not similar so cant find a basis in which one looks like the other")
    raise NotImplementedError
    */
}

#[cfg(test)]
mod tests {
    use crate::number::algebraic::complex::ComplexAlgebraic;

    use super::*;

    #[test]
    fn jordan_normal_form() {
        let mat = Matrix::from_rows(vec![
            vec![
                Rational::from(3),
                Rational::from(1),
                Rational::from(0),
                Rational::from(1),
            ],
            vec![
                Rational::from(-1),
                Rational::from(5),
                Rational::from(4),
                Rational::from(1),
            ],
            vec![
                Rational::from(0),
                Rational::from(0),
                Rational::from(2),
                Rational::from(0),
            ],
            vec![
                Rational::from(0),
                Rational::from(0),
                Rational::from(0),
                Rational::from(4),
            ],
        ]);

        mat.pprint();
        for root in
            MatrixStructure::new(ComplexAlgebraic::structure()).eigenvalues_list(mat.clone())
        {
            println!("{}", root);
        }

        let (j, b) = MatrixStructure::new(ComplexAlgebraic::structure()).jordan_algorithm(&mat);
        println!("{:?}", j);
        j.matrix().pprint();
        b.pprint();

        let mat = Matrix::<Rational>::from_rows(vec![vec![1, 0, 0], vec![0, 0, -1], vec![0, 2, 0]]);

        mat.pprint();
        for root in
            MatrixStructure::new(ComplexAlgebraic::structure()).eigenvalues_list(mat.clone())
        {
            println!("{}", root);
        }

        let (j, b) = MatrixStructure::new(ComplexAlgebraic::structure()).jordan_algorithm(&mat);
        println!("{:?}", j);
        j.matrix().pprint();
        b.pprint();
    }
}
