use crate::linear::{
    finitely_free_modules::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
};

use super::*;

#[derive(Debug, Clone)]
pub struct JordanBlock<FS: AlgebraicClosureSignature>
where
    PolynomialStructure<FS::BFS, FS::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<FS::BFS as SetSignature>::Set>>,
{
    eigenvalue: FS::Set,
    blocksize: usize,
}

impl<FS: AlgebraicClosureSignature> JordanBlock<FS>
where
    PolynomialStructure<FS::BFS, FS::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<FS::BFS as SetSignature>::Set>>,
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
pub struct JordanNormalForm<FS: AlgebraicClosureSignature>
where
    PolynomialStructure<FS::BFS, FS::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<FS::BFS as SetSignature>::Set>>,
{
    field: FS,
    blocks: Vec<JordanBlock<FS>>,
}

impl<FS: AlgebraicClosureSignature> JordanNormalForm<FS>
where
    PolynomialStructure<FS::BFS, FS::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<FS::BFS as SetSignature>::Set>>,
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

impl<FS: AlgebraicClosureSignature, FSB: BorrowedStructure<FS>> MatrixStructure<FS, FSB>
where
    PolynomialStructure<FS::BFS, FS::BFS>:
        FactorableSignature + SetSignature<Set = Polynomial<<FS::BFS as SetSignature>::Set>>,
{
    pub fn eigenvalues_list(&self, mat: Matrix<<FS::BFS as SetSignature>::Set>) -> Vec<FS::Set> {
        let base_field_mat_structure = MatrixStructure::new(self.ring().base_field().clone());
        self.ring()
            .all_roots_list(
                &base_field_mat_structure
                    .characteristic_polynomial(mat)
                    .unwrap(),
            )
            .unwrap()
    }

    pub fn eigenvalues_unique(&self, mat: Matrix<<FS::BFS as SetSignature>::Set>) -> Vec<FS::Set> {
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
        mat: Matrix<<FS::BFS as SetSignature>::Set>,
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
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
        eigenvalue: &FS::Set,
        k: usize,
    ) -> FinitelyFreeSubmodule<FS::Set> {
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
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
        eigenvalue: &FS::Set,
        k: usize,
    ) -> FinitelyFreeSubmodule<FS::Set> {
        self.generalized_col_eigenspace(&mat.transpose_ref(), eigenvalue, k)
    }

    pub fn col_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
        eigenvalue: &FS::Set,
    ) -> FinitelyFreeSubmodule<FS::Set> {
        self.generalized_col_eigenspace(mat, eigenvalue, 1)
    }

    pub fn row_eigenspace(
        &self,
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
        eigenvalue: &FS::Set,
    ) -> FinitelyFreeSubmodule<FS::Set> {
        self.generalized_row_eigenspace(mat, eigenvalue, 1)
    }

    //return the jordan normal form F of the matrix M and a basis matrix B such that
    // B^-1 M B = J
    pub fn jordan_algorithm(
        &self,
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
    ) -> (JordanNormalForm<FS>, Matrix<FS::Set>) {
        let n = mat.rows();
        assert_eq!(n, mat.cols());

        let ac_field = self.ring();
        let ac_mat_structure = self;

        let ac_mat = mat.apply_map(|x| self.ring().base_field_inclusion(x));

        let mut basis = vec![];
        let mut eigenvalues = vec![]; //store (gesp_basis, eigenvalue, multiplicity)
        for (eigenvalue, multiplicity) in self.eigenvalues_powers(mat.clone()) {
            let eigenspace = self.generalized_col_eigenspace(mat, &eigenvalue, multiplicity);
            debug_assert_eq!(eigenspace.rank(), multiplicity);
            basis.append(&mut eigenspace.basis());
            eigenvalues.push((eigenvalue, multiplicity));
        }

        //b = direct sum of generalized eigenspace
        let gesp_basis = Matrix::construct(mat.rows(), basis.len(), |r, c| basis[c][r].clone());
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

                    let mat_s_pow_kers = mat_s_pows
                        .into_iter()
                        .map(|s_mat_pow| ac_mat_structure.col_kernel(s_mat_pow))
                        .collect_vec();
                    // ker(S) in ker(S^2) in ker(S^3) in ...

                    let module = FinitelyFreeModuleStructure::<FS, _>::new(ac_field, m);
                    let mut accounted = module
                        .submodules()
                        .zero_submodule();

                    let mut jordan_block_bases = vec![];
                    for k in (0..m).rev() {
                        //extend the basis by stuff in ker(S^{k+1}) but not in ker(S^k) and their images under S, and which are not already accounted for
                        let ker_ext = module.submodules().from_span(
                            module.submodules().extension_basis(
                                &mat_s_pow_kers[k],
                                &mat_s_pow_kers[k + 1],
                            )
                            .iter()
                            .collect(),
                        );

                        let unaccounted_ker_ext_basis = module.submodules().extension_basis(
                            &module.submodules().intersect(&accounted, &ker_ext),
                            &ker_ext,
                        );

                        for ukeb in unaccounted_ker_ext_basis {
                            //one new jordan block for each ukeb

                            let mut jb_basis = vec![ukeb];
                            for _i in 0..k {
                                let ukeb_img = ac_mat_structure
                                    .mul(
                                        &mat_s,
                                        &Matrix::construct(m, 1, |r, _| {
                                            jb_basis.last().unwrap()[r].clone()
                                        }),
                                    )
                                    .unwrap();
                                assert_eq!(ukeb_img.cols(), 1);
                                assert_eq!(ukeb_img.rows(), m);
                                // ac_mat_structure.pprint(&ukeb_img);
                                jb_basis.push(ukeb_img.get_col(0));
                            }

                            accounted = module.submodules().add(
                                &accounted,
                                &module.submodules().from_span(
                                    jb_basis.iter().collect(),
                                ),
                            );
                            jordan_block_bases.push(jb_basis.into_iter().rev().collect_vec());
                        }
                    }

                    jordan_block_bases
                };

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
            jnf_basis_rel_gesp_basis.push(Matrix::join_cols(
                mult,
                eigenblock_basis
                    .into_iter()
                    .map(|col| Matrix::from_cols(vec![col]))
                    .collect(),
            ));
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
        mat: &Matrix<<FS::BFS as SetSignature>::Set>,
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
    use super::*;
    use crate::rings::isolated_algebraic::complex::ComplexAlgebraic;

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
