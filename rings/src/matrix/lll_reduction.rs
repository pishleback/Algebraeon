use crate::{
    isolated_algebraic::polynomial::Interleave,
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    matrix::{Matrix, SymmetricMatrix},
};
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::{MetaType, SetSignature};

fn is_positive_definite(mat: &SymmetricMatrix<Integer>) -> bool {
    let n = mat.n();
    let mat = Matrix::construct(n, n, |r, c| mat.get(r, c).unwrap().clone());
    for i in 1..n {
        if mat
            .submatrix((0..i).collect(), (0..i).collect())
            .det()
            .unwrap()
            <= Integer::ZERO
        {
            return false;
        }
    }
    true
}

fn integral_lll_reduction(
    row_basis: Matrix<Integer>,
    gram_matrix: &SymmetricMatrix<Integer>,
) -> Result<(), String> {
    // source: Cohen H - a course in computational algebraic number theory
    debug_assert!(is_positive_definite(gram_matrix));
    let n = row_basis.rows(); // number of vectors in the basis
    let m = row_basis.cols(); // dimension of the ambient space
    assert_eq!(m, gram_matrix.n());
    debug_assert_eq!(row_basis.rank(), n);

    // keep track of the transformatinos performed on the lattice basis
    // the columns express the linear combination of the input basis required to obtain the output basis
    let mut h = Matrix::<Integer>::ident(n);
    let mut k = 1;
    let mut k_max = 0;
    let mut d = vec![Integer::ZERO; n];
    d[0] = todo!("row_basis.get_row(0) dot row_basis.get_row(0)");
    let mut lambda = todo!();

    loop {
        if k > k_max {
            k_max = k;
            for j in 0..(k + 1) {
                let mut u: Integer = todo!("row_basis.get_row(k) dot row_basis.get_row(j)");
                for i in 0..j {
                    u = (d[i] * u - lambda[k][i] * lambda[j][i])
                        .try_divide(d[i - 1])
                        .unwrap();
                }
            }
        }
    }

    todo!()
}

#[cfg(test)]
mod tests {
    use crate::matrix::{
        Matrix, SymmetricMatrix,
        lll_reduction::{integral_lll_reduction, is_positive_definite},
    };
    use algebraeon_nzq::Integer;

    #[test]
    fn test_is_positive_definite() {
        let mat = Matrix::from_rows(vec![
            vec![Integer::from(2), Integer::from(-1), Integer::from(0)],
            vec![Integer::from(-1), Integer::from(2), Integer::from(-1)],
            vec![Integer::from(0), Integer::from(-1), Integer::from(0)],
        ]);
        assert!(is_positive_definite(&mat.try_into().unwrap()));
    }

    #[test]
    fn test() {
        let gram_matrix =
            SymmetricMatrix::construct(3, |r, c| if r == c { Integer::ONE } else { Integer::ZERO });

        let h = Matrix::from_rows(vec![
            vec![Integer::from(100), Integer::from(69), Integer::from(0)],
            vec![Integer::from(101), Integer::from(420), Integer::from(0)],
        ]);

        integral_lll_reduction(h, &gram_matrix).unwrap();

        println!("Hii");
    }
}
