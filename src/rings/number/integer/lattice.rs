use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::integer::Integer;
use malachite_q::Rational;

use crate::rings::linear::matrix::*;

impl Matrix<Rational> {
    pub fn gram_schmidt_col_orthogonalization(mut self) -> Self {
        //make the columns orthogonal by gram_schmidt gram_schmidt_col_orthogonalization
        for i in 0..self.cols() {
            for j in 0..i {
                //col_i = col_i - (projection of col_j onto col_i)
                let lambda = Self::dot(&self.get_col(i), &self.get_col(j))
                    / Self::dot(&self.get_col(j), &self.get_col(j));
                //col_i -= lambda col_j
                for r in 0..self.rows() {
                    *self.at_mut(r, i).unwrap() =
                        self.at(r, i).unwrap() - &lambda * self.at(r, j).unwrap();
                }
            }
        }
        for i in 0..self.cols() {
            for j in (i + 1)..self.cols() {
                debug_assert_eq!(
                    Self::dot(&self.get_col(i), &self.get_col(j)),
                    Rational::ZERO
                );
            }
        }
        self
    }

    pub fn gram_schmidt_row_orthogonalization(self) -> Self {
        self.transpose()
            .gram_schmidt_col_orthogonalization()
            .transpose()
    }
}

impl Matrix<Integer> {
    //applications TODO:
    // algorithm for the principal ideal problem
    // integer polynomial factorisation

    pub fn lll_col_algorithm(self, delta: &Rational) -> Self {
        //https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm

        // The algorithm is valid for 1/4 <= delta <= 1
        assert!(
            &Rational::from_integers(Integer::from(1), Integer::from(4)) <= delta
                && delta <= &Rational::ONE
        );
        //columns are a basis of a lattice
        debug_assert_eq!(self.rank(), self.cols());

        let mut a = self;
        let mut b = a
            .apply_map(|elem| Rational::from(elem))
            .gram_schmidt_col_orthogonalization();

        a.pprint();
        b.pprint();

        todo!()
    }

    pub fn lll_row_algorithm(self, delta: &Rational) -> Self {
        self.transpose().lll_col_algorithm(delta).transpose()
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use super::*;

    #[test]
    fn run() {
        let mat = Matrix::from_rows(vec![
            vec![Integer::from(1), Integer::from(-1), Integer::from(3)],
            vec![Integer::from(1), Integer::from(0), Integer::from(5)],
            vec![Integer::from(1), Integer::from(2), Integer::from(6)],
        ]);

        mat.lll_col_algorithm(&Rational::from_integers(Integer::from(3), Integer::from(4)));
    }
}
