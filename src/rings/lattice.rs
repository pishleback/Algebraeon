#![allow(dead_code)]

use super::matrix::*;
use super::ring::*;

pub struct LinearLattice<R: PrincipalIdealDomain> {
    //matrix whose rows are a basis of the linear lattice
    metamatrix: Matrix<R>,
    //each row represents a matrix of this shape
    rows: usize,
    cols: usize,
}

//from matrix coords to meta row index
fn rc_to_idx(rows: usize, cols: usize, r: usize, c: usize) -> usize {
    if rows <= r || cols <= c {
        panic!();
    }
    c + r * cols
}

//from meta row index to matrix coords
fn idx_to_rc(rows: usize, cols: usize, idx: usize) -> (usize, usize) {
    if rows * cols <= idx {
        panic!();
    }
    (idx / cols, idx % cols)
}

// fn metamatrix_join_matrix<R: PrincipalIdealDomain>(
//     metamatrix: &Matrix<R>,
//     mat: &Matrix<R>,
// ) -> Matrix<R> {
//     assert_eq!(metamatrix.cols(), mat.rows() * mat.cols());
//     let mut extended_metamatrix = Matrix::zero(metamatrix.rows() + 1, metamatrix.cols());
//     for r in 0..metamatrix.rows() {
//         for c in 0..metamatrix.cols() {
//             *extended_metamatrix.at_mut(r, c).unwrap() = metamatrix.at(r, c).unwrap().clone();
//         }
//     }
//     for c in 0..metamatrix.cols() {
//         let (mr, mc) = idx_to_rc(mat.rows(), mat.cols(), c);
//         *extended_metamatrix.at_mut(metamatrix.rows(), c).unwrap() =
//             mat.at(mr, mc).unwrap().clone();
//     }
//     extended_metamatrix
// }

impl<R: PrincipalIdealDomain + std::fmt::Display> LinearLattice<R> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if self.rows * self.cols != self.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if self.metamatrix.clone().rank() != self.metamatrix.rows() {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn from_span(mats: Vec<Matrix<R>>) -> Self {
        if mats.len() == 0 {
            panic!();
        } else {
            let rows = mats[0].rows();
            let cols = mats[0].cols();
            for mat in &mats {
                if mat.rows() != rows {
                    panic!();
                }
                if mat.cols() != cols {
                    panic!();
                }
            }
            let mut spanning_meta_matrix: Matrix<R> = Matrix::zero(mats.len(), rows * cols);
            for (r, mat) in mats.into_iter().enumerate() {
                for mr in 0..rows {
                    for mc in 0..cols {
                        *spanning_meta_matrix
                            .at_mut(r, rc_to_idx(rows, cols, mr, mc))
                            .unwrap() = mat.at(mr, mc).unwrap().clone();
                    }
                }
            }
            let (h, _u, pivs) = spanning_meta_matrix.row_hermite_algorithm();
            let meta_matrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
            let lattice = Self {
                metamatrix: meta_matrix,
                rows,
                cols,
            };

            debug_assert!(lattice.check_invariants().is_ok());
            lattice
        }
    }

    pub fn rank(&self) -> usize {
        self.metamatrix.rows()
    }

    pub fn basis_matrix(&self, r: usize) -> Matrix<R> {
        if self.rank() <= r {
            panic!();
        }
        let mut mat = Matrix::zero(self.rows, self.cols);
        for mr in 0..self.rows {
            for mc in 0..self.cols {
                *mat.at_mut(mr, mc).unwrap() = self
                    .metamatrix
                    .at(r, rc_to_idx(self.rows, self.cols, mr, mc))
                    .unwrap()
                    .clone();
            }
        }
        mat
    }

    pub fn basis_matrix_element(&self, basis_num: usize, r: usize, c: usize) -> &R {
        self.metamatrix
            .at(basis_num, rc_to_idx(self.rows, self.cols, r, c))
            .unwrap()
    }
}

impl<R: PrincipalIdealDomain + std::fmt::Display> LinearLattice<R> {
    pub fn contains(&self, mat: &Matrix<R>) -> bool {
        todo!();
        // assert_eq!(self.rows, mat.rows());
        // assert_eq!(self.cols, mat.cols());
        // let ext_metamat = metamatrix_join_matrix(&self.metamatrix, mat);
        // println!();
        // self.metamatrix.pprint();
        // println!("{}", self.metamatrix.clone().rank());
        // ext_metamat.pprint();
        // println!("{}", ext_metamat.clone().rank());
        // let (h, u, pivs) = ext_metamat.clone().row_reduced_hermite_algorithm();
        // h.pprint();
        // println!("{}", h.clone().rank());

        // return ext_metamat.rank() == self.rank();
    }
}

impl<R: PrincipalIdealDomain> PartialEq for LinearLattice<R> {
    fn eq(&self, other: &Self) -> bool {
        todo!();
    }
}
impl<R: PrincipalIdealDomain> Eq for LinearLattice<R> {}

impl<R: PrincipalIdealDomain + std::fmt::Display> LinearLattice<R> {
    pub fn pprint(&self) {
        println!("Start Linear Span");
        for r in 0..self.metamatrix.rows() {
            self.basis_matrix(r).pprint();
        }
        println!("End Linear Span");
    }
}

// pub fn linlat_sum<R: PrincipalIdealDomain>(
//     lattices: Vec<LinearLattice<R>>,
// ) -> LinearLattice<R> {
//     todo!();
// }

// pub fn linlat_intersection<R: PrincipalIdealDomain>(
//     lattices: Vec<LinearLattice<R>>,
// ) -> LinearLattice<R> {
//     todo!();
// }

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use super::*;

    #[test]
    fn invariant() {
        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(2),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
            ]),
            rows: 2,
            cols: 2,
        };
        lattice.check_invariants().unwrap();

        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(3), Integer::from(0)],
                vec![Integer::from(2), Integer::from(0), Integer::from(1)],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(lattice.check_invariants().is_err());

        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![
                    Integer::from(6),
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(2),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(lattice.check_invariants().is_err());
    }

    #[test]
    fn index_conversions() {
        let rows = 5;
        let cols = 7;

        for idx in 0..rows * cols {
            let (r, c) = idx_to_rc(rows, cols, idx);
            assert_eq!(idx, rc_to_idx(rows, cols, r, c));
        }
        for r in 0..rows {
            for c in 0..cols {
                let idx = rc_to_idx(rows, cols, r, c);
                assert_eq!((r, c), idx_to_rc(rows, cols, idx));
            }
        }
    }

    // #[test]
    fn containment() {
        let lattice = LinearLattice::from_span(vec![
            Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(3)],
                vec![Integer::from(0), Integer::from(0)],
            ]),
            Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1)],
            ]),
        ]);

        assert_eq!(
            true,
            lattice.contains(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(4)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(1), Integer::from(1)],
            ]))
        );
    }
}
