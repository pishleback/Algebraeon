#![allow(dead_code)]

use std::borrow::Borrow;

use super::matrix::*;
use super::ring::*;

//return a metamatrix whose rows are a basis for the joint row span of all the passed metamatricies
fn metamatrix_row_sum<R: PrincipalIdealDomain + std::fmt::Display, MetaMatT: Borrow<Matrix<R>>>(
    cols: usize,
    metamats: Vec<MetaMatT>,
) -> Matrix<R> {
    for metamat in &metamats {
        assert_eq!(metamat.borrow().cols(), cols);
    }
    let joined_metamat = join_rows(cols, metamats);
    let (h, _u, pivs) = joined_metamat.row_hermite_algorithm();
    h.submatrix((0..pivs.len()).collect(), (0..cols).collect()) //return the top non-zero and linearly independent rows from h
}

//return a metamatrix whose rows are a basis for the intersection of the row spans of the passed metamatricies
fn metamatrix_row_intersection<
    R: PrincipalIdealDomain + std::fmt::Display,
    MetaMatT: Borrow<Matrix<R>>,
>(
    cols: usize,
    mut metamat1: MetaMatT,
    mut metamat2: MetaMatT,
) -> Matrix<R> {
    assert_eq!(metamat1.borrow().cols(), cols);
    assert_eq!(metamat2.borrow().cols(), cols);
    //metamats should have linearly independent rows
    debug_assert_eq!(metamat1.borrow().clone().rank(), metamat1.borrow().rows());
    debug_assert_eq!(metamat2.borrow().clone().rank(), metamat2.borrow().rows());
    if metamat1.borrow().rows() > metamat2.borrow().rows() {
        //optimize for when we take linear combinations of rows from metamat1 rather than metamat2 later
        (metamat1, metamat2) = (metamat2, metamat1);
    }
    let joined_metamat = join_rows(
        cols,
        vec![metamat1.borrow(), &metamat2.borrow().clone().neg()],
    );
    //the row kernel of joined_metamat tells us about which linear combinations of rows of metamat1 are equal to which linear combinations of rows of metamat2
    let row_ker = joined_metamat.row_kernel();
    //the rows in row_ker are in two halves
    //the first represents a linear combination of metamat1 rows
    //the second represents a linear combination of metamat2 rows
    //take without loss of generality all the linear combinations of metamat1 rows
    //the rows of the resulting matric are linearly independent
    //    because projection from (linear combinations of rows of metamat1 and metamat2) to (linear combinations of rows of metamat1) is injective
    //    because if a linear combination of rows of metamat1 and metamat2 is such that the metamat1 part is zero, then a linear combination of rows of metamat2 is zero
    //    but metamat2 is linearly independent, so the whole linear combination is zero
    join_rows(
        cols,
        (0..row_ker.rank())
            .map(|i| {
                Matrix::mul_refs(
                    &row_ker
                        .basis_matrix(i)
                        .submatrix(vec![0], (0..metamat1.borrow().rows()).collect()),
                    metamat1.borrow(),
                )
                .unwrap()
            })
            .collect(),
    )
}

#[derive(Debug)]
pub struct LinearLattice<R: PrincipalIdealDomain> {
    //matrix whose rows are a basis of the linear lattice
    //NOT necessarily in row hermite normal form
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

fn mats_to_rows<R: PrincipalIdealDomain + std::fmt::Display, MatT: Borrow<Matrix<R>>>(
    rows: usize,
    cols: usize,
    mats: Vec<MatT>,
) -> Matrix<R> {
    for mat in &mats {
        assert_eq!(mat.borrow().rows(), rows);
        assert_eq!(mat.borrow().cols(), cols);
    }
    let mut mats_as_rows: Matrix<R> = Matrix::zero(mats.len(), rows * cols);
    for (r, mat) in mats.into_iter().enumerate() {
        for mr in 0..rows {
            for mc in 0..cols {
                *mats_as_rows
                    .at_mut(r, rc_to_idx(rows, cols, mr, mc))
                    .unwrap() = mat.borrow().at(mr, mc).unwrap().clone();
            }
        }
    }
    mats_as_rows
}

impl<R: PrincipalIdealDomain + std::fmt::Display> LinearLattice<R> {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if self.rows * self.cols != self.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if self.metamatrix.clone().rank() != self.metamatrix.rows() {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn from_span<MatT: Borrow<Matrix<R>>>(rows: usize, cols: usize, mats: Vec<MatT>) -> Self {
        let spanning_meta_matrix = mats_to_rows(rows, cols, mats);
        let (h, _u, pivs) = spanning_meta_matrix.row_hermite_algorithm();
        let metamatrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
        let lattice = Self {
            metamatrix,
            rows,
            cols,
        };
        debug_assert!(lattice.check_invariants().is_ok());
        lattice
    }

    pub fn from_basis<MatT: Borrow<Matrix<R>>>(rows: usize, cols: usize, mats: Vec<MatT>) -> Self {
        let metamatrix = mats_to_rows(rows, cols, mats);
        let lattice = Self {
            metamatrix: metamatrix,
            rows,
            cols,
        };
        debug_assert!(lattice.check_invariants().is_ok());
        lattice
    }

    pub fn rank(&self) -> usize {
        self.metamatrix.rows()
    }

    fn basis_row(&self, basis_num: usize) -> Matrix<R> {
        self.metamatrix
            .submatrix(vec![basis_num], (0..self.metamatrix.cols()).collect())
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

    fn contains_row<MatT: Borrow<Matrix<R>>>(&self, mat_as_row: MatT) -> bool {
        match self.metamatrix.row_solve(mat_as_row) {
            Some(_taps) => true,
            None => false,
        }
    }

    pub fn contains_point<MatT: Borrow<Matrix<R>>>(&self, mat: MatT) -> bool {
        self.contains_row(mats_to_rows(self.rows, self.cols, vec![mat]))
    }

    //is lat a subset of self?
    pub fn contains_lattice<LatT: Borrow<LinearLattice<R>>>(&self, lat: LatT) -> bool {
        assert_eq!(self.metamatrix.cols(), lat.borrow().metamatrix.cols());
        for basis_num in 0..lat.borrow().rank() {
            if !self.contains_row(lat.borrow().basis_row(basis_num)) {
                return false;
            }
        }
        true
    }

    pub fn sum<LatT: Borrow<LinearLattice<R>>>(rows: usize, cols: usize, lats: Vec<LatT>) -> Self {
        for lat in &lats {
            assert_eq!(lat.borrow().rows, rows);
            assert_eq!(lat.borrow().cols, cols);
        }
        Self {
            rows,
            cols,
            metamatrix: metamatrix_row_sum(
                rows * cols,
                lats.iter().map(|lat| &lat.borrow().metamatrix).collect(),
            ),
        }
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        Self::sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect_pair<LatT: Borrow<LinearLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        assert_eq!(lat1.borrow().rows, rows);
        assert_eq!(lat1.borrow().cols, cols);
        assert_eq!(lat2.borrow().rows, rows);
        assert_eq!(lat2.borrow().cols, cols);
        let intersection_lattice = Self {
            rows,
            cols,
            metamatrix: metamatrix_row_intersection(
                rows * cols,
                &lat1.borrow().metamatrix,
                &lat2.borrow().metamatrix,
            ),
        };
        debug_assert!(intersection_lattice.check_invariants().is_ok());
        intersection_lattice
    }
}

impl<R: PrincipalIdealDomain + std::fmt::Display> PartialEq for LinearLattice<R> {
    fn eq(&self, other: &Self) -> bool {
        self.contains_lattice(other) && other.contains_lattice(self)
    }
}
impl<R: PrincipalIdealDomain + std::fmt::Display> Eq for LinearLattice<R> {}

impl<R: PrincipalIdealDomain + std::fmt::Display> LinearLattice<R> {
    pub fn pprint(&self) {
        println!("Start Linear Lattice");
        for r in 0..self.metamatrix.rows() {
            self.basis_matrix(r).pprint();
        }
        println!("End Linear Lattice");
    }
}

#[derive(Debug)]
enum AffineLatticeElements<R: PrincipalIdealDomain> {
    Empty(),
    NonEmpty {
        offset: Matrix<R>,        //offset.rows == 1 and offset.cols == self.cols
        linlat: LinearLattice<R>, //linlat.rows == self.rows and linlat.cols == self.cols
    },
}

#[derive(Debug)]
pub struct AffineLattice<R: PrincipalIdealDomain> {
    rows: usize,
    cols: usize,
    elems: AffineLatticeElements<R>,
}

impl<R: PrincipalIdealDomain + std::fmt::Display> AffineLattice<R> {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match &self.elems {
            AffineLatticeElements::Empty() => {}
            AffineLatticeElements::NonEmpty { offset, linlat } => match linlat.check_invariants() {
                Ok(()) => {
                    if offset.rows() != self.rows {
                        return Err("offset rows doesnt match self rows");
                    }
                    if offset.cols() != self.cols {
                        return Err("offset columns doesnt match self columns");
                    }
                    if linlat.rows() != self.rows {
                        return Err("linlat rows doesnt match self rows");
                    }
                    if linlat.cols() != self.cols {
                        return Err("linlat columns doesnt match self columns");
                    }
                }
                Err(msg) => {
                    return Err(msg);
                }
            },
        }
        Ok(())
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn rank(&self) -> Option<usize> {
        match &self.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty {
                offset: _offset,

                linlat,
            } => Some(linlat.rank()),
        }
    }

    pub fn empty(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            elems: AffineLatticeElements::Empty(),
        }
    }

    pub fn from_offset_and_linear_lattice(
        rows: usize,
        cols: usize,
        offset: Matrix<R>,
        linlat: LinearLattice<R>,
    ) -> Self {
        assert_eq!(offset.rows(), rows);
        assert_eq!(offset.cols(), cols);
        assert_eq!(linlat.rows(), rows);
        assert_eq!(linlat.cols(), cols);
        let afflat = Self {
            rows,
            cols,
            elems: AffineLatticeElements::NonEmpty { offset, linlat },
        };
        debug_assert!(afflat.check_invariants().is_ok());
        afflat
    }
}

impl<R: PrincipalIdealDomain + std::fmt::Display> AffineLattice<R> {
    pub fn pprint(&self) {
        println!("Start Affine Lattice");
        match &self.elems {
            AffineLatticeElements::Empty() => println!("Empty"),
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                println!("Offset");
                offset.pprint();
                linlat.pprint();
            }
        }
        println!("End Affine Lattice");
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::*;

    #[test]
    fn linear_lattice_invariant() {
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

    #[test]
    fn containment() {
        let lattice = LinearLattice::from_span(
            2,
            2,
            vec![
                &Matrix::from_rows(vec![
                    vec![Integer::from(0), Integer::from(3)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
                &Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(1)],
                ]),
            ],
        );

        assert_eq!(
            true,
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(4)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(1), Integer::from(1)],
            ]))
        );

        assert_ne!(
            LinearLattice::from_span(
                2,
                3,
                vec![
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            ),
            LinearLattice::from_span(
                2,
                3,
                vec![
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            )
        );
    }

    #[test]
    fn linear_lattice_sum_and_intersection() {
        {
            //standard basis sum and intersection
            let a = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            println!();
            println!("a");
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int = LinearLattice::intersect_pair(3, 1, a.col_span(), b.col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.col_span(), b.col_span());
            sum.pprint();

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }

        {
            //sum and intersection as gcd and lcm
            let a = Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(5), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(7)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![Integer::from(6), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(20), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(7)],
            ]);

            println!();
            println!("a");
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int = LinearLattice::intersect_pair(3, 1, a.col_span(), b.col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.col_span(), b.col_span());
            sum.pprint();

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }

        {
            //complex example
            let a = Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(9), Integer::from(27)],
                vec![Integer::from(-4), Integer::from(6), Integer::from(-100)],
                vec![Integer::from(2), Integer::from(8), Integer::from(7)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(12), Integer::from(-1), Integer::from(18)],
                vec![Integer::from(-5), Integer::from(12), Integer::from(-24)],
                vec![Integer::from(1), Integer::from(2), Integer::from(14)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![Integer::from(21), Integer::from(3852), Integer::from(3315300)],
                vec![Integer::from(-252), Integer::from(-46214), Integer::from(-39775000)],
                vec![Integer::from(-42), Integer::from(-4454), Integer::from(-3833450)],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(-12), Integer::from(1), Integer::from(0)],
                vec![Integer::from(-2), Integer::from(325), Integer::from(1)],
            ]);

            println!();
            println!("a");
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int = LinearLattice::intersect_pair(3, 1, a.col_span(), b.col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.col_span(), b.col_span());
            sum.pprint();

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }
    }

    #[test]
    fn affine_lattice_invariants() {
        let afflat = AffineLattice::<Integer> {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::Empty(),
        };
        assert!(afflat.check_invariants().is_ok());

        let afflat = AffineLattice {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::NonEmpty {
                linlat: LinearLattice::from_basis(
                    2,
                    2,
                    vec![
                        Matrix::from_rows(vec![
                            vec![Integer::from(1), Integer::from(0)],
                            vec![Integer::from(0), Integer::from(1)],
                        ]),
                        Matrix::from_rows(vec![
                            vec![Integer::from(0), Integer::from(1)],
                            vec![Integer::from(1), Integer::from(0)],
                        ]),
                    ],
                ),
                offset: Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
            },
        };
        assert!(afflat.check_invariants().is_ok());
    }
}
