#![allow(dead_code)]

use std::borrow::Borrow;

use super::matrix::*;
use super::ring::*;

//return a metamatrix whose rows are a basis for the joint row span of all the passed metamatricies
fn metamatrix_row_sum<R: PrincipalIdealDomain, MetaMatT: Borrow<Matrix<R::ElemT>>>(
    ring: R,
    cols: usize,
    metamats: Vec<MetaMatT>,
) -> Matrix<R::ElemT> {
    let mat_struct = MatrixStructure { ring };
    for metamat in &metamats {
        assert_eq!(metamat.borrow().cols(), cols);
    }
    let joined_metamat = Matrix::join_rows(cols, metamats);
    let (h, _u, _u_det, pivs) = mat_struct.row_hermite_algorithm(joined_metamat);
    h.submatrix((0..pivs.len()).collect(), (0..cols).collect()) //return the top non-zero and linearly independent rows from h
}

//return a metamatrix whose rows are a basis for the intersection of the row spans of the passed metamatricies
fn metamatrix_row_intersection<R: PrincipalIdealDomain, MetaMatT: Borrow<Matrix<R::ElemT>>>(
    ring: R,
    cols: usize,
    mut metamat1: MetaMatT,
    mut metamat2: MetaMatT,
) -> Matrix<R::ElemT> {
    let mat_struct = MatrixStructure { ring };
    assert_eq!(metamat1.borrow().cols(), cols);
    assert_eq!(metamat2.borrow().cols(), cols);
    //metamats should have linearly independent rows
    debug_assert_eq!(
        mat_struct.rank(metamat1.borrow().clone()),
        metamat1.borrow().rows()
    );
    debug_assert_eq!(
        mat_struct.rank(metamat2.borrow().clone()),
        metamat2.borrow().rows()
    );
    if metamat1.borrow().rows() > metamat2.borrow().rows() {
        //optimize for when we take linear combinations of rows from metamat1 rather than metamat2 later
        (metamat1, metamat2) = (metamat2, metamat1);
    }
    let joined_metamat = Matrix::join_rows(
        cols,
        vec![
            metamat1.borrow(),
            &mat_struct.neg(metamat2.borrow().clone()),
        ],
    );
    //the row kernel of joined_metamat tells us about which linear combinations of rows of metamat1 are equal to which linear combinations of rows of metamat2
    let row_ker = mat_struct.row_kernel(joined_metamat);
    //the rows in row_ker are in two halves
    //the first represents a linear combination of metamat1 rows
    //the second represents a linear combination of metamat2 rows
    //take without loss of generality all the linear combinations of metamat1 rows
    //the rows of the resulting matric are linearly independent
    //    because projection from (linear combinations of rows of metamat1 and metamat2) to (linear combinations of rows of metamat1) is injective
    //    because if a linear combination of rows of metamat1 and metamat2 is such that the metamat1 part is zero, then a linear combination of rows of metamat2 is zero
    //    but metamat2 is linearly independent, so the whole linear combination is zero
    Matrix::join_rows(
        cols,
        (0..row_ker.rank())
            .map(|i| {
                mat_struct
                    .mul_refs(
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

#[derive(Debug, Clone)]
pub struct LinearLattice<ElemT: Clone + PartialEq + Eq> {
    //matrix whose rows are a basis of the linear lattice
    //NOT necessarily in row hermite normal form
    metamatrix: Matrix<ElemT>,
    //each row represents a matrix of this shape
    rows: usize,
    cols: usize,
}

impl<ElemT: Clone + PartialEq + Eq> LinearLattice<ElemT> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }
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

fn mats_to_rows<ElemT: Clone, MatT: Borrow<Matrix<ElemT>>>(
    rows: usize,
    cols: usize,
    mats: Vec<MatT>,
) -> Matrix<ElemT> {
    for mat in &mats {
        assert_eq!(mat.borrow().rows(), rows);
        assert_eq!(mat.borrow().cols(), cols);
    }
    Matrix::construct(mats.len(), rows * cols, |r, c| {
        let (mr, mc) = idx_to_rc(rows, cols, c);
        mats[r].borrow().at(mr, mc).unwrap().clone()
    })
}

#[derive(Debug, Clone)]
pub struct LinearLatticeStructure<R: PrincipalIdealDomain> {
    ring: R,
}

impl<R: PrincipalIdealDomain> LinearLatticeStructure<R> {
    pub fn matrix_structure(&self) -> MatrixStructure<R> {
        MatrixStructure { ring: self.ring }
    }

    pub fn check_invariants(&self, lat: LinearLattice<R::ElemT>) -> Result<(), &'static str> {
        if lat.rows * lat.cols != lat.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if self.matrix_structure().rank(lat.metamatrix.clone()) != lat.metamatrix.rows() {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn from_span<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R::ElemT> {
        let spanning_meta_matrix = mats_to_rows(rows, cols, mats);
        let (h, _u, _u_det, pivs) = self
            .matrix_structure()
            .row_hermite_algorithm(spanning_meta_matrix);
        let metamatrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
        let lattice = LinearLattice {
            metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(lattice).is_ok());
        lattice
    }

    pub fn from_basis<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R::ElemT> {
        let metamatrix = mats_to_rows(rows, cols, mats);
        let lattice = LinearLattice {
            metamatrix: metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(lattice).is_ok());
        lattice
    }

    pub fn rank(&self, lat: &LinearLattice<R::ElemT>) -> usize {
        lat.metamatrix.rows()
    }

    fn basis_row(&self, lat: &LinearLattice<R::ElemT>, basis_num: usize) -> Matrix<R::ElemT> {
        lat.metamatrix
            .submatrix(vec![basis_num], (0..lat.metamatrix.cols()).collect())
    }

    pub fn basis_matrix(&self, lat: &LinearLattice<R::ElemT>, r: usize) -> Matrix<R::ElemT> {
        if self.rank(lat) <= r {
            panic!();
        }
        Matrix::construct(lat.rows, lat.cols, |mr, mc| {
            lat.metamatrix
                .at(r, rc_to_idx(lat.rows, lat.cols, mr, mc))
                .unwrap()
                .clone()
        })
    }

    pub fn basis_matrix_element(
        &self,
        lat: &LinearLattice<R::ElemT>,
        basis_num: usize,
        r: usize,
        c: usize,
    ) -> &R::ElemT {
        lat.metamatrix
            .at(basis_num, rc_to_idx(lat.rows, lat.cols, r, c))
            .unwrap()
    }

    fn contains_row<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        mat_as_row: MatT,
    ) -> bool {
        match self
            .matrix_structure()
            .row_solve(&lat.metamatrix, mat_as_row)
        {
            Some(_taps) => true,
            None => false,
        }
    }

    pub fn contains_point<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        mat: MatT,
    ) -> bool {
        self.contains_row(lat, mats_to_rows(lat.rows, lat.cols, vec![mat]))
    }

    //is lat a subset of self?
    pub fn contains_sublattice<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        sublat: LatT,
    ) -> bool {
        assert_eq!(lat.metamatrix.cols(), sublat.borrow().metamatrix.cols());
        for basis_num in 0..self.rank(sublat.borrow()) {
            if !self.contains_row(lat, self.basis_row(sublat.borrow(), basis_num)) {
                return false;
            }
        }
        true
    }

    pub fn sum<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R::ElemT> {
        for lat in &lats {
            assert_eq!(lat.borrow().rows, rows);
            assert_eq!(lat.borrow().cols, cols);
        }
        LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_sum(
                self.ring,
                rows * cols,
                lats.iter().map(|lat| &lat.borrow().metamatrix).collect(),
            ),
        }
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<R::ElemT> {
        self.sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R::ElemT> {
        if lats.len() == 0 {
            LinearLattice {
                rows,
                cols,
                metamatrix: self.matrix_structure().ident(rows * cols),
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = self.intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = self.intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    pub fn intersect_pair<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<R::ElemT> {
        assert_eq!(lat1.borrow().rows, rows);
        assert_eq!(lat1.borrow().cols, cols);
        assert_eq!(lat2.borrow().rows, rows);
        assert_eq!(lat2.borrow().cols, cols);
        let intersection_lattice = LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_intersection(
                self.ring,
                rows * cols,
                &lat1.borrow().metamatrix,
                &lat2.borrow().metamatrix,
            ),
        };
        debug_assert!(self.check_invariants(intersection_lattice).is_ok());
        intersection_lattice
    }
}

// impl<R: PrincipalIdealDomain> PartialEq for LinearLatticeStructure<R> {
//     fn eq(&self, other: &Self) -> bool {
//         self.contains_sublattice(other) && other.contains_sublattice(self)
//     }
// }
// impl<R: PrincipalIdealDomain> Eq for LinearLatticeStructure<R> {}

impl<R: PrincipalIdealDomain> LinearLatticeStructure<R> {
    pub fn pprint(&self, lat : LinearLattice<R::ElemT>) {
        println!("Start Linear Lattice");
        for r in 0..self.metamatrix.rows() {
            self.basis_matrix(lat, r).pprint();
        }
        println!("End Linear Lattice");
    }
}

#[derive(Debug, Clone)]
enum AffineLatticeElements<R: PrincipalIdealDomain> {
    Empty(),
    NonEmpty {
        offset: MatrixStructure<R>, //offset.rows == 1 and offset.cols == self.cols
        linlat: LinearLattice<R>,   //linlat.rows == self.rows and linlat.cols == self.cols
    },
}

#[derive(Debug, Clone)]
pub struct AffineLattice<R: PrincipalIdealDomain> {
    rows: usize,
    cols: usize,
    elems: AffineLatticeElements<R>,
}

impl<R: PrincipalIdealDomain> AffineLattice<R> {
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

    pub fn contains_point<MatT: Borrow<MatrixStructure<R>>>(&self, mat: MatT) -> bool {
        match &self.elems {
            AffineLatticeElements::Empty() => false,
            AffineLatticeElements::NonEmpty { offset, linlat } => linlat.contains_point(
                MatrixStructure::add_ref(offset.clone().neg(), mat.borrow()).unwrap(),
            ),
        }
    }

    //is other a subset of self?
    pub fn contains_sublattice<LatT: Borrow<AffineLattice<R>>>(&self, other: LatT) -> bool {
        match &other.borrow().elems {
            AffineLatticeElements::Empty() => true,
            AffineLatticeElements::NonEmpty {
                offset: other_offset,
                linlat: other_linlat,
            } => match &self.elems {
                AffineLatticeElements::Empty() => false,
                AffineLatticeElements::NonEmpty {
                    offset: _self_offset,
                    linlat: _self_linlat,
                } => {
                    for bn in 0..other_linlat.borrow().rank() {
                        if !self.contains_point(
                            MatrixStructure::add_ref(
                                other_linlat.borrow().basis_matrix(bn),
                                other_offset,
                            )
                            .unwrap(),
                        ) {
                            return false;
                        }
                    }
                    true
                }
            },
        }
    }

    pub fn sum<LatT: Borrow<AffineLattice<R>>>(rows: usize, cols: usize, lats: Vec<LatT>) -> Self {
        let mut sum_offset = MatrixStructure::zero(rows, cols);
        let mut sum_linlats = vec![];
        for lat in &lats {
            assert_eq!(lat.borrow().rows(), rows);
            assert_eq!(lat.borrow().cols(), cols);
            match &lat.borrow().elems {
                AffineLatticeElements::Empty() => {
                    return AffineLattice::empty(rows, cols);
                }
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    sum_offset.add_mut(offset).unwrap();
                    sum_linlats.push(linlat);
                }
            }
        }
        AffineLattice {
            rows: rows,
            cols: cols,
            elems: AffineLatticeElements::NonEmpty {
                offset: sum_offset,
                linlat: LinearLattice::sum(rows, cols, sum_linlats),
            },
        }
    }

    pub fn sum_pair<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        Self::sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> Self {
        if lats.len() == 0 {
            Self {
                rows,
                cols,
                elems: AffineLatticeElements::NonEmpty {
                    offset: MatrixStructure::zero(rows, cols),
                    linlat: LinearLattice {
                        rows,
                        cols,
                        metamatrix: MatrixStructure::ident(rows * cols),
                    },
                },
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = Self::intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = Self::intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    pub fn intersect_pair<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        assert_eq!(lat1.borrow().rows(), rows);
        assert_eq!(lat1.borrow().cols(), cols);
        assert_eq!(lat2.borrow().rows(), rows);
        assert_eq!(lat2.borrow().cols(), cols);
        match &lat1.borrow().elems {
            AffineLatticeElements::Empty() => AffineLattice::empty(rows, cols),
            AffineLatticeElements::NonEmpty {
                offset: offset1,
                linlat: linlat1,
            } => match &lat2.borrow().elems {
                AffineLatticeElements::Empty() => AffineLattice::empty(rows, cols),
                AffineLatticeElements::NonEmpty {
                    offset: offset2,
                    linlat: linlat2,
                } => {
                    //model an affine lattice as the intersection of a linear lattice (henceforth called the hyperlattice) living in one higher dimension with the plane (1, *, *, ..., *)
                    //take the intersection of the linear lattices and row reduce and get something like
                    // a * * * *
                    // 0 0 b * *
                    // 0 0 0 c *
                    //if a=1 then the rest of the top row is the offset and the bottom right submatrix is the basis of the linear lattice

                    fn offset_linlat_to_metamat<R: PrincipalIdealDomain>(
                        rows: usize,
                        cols: usize,
                        offset: &MatrixStructure<R>,
                        linlat: &LinearLattice<R>,
                    ) -> MatrixStructure<R> {
                        let mut metamat = MatrixStructure::zero(1 + linlat.rank(), 1 + rows * cols);
                        *metamat.at_mut(0, 0).unwrap() = R::one();
                        for idx in 0..rows * cols {
                            let (r, c) = idx_to_rc(rows, cols, idx);
                            println!("rows={} cols={} r={} c={} idx={}", rows, cols, r, c, idx);
                            *metamat.at_mut(0, 1 + idx).unwrap() = offset.at(r, c).unwrap().clone();
                        }
                        for bn in 0..linlat.rank() {
                            for idx in 0..rows * cols {
                                let (r, c) = idx_to_rc(rows, cols, idx);
                                *metamat.at_mut(0 + 1 + bn, 1 + idx).unwrap() =
                                    linlat.basis_matrix_element(bn, r, c).clone();
                            }
                        }
                        metamat
                    }
                    let metamat1: MatrixStructure<R> =
                        offset_linlat_to_metamat(rows, cols, offset1, linlat1);
                    let metamat2: MatrixStructure<R> =
                        offset_linlat_to_metamat(rows, cols, offset2, linlat2);
                    let int_metamat =
                        metamatrix_row_intersection(1 + rows * cols, metamat1, metamat2);

                    if int_metamat.rows() == 0 {
                        //the hyperlattice is just the origin, so the coresponding affine lattice - the intersection with the plane (1, *, ..., *) - is empty.
                        AffineLattice::empty(rows, cols)
                    } else {
                        let (int_metamat_h, _u, _u_det, pivs) = int_metamat.row_hermite_algorithm();
                        int_metamat_h.pprint();
                        if R::is_unit(int_metamat_h.at(0, 0).unwrap().clone()) {
                            debug_assert_eq!(int_metamat_h.at(0, 0).unwrap(), &R::one());
                        }
                        if int_metamat_h.at(0, 0).unwrap() == &R::one() {
                            let mut int_offset = MatrixStructure::zero(rows, cols);
                            for idx in 0..rows * cols {
                                let (r, c) = idx_to_rc(rows, cols, idx);
                                *int_offset.at_mut(r, c).unwrap() =
                                    int_metamat_h.at(0, 1 + idx).unwrap().clone();
                            }
                            let int_basis_mats: Vec<MatrixStructure<R>> = (0..pivs.len() - 1)
                                .map(|bn| {
                                    debug_assert_eq!(
                                        int_metamat_h.at(1 + bn, 0).unwrap(),
                                        &R::zero()
                                    );
                                    let mut basis_mat = MatrixStructure::zero(rows, cols);
                                    for idx in 0..rows * cols {
                                        let (r, c) = idx_to_rc(rows, cols, idx);
                                        *basis_mat.at_mut(r, c).unwrap() =
                                            int_metamat_h.at(1 + bn, 1 + idx).unwrap().clone();
                                    }
                                    basis_mat
                                })
                                .collect();
                            int_offset.pprint();
                            for basis_mat in &int_basis_mats {
                                basis_mat.pprint();
                            }
                            Self::from_offset_and_linear_lattice(
                                rows,
                                cols,
                                int_offset,
                                LinearLattice::from_basis(rows, cols, int_basis_mats),
                            )
                        } else {
                            //the hyperlattice does not intersect the plane (1, *, ..., *) because int_metamat_h(0, 0) is not a unit
                            AffineLattice::empty(rows, cols)
                        }
                    }
                }
            },
        }
    }
}

impl<R: PrincipalIdealDomain> PartialEq for AffineLattice<R> {
    fn eq(&self, other: &Self) -> bool {
        self.contains_sublattice(other) && other.contains_sublattice(self)
    }
}
impl<R: PrincipalIdealDomain> Eq for AffineLattice<R> {}

impl<R: PrincipalIdealDomain> AffineLattice<R> {
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

    use super::super::nzq::*;
    use super::*;

    #[test]
    fn linear_lattice_invariant() {
        let lattice = LinearLattice::<IntegerRing> {
            metamatrix: MatrixStructure::from_rows(vec![
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

        let lattice = LinearLattice::<IntegerRing> {
            metamatrix: MatrixStructure::from_rows(vec![
                vec![Integer::from(0), Integer::from(3), Integer::from(0)],
                vec![Integer::from(2), Integer::from(0), Integer::from(1)],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(lattice.check_invariants().is_err());

        let lattice = LinearLattice::<IntegerRing> {
            metamatrix: MatrixStructure::from_rows(vec![
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
        let lattice = LinearLattice::<IntegerRing>::from_span(
            2,
            2,
            vec![
                &MatrixStructure::from_rows(vec![
                    vec![Integer::from(0), Integer::from(3)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
                &MatrixStructure::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(1)],
                ]),
            ],
        );

        assert_eq!(
            true,
            lattice.contains_point(&MatrixStructure::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&MatrixStructure::from_rows(vec![
                vec![Integer::from(2), Integer::from(4)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&MatrixStructure::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(1), Integer::from(1)],
            ]))
        );

        assert_ne!(
            LinearLattice::from_span(
                2,
                3,
                vec![
                    &MatrixStructure::<IntegerRing>::from_rows(vec![
                        vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &MatrixStructure::<IntegerRing>::from_rows(vec![
                        vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            ),
            LinearLattice::from_span(
                2,
                3,
                vec![
                    &MatrixStructure::from_rows(vec![
                        vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &MatrixStructure::from_rows(vec![
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
            let a = MatrixStructure::<IntegerRing>::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let b = MatrixStructure::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            let c = MatrixStructure::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let d = MatrixStructure::from_rows(vec![
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
            let a = MatrixStructure::<IntegerRing>::from_rows(vec![
                vec![Integer::from(3), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(5), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(7)],
            ]);

            let b = MatrixStructure::from_rows(vec![
                vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let c = MatrixStructure::from_rows(vec![
                vec![Integer::from(6), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(20), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let d = MatrixStructure::from_rows(vec![
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
            //triple intersection
            let a = MatrixStructure::<IntegerRing>::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
            ]);

            let b = MatrixStructure::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
            ]);

            let c = MatrixStructure::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
            ]);

            let int =
                LinearLattice::intersect(4, 1, vec![a.col_span(), b.col_span(), c.col_span()]);

            assert_eq!(
                int,
                MatrixStructure::from_rows(vec![
                    vec![
                        Integer::from(1),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                ])
                .col_span()
            );
        }

        {
            //complex example
            let a = MatrixStructure::<IntegerRing>::from_rows(vec![
                vec![Integer::from(3), Integer::from(9), Integer::from(27)],
                vec![Integer::from(-4), Integer::from(6), Integer::from(-100)],
                vec![Integer::from(2), Integer::from(8), Integer::from(7)],
            ]);

            let b = MatrixStructure::from_rows(vec![
                vec![Integer::from(12), Integer::from(-1), Integer::from(18)],
                vec![Integer::from(-5), Integer::from(12), Integer::from(-24)],
                vec![Integer::from(1), Integer::from(2), Integer::from(14)],
            ]);

            let c = MatrixStructure::from_rows(vec![
                vec![
                    Integer::from(21),
                    Integer::from(3852),
                    Integer::from(3315300),
                ],
                vec![
                    Integer::from(-252),
                    Integer::from(-46214),
                    Integer::from(-39775000),
                ],
                vec![
                    Integer::from(-42),
                    Integer::from(-4454),
                    Integer::from(-3833450),
                ],
            ]);

            let d = MatrixStructure::from_rows(vec![
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
        let afflat = AffineLattice::<IntegerRing> {
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
                        MatrixStructure::<IntegerRing>::from_rows(vec![
                            vec![Integer::from(1), Integer::from(0)],
                            vec![Integer::from(0), Integer::from(1)],
                        ]),
                        MatrixStructure::from_rows(vec![
                            vec![Integer::from(0), Integer::from(1)],
                            vec![Integer::from(1), Integer::from(0)],
                        ]),
                    ],
                ),
                offset: MatrixStructure::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
            },
        };
        assert!(afflat.check_invariants().is_ok());
    }

    fn affine_lattice_sum_and_intersection() {
        let a1 = MatrixStructure::<IntegerRing>::from_rows(vec![
            vec![Integer::from(3), Integer::from(1), Integer::from(0)],
            vec![Integer::from(3), Integer::from(1), Integer::from(0)],
            vec![Integer::from(3), Integer::from(1), Integer::from(1)],
        ]);

        let y1 = MatrixStructure::from_rows(vec![
            vec![Integer::from(1)],
            vec![Integer::from(1)],
            vec![Integer::from(1)],
        ]);

        let a2 = MatrixStructure::from_rows(vec![
            vec![Integer::from(3), Integer::from(5), Integer::from(0)],
            vec![Integer::from(3), Integer::from(5), Integer::from(0)],
            vec![Integer::from(0), Integer::from(0), Integer::from(1)],
        ]);

        let y2 = MatrixStructure::from_rows(vec![
            vec![Integer::from(1)],
            vec![Integer::from(1)],
            vec![Integer::from(1)],
        ]);

        let alat1 = a1.col_solution_lattice(y1);
        let alat2 = a2.col_solution_lattice(y2);

        alat1.pprint();
        println!();
        alat2.pprint();

        let alat3 = AffineLattice::sum(3, 1, vec![alat1, alat2]);
        println!();
        alat3.pprint();

        let expected_alat3 = AffineLattice::from_offset_and_linear_lattice(
            3,
            1,
            MatrixStructure::from_rows(vec![
                vec![Integer::from(2)],
                vec![Integer::from(0)],
                vec![Integer::from(1)],
            ]),
            LinearLattice::from_span(
                3,
                1,
                vec![
                    MatrixStructure::from_rows(vec![
                        vec![Integer::from(1)],
                        vec![Integer::from(-3)],
                        vec![Integer::from(0)],
                    ]),
                    MatrixStructure::from_rows(vec![
                        vec![Integer::from(5)],
                        vec![Integer::from(-3)],
                        vec![Integer::from(0)],
                    ]),
                ],
            ),
        );

        assert_eq!(alat3, expected_alat3);
    }
}
