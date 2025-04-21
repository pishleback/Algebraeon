use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix,
};
use crate::{linear::matrix::*, structure::*};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmodule<Ring: ReducedHermiteAlgorithmSignature> {
    module: FinitelyFreeModuleStructure<Ring>,
    // a matrix in reduced hermite normal form with all non-zero rows whose rows form a basis for the submodule
    row_basis: Matrix<Ring::Set>,
    // the columns of the pivots of row_basis
    pivots: Vec<usize>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmodule<Ring> {
    pub fn ring(&self) -> &Ring {
        self.module().ring()
    }

    pub fn module_rank(&self) -> usize {
        self.row_basis.cols()
    }

    pub fn submodule_rank(&self) -> usize {
        self.row_basis.rows()
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }

    pub fn submodule(&self) -> FinitelyFreeModuleStructure<Ring> {
        FinitelyFreeModuleStructure::new(self.ring().clone(), self.submodule_rank())
    }

    pub fn into_row_basis_matrix(self) -> Matrix<Ring::Set> {
        self.row_basis
    }

    pub fn into_col_basis_matrix(self) -> Matrix<Ring::Set> {
        self.row_basis.transpose()
    }

    pub(crate) fn row_basis_matrix(&self) -> &Matrix<Ring::Set> {
        &self.row_basis
    }

    pub fn basis(&self) -> Vec<Vec<Ring::Set>> {
        (0..self.row_basis.rows())
            .map(|r| {
                (0..self.module().rank())
                    .map(|c| self.row_basis.at(r, c).unwrap().clone())
                    .collect()
            })
            .collect()
    }

    pub fn matrix_row_span_and_basis(
        ring: Ring,
        matrix: Matrix<Ring::Set>,
    ) -> (Self, Matrix<Ring::Set>) {
        let rows = matrix.rows();
        let cols = matrix.cols();
        let (h, u, _det, pivots) =
            MatrixStructure::new(ring.clone()).row_reduced_hermite_algorithm(matrix);
        let row_basis = h.submatrix((0..pivots.len()).collect(), (0..cols).collect());
        let u = u.submatrix((0..pivots.len()).collect(), (0..rows).collect());
        (
            FinitelyFreeSubmodule {
                module: FinitelyFreeModuleStructure::new(ring, row_basis.cols()),
                row_basis,
                pivots,
            },
            u,
        )
    }

    pub fn matrix_col_span_and_basis(
        ring: Ring,
        matrix: Matrix<Ring::Set>,
    ) -> (Self, Matrix<Ring::Set>) {
        let (s, u) = Self::matrix_row_span_and_basis(ring, matrix.transpose());
        (s, u.transpose())
    }

    pub fn matrix_row_span(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        Self::matrix_row_span_and_basis(ring, matrix).0
    }

    pub fn matrix_col_span(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        Self::matrix_col_span_and_basis(ring, matrix).0
    }

    pub fn from_span(
        module: FinitelyFreeModuleStructure<Ring>,
        span: Vec<&Vec<Ring::Set>>,
    ) -> Self {
        for v in &span {
            debug_assert_eq!(v.len(), module.rank());
        }
        Self::matrix_row_span(
            module.ring().clone(),
            Matrix::construct(span.len(), module.rank(), |r, c| span[r][c].clone()),
        )
    }

    pub fn matrix_row_kernel(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        let rows = matrix.rows();
        let (_h, u, _u_det, pivs) =
            MatrixStructure::new(ring.clone()).row_hermite_algorithm(matrix);
        debug_assert_eq!(rows, u.rows());
        debug_assert_eq!(rows, u.cols());
        let ker = u.submatrix((pivs.len()..rows).collect(), (0..rows).collect());
        Self::matrix_row_span(ring, ker)
    }

    pub fn matrix_col_kernel(ring: Ring, matrix: Matrix<Ring::Set>) -> Self {
        Self::matrix_row_kernel(ring, matrix.transpose())
    }

    pub fn zero_submodule(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        let cols = module.rank();
        Self {
            module: module,
            row_basis: Matrix::construct(0, cols, |_, _| unreachable!()),
            pivots: vec![],
        }
    }

    pub fn reduce_element(&self, element: &Vec<Ring::Set>) -> (Vec<Ring::Set>, Vec<Ring::Set>) {
        debug_assert!(self.module().is_element(&element));
        let mut reduced_element = element.clone();
        let mut coset = vec![];
        for (r, &c) in self.pivots.iter().enumerate() {
            let quo = self
                .ring()
                .quo(&reduced_element[c], self.row_basis.at(r, c).unwrap())
                .unwrap();
            for c2 in 0..self.module_rank() {
                reduced_element[c2] = self.ring().add(
                    &reduced_element[c2],
                    &self
                        .ring()
                        .neg(&self.ring().mul(&quo, &self.row_basis.at(r, c2).unwrap())),
                );
            }
            coset.push(quo);
        }
        (coset, reduced_element)
    }

    pub fn equal_slow(x: &Self, y: &Self) -> bool {
        debug_assert_eq!(x.module(), y.module());
        Self::contains(x, y) && Self::contains(y, x)
    }

    pub fn contains_element(&self, element: &Vec<Ring::Set>) -> bool {
        debug_assert!(self.module().is_element(&element));
        let (_offset, element_reduced) = self.reduce_element(element);
        element_reduced
            .iter()
            .all(|coeff| self.ring().is_zero(coeff))
    }

    pub fn contains(x: &Self, y: &Self) -> bool {
        debug_assert_eq!(x.module(), y.module());
        for b in y.basis() {
            if !Self::contains_element(&x, &b) {
                return false;
            }
        }
        true
    }

    pub fn add(x: &Self, y: &Self) -> Self {
        debug_assert_eq!(x.module(), y.module());
        let ring = x.ring();
        debug_assert_eq!(ring, y.ring());
        let cols = x.module_rank();
        debug_assert_eq!(cols, y.module_rank());
        Self::matrix_row_span(
            ring.clone(),
            Matrix::join_rows(
                cols,
                vec![
                    x.clone().into_row_basis_matrix(),
                    y.clone().into_row_basis_matrix(),
                ],
            ),
        )
    }

    pub fn intersect(x: &Self, y: &Self) -> Self {
        debug_assert_eq!(x.module(), y.module());
        let ring = x.ring();
        debug_assert_eq!(ring, y.ring());
        let cols = x.module_rank();
        debug_assert_eq!(cols, y.module_rank());
        let x_rows = x.clone().into_row_basis_matrix();
        let y_rows = y.clone().into_row_basis_matrix();
        let matrix = Matrix::join_rows(cols, vec![&x_rows, &y_rows]);
        let matrix_ker = Self::matrix_row_kernel(ring.clone(), matrix).into_row_basis_matrix();
        let matrix_ker_first_part = matrix_ker.submatrix(
            (0..matrix_ker.rows()).collect(),
            (0..x_rows.rows()).collect(),
        );
        Self::matrix_row_span(
            ring.clone(),
            MatrixStructure::new(ring.clone())
                .mul(&matrix_ker_first_part, &x_rows)
                .unwrap(),
        )
    }

    //given a contained in b, find rank(b) - rank(a) basis vectors needed to extend a to b
    pub fn extension_basis(lat_a: &Self, lat_b: &Self) -> Vec<Vec<Ring::Set>> {
        //https://math.stackexchange.com/questions/2554408/how-to-find-the-basis-of-a-quotient-space
        let module =
            common_structure::<FinitelyFreeModuleStructure<Ring>>(lat_a.module(), lat_b.module());
        debug_assert!(Self::contains(lat_b, lat_a));

        let n = module.rank();
        // form matrix of all vectors from [other self]
        // row reduce and get pivots - take cols from orig to form basies of the quotient space

        let mat_structure = MatrixStructure::new(module.ring().clone());
        let row_span = Matrix::join_rows(n, vec![&lat_a.row_basis, &lat_b.row_basis]);
        let (_h, _u, _u_det, pivs) = mat_structure.col_hermite_algorithm(row_span.clone());

        let mut extension_basis = vec![];
        for r in pivs {
            //dont take vectors which form a basis of lat_a
            if r >= lat_a.submodule_rank() {
                extension_basis.push(row_span.get_row(r));
            }
        }

        debug_assert_eq!(
            lat_a.submodule_rank() + extension_basis.len(),
            lat_b.submodule_rank()
        );

        extension_basis
    }

    pub fn coset(&self, offset: &Vec<Ring::Set>) -> FinitelyFreeSubmoduleCoset<Ring> {
        FinitelyFreeSubmoduleCoset::from_offset_and_module(offset, self.clone())
    }

    pub fn into_coset(self) -> FinitelyFreeSubmoduleCoset<Ring> {
        FinitelyFreeSubmoduleCoset::from_offset_and_module(&self.module().zero(), self)
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature> FinitelyFreeSubmodule<Ring> {
    pub fn equal(x: &Self, y: &Self) -> bool {
        debug_assert_eq!(x.module(), y.module());
        let ring = x.ring();
        debug_assert_eq!(ring, y.ring());
        MatrixStructure::new(ring.clone()).equal(&x.row_basis, &y.row_basis)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::{Integer, Rational};

    #[test]
    fn test_finitely_free_submodule_reduction() {
        let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
        a.pprint();
        let a_reduced =
            FinitelyFreeSubmodule::matrix_row_span(Integer::structure(), a).into_row_basis_matrix();
        let a_expected = Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]);
        a_reduced.pprint();
        a_expected.pprint();
        assert_eq!(a_reduced, a_expected);
    }

    #[test]
    fn test_finitely_free_submodule_unreduced_equal() {
        assert!(FinitelyFreeSubmodule::equal(
            &FinitelyFreeSubmodule::matrix_row_span(
                Integer::structure(),
                Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]])
            ),
            &FinitelyFreeSubmodule::matrix_row_span(
                Integer::structure(),
                Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]])
            )
        ))
    }

    #[test]
    fn test_finitely_free_submodule_intersect() {
        let a = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![
                vec![2, 0, 0, 0],
                vec![0, 2, 0, 0],
                vec![0, 0, 2, 0],
                vec![0, 0, 0, 2],
            ]),
        );

        let b = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![3, 3, 3, 3], vec![3, 3, -3, -3]]),
        );

        let c = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![6, 6, 0, 0], vec![0, 0, 6, 6]]),
        );

        let s = FinitelyFreeSubmodule::intersect(&a, &b);

        s.clone().into_row_basis_matrix().pprint();

        assert!(FinitelyFreeSubmodule::equal(&c, &s));
    }

    #[test]
    fn test_finitely_free_submodule_element_reduction() {
        let a = FinitelyFreeSubmodule::matrix_row_span(
            Integer::structure(),
            Matrix::from_rows(vec![vec![3, 2, 0, 3], vec![0, 14, 3, 1], vec![0, 0, 0, 10]]),
        );
        a.clone().into_row_basis_matrix().pprint();

        let element = vec![20, 20, 20, 20]
            .into_iter()
            .map(|x| Integer::from(x))
            .collect::<Vec<_>>();
        println!("element = {:?}", element);

        let (offset, reduced_element) = a.reduce_element(&element);

        println!("offset = {:?}", offset);
        println!("reduced_element = {:?}", reduced_element);
        debug_assert_eq!(
            reduced_element,
            vec![2, 8, 20, 2]
                .into_iter()
                .map(|x| Integer::from(x))
                .collect::<Vec<_>>()
        )
    }

    #[test]
    fn test_finitely_free_submodule_extension_basis() {
        let a = Matrix::<Rational>::from_rows(vec![vec![1, 0, 0], vec![1, 0, 0], vec![-1, 0, 0]]);

        let b = Matrix::<Rational>::from_rows(vec![vec![1, 1, 0], vec![1, 1, 0], vec![1, -1, 0]]);

        println!("a");
        a.pprint();
        println!("b");
        b.pprint();

        let ext = FinitelyFreeSubmodule::extension_basis(&a.col_span(), &b.col_span());

        println!("ext");
        for v in ext {
            println!("{:?}", v);
        }
    }
}
