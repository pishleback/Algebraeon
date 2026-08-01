use super::{finitely_free_cosets::*, finitely_free_module::*};
use crate::{matrix::*, structure::*};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use std::{borrow::Borrow, fmt::Debug, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmodule<Elem: Clone + Debug> {
    // a matrix in reduced hermite normal form with all non-zero rows whose rows form a basis for the submodule
    row_basis: Matrix<Elem>,
    // the columns of the pivots of row_basis
    pivots: Vec<usize>,
}

impl<Elem: Clone + Debug> FinitelyFreeSubmodule<Elem> {
    pub fn into_row_basis_matrix(self) -> Matrix<Elem> {
        self.row_basis
    }

    pub fn into_col_basis_matrix(self) -> Matrix<Elem> {
        self.row_basis.transpose()
    }

    pub(crate) fn row_basis_matrix(&self) -> &Matrix<Elem> {
        &self.row_basis
    }

    pub fn basis(&self) -> Vec<Vec<Elem>> {
        (0..self.row_basis.rows())
            .map(|r| {
                (0..self.row_basis.cols())
                    .map(|c| self.row_basis.at(r, c).unwrap().clone())
                    .collect()
            })
            .collect()
    }

    pub fn rank(&self) -> usize {
        self.row_basis.rows()
    }

    pub fn module_rank(&self) -> usize {
        self.row_basis.cols()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmodulesStructure<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> {
    _set: PhantomData<Set>,
    _ring: PhantomData<Ring>,
    _module: PhantomData<Module>,
    module: ModuleB,
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> FinitelyFreeSubmodulesStructure<Set, Ring, Module, ModuleB>
{
    pub fn new(module: ModuleB) -> Self {
        Self {
            _set: PhantomData,
            _ring: PhantomData,
            _module: PhantomData,
            module,
        }
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> Signature for FinitelyFreeSubmodulesStructure<Set, Ring, Module, ModuleB>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> SetSignature for FinitelyFreeSubmodulesStructure<Set, Ring, Module, ModuleB>
{
    type Elem = FinitelyFreeSubmodule<Ring::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if x.row_basis.cols() != self.module().rank() {
            return Err("dimensions don't match".to_string());
        }
        // TODO: more checks
        Ok(())
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> FinitelyFreeSubmodulesStructure<Set, Ring, Module, ModuleB>
{
    pub fn module(&self) -> &Module {
        self.module.borrow()
    }

    pub fn ring(&self) -> &Ring {
        self.module().ring()
    }

    pub fn matrix_row_span_and_basis(
        &self,
        matrix: Matrix<Ring::Elem>,
    ) -> (FinitelyFreeSubmodule<Ring::Elem>, Matrix<Ring::Elem>) {
        let ring = self.ring();
        let rows = matrix.rows();
        let cols = matrix.cols();
        debug_assert_eq!(cols, self.module().rank());
        let (h, u, _det, pivots) =
            MatrixStructure::new(ring.clone()).row_reduced_hermite_algorithm(matrix);
        let row_basis = h.submatrix((0..pivots.len()).collect(), (0..cols).collect());
        let u = u.submatrix((0..pivots.len()).collect(), (0..rows).collect());
        (FinitelyFreeSubmodule { row_basis, pivots }, u)
    }

    pub fn matrix_col_span_and_basis(
        &self,
        matrix: Matrix<Ring::Elem>,
    ) -> (FinitelyFreeSubmodule<Ring::Elem>, Matrix<Ring::Elem>) {
        let (s, u) = self.matrix_row_span_and_basis(matrix.transpose());
        (s, u.transpose())
    }

    pub fn matrix_row_span(&self, matrix: Matrix<Ring::Elem>) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.matrix_row_span_and_basis(matrix).0
    }

    pub fn matrix_col_span(&self, matrix: Matrix<Ring::Elem>) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.matrix_col_span_and_basis(matrix).0
    }

    pub fn span(&self, span: Vec<impl Borrow<Module::Elem>>) -> FinitelyFreeSubmodule<Ring::Elem> {
        for v in &span {
            debug_assert!(self.module().is_element(v.borrow()));
        }
        let span = span
            .into_iter()
            .map(|v| self.module().to_vec(v.borrow()))
            .collect::<Vec<_>>();
        self.matrix_row_span(Matrix::construct(
            span.len(),
            self.module().rank(),
            |r, c| span[r][c].clone(),
        ))
    }

    pub fn kernel(
        &self,
        items: Vec<impl Borrow<Vec<Ring::Elem>>>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        let n = self.module().rank();
        debug_assert_eq!(items.len(), n);
        if n == 0 {
            self.zero_submodule()
        } else {
            let cols = items.first().unwrap().borrow().len();
            for v in &items[1..] {
                assert_eq!(v.borrow().len(), cols);
            }
            self.matrix_row_kernel(Matrix::construct(n, cols, |r, c| {
                items[r].borrow()[c].clone()
            }))
        }
    }

    pub fn matrix_row_preimage(
        &self,
        matrix: &Matrix<Ring::Elem>,
        space: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        debug_assert_eq!(matrix.rows(), self.module().rank());
        debug_assert_eq!(matrix.cols(), space.module_rank());
        /*
        Let M be the matrix and S a row basis matrix for the space.
        Need to solve xM=yS for all possible x
        So concat the rows of M over the rows of S to get a matrix A
        Solve for the kernel of A and take only the first bits corresponding to x
        */
        let a = Matrix::join_rows(matrix.cols(), vec![&matrix, space.row_basis_matrix()]);
        // Find a row U such that uA is in HNF
        // The kernel of UA is standard basis vectors pivs.len()..a.rows()
        // So the kernel of A is the rows pivs.len()..a.rows() of U
        // And we are after the span of the first matrix.rows() of them
        let (h, u, _u_det, pivs) =
            MatrixStructure::<Ring, _>::new(self.ring()).row_hermite_algorithm(a);

        self.matrix_row_span(u.submatrix(
            (pivs.len()..h.rows()).collect(),
            (0..matrix.rows()).collect(),
        ))
    }

    pub fn matrix_col_preimage(
        &self,
        matrix: &Matrix<Ring::Elem>,
        space: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.matrix_row_preimage(&matrix.transpose_ref(), space)
    }

    pub fn matrix_row_kernel(
        &self,
        matrix: Matrix<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        debug_assert_eq!(matrix.rows(), self.module().rank());
        let rows = matrix.rows();
        let (_h, u, _u_det, pivs) =
            MatrixStructure::<Ring, _>::new(self.ring()).row_hermite_algorithm(matrix);
        debug_assert_eq!(rows, u.rows());
        debug_assert_eq!(rows, u.cols());
        let ker = u.submatrix((pivs.len()..rows).collect(), (0..rows).collect());
        self.matrix_row_span(ker)
    }

    pub fn matrix_col_kernel(
        &self,
        matrix: Matrix<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        debug_assert_eq!(matrix.cols(), self.module().rank());
        self.matrix_row_kernel(matrix.transpose())
    }

    pub fn full_submodule(&self) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.matrix_row_span(
            MatrixStructure::<Ring, _>::new(self.ring().clone()).ident(self.module().rank()),
        )
    }

    pub fn zero_submodule(&self) -> FinitelyFreeSubmodule<Ring::Elem> {
        let cols = self.module().rank();
        FinitelyFreeSubmodule {
            row_basis: Matrix::construct(0, cols, |_, _| unreachable!()),
            pivots: vec![],
        }
    }

    pub fn reduce_element(
        &self,
        submodule: &FinitelyFreeSubmodule<Ring::Elem>,
        element: &Module::Elem,
    ) -> (Vec<Ring::Elem>, Module::Elem) {
        debug_assert!(self.validate_element(submodule).is_ok());
        debug_assert!(self.module().validate_element(element).is_ok());
        let mut reduced_element = self.module().to_vec(element);
        let mut coset = vec![];
        for (r, &c) in submodule.pivots.iter().enumerate() {
            let quo = self
                .ring()
                .quo(&reduced_element[c], submodule.row_basis.at(r, c).unwrap())
                .unwrap();
            #[allow(clippy::needless_range_loop)]
            for c2 in 0..self.module().rank() {
                reduced_element[c2] = self.ring().add(
                    &reduced_element[c2],
                    &self.ring().neg(
                        &self
                            .ring()
                            .mul(&quo, submodule.row_basis.at(r, c2).unwrap()),
                    ),
                );
            }
            coset.push(quo);
        }

        (coset, self.module().from_vec(reduced_element))
    }

    pub fn equal_slow(
        &self,
        x: &FinitelyFreeSubmodule<Ring::Elem>,
        y: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        self.contains(x, y) && self.contains(y, x)
    }

    pub fn contains_element(
        &self,
        submodule: &FinitelyFreeSubmodule<Ring::Elem>,
        element: &Module::Elem,
    ) -> bool {
        debug_assert!(self.validate_element(submodule).is_ok());
        debug_assert!(self.module().validate_element(element).is_ok());
        let (_offset, element_reduced) = self.reduce_element(submodule, element);
        self.module()
            .to_vec(&element_reduced)
            .into_iter()
            .all(|coeff| self.ring().is_zero(&coeff))
    }

    pub fn contains(
        &self,
        x: &FinitelyFreeSubmodule<Ring::Elem>,
        y: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        for b in y.basis() {
            if !self.contains_element(x, &self.module().from_vec(b)) {
                return false;
            }
        }
        true
    }

    pub fn add(
        &self,
        x: FinitelyFreeSubmodule<Ring::Elem>,
        y: FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        self.sum(vec![x, y])
    }

    pub fn sum(
        &self,
        xs: Vec<FinitelyFreeSubmodule<Ring::Elem>>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        for x in &xs {
            debug_assert!(self.validate_element(x).is_ok());
        }
        let cols = self.module().rank();
        self.matrix_row_span(Matrix::join_rows(
            cols,
            xs.into_iter().map(|x| x.into_row_basis_matrix()).collect(),
        ))
    }

    pub fn intersect(
        &self,
        x: FinitelyFreeSubmodule<Ring::Elem>,
        y: FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        debug_assert!(self.validate_element(&x).is_ok());
        debug_assert!(self.validate_element(&y).is_ok());

        let cols = self.module().rank();
        let x_rows = x.into_row_basis_matrix();
        let y_rows = y.into_row_basis_matrix();
        let matrix = Matrix::join_rows(cols, vec![&x_rows, &y_rows]);
        let matrix_ker = self
            .ring()
            .free_module(EnumeratedFiniteSetStructure::new(matrix.rows()))
            .submodules()
            .matrix_row_kernel(matrix)
            .into_row_basis_matrix();
        let matrix_ker_first_part = matrix_ker.submatrix(
            (0..matrix_ker.rows()).collect(),
            (0..x_rows.rows()).collect(),
        );
        self.matrix_row_span(
            MatrixStructure::<Ring, _>::new(self.ring())
                .mul(&matrix_ker_first_part, &x_rows)
                .unwrap(),
        )
    }

    pub fn intersect_list(
        &self,
        mut xs: Vec<FinitelyFreeSubmodule<Ring::Elem>>,
    ) -> FinitelyFreeSubmodule<Ring::Elem> {
        if let Some(a) = xs.pop() {
            if let Some(b) = xs.pop() {
                let mut i = self.intersect(a, b);
                for c in xs {
                    i = self.intersect(i, c);
                }
                i
            } else {
                a
            }
        } else {
            self.full_submodule()
        }
    }

    //given x contained in y, find rank(y) - rank(x) basis vectors needed to extend x to y
    pub fn extension_basis(
        &self,
        x: &FinitelyFreeSubmodule<Ring::Elem>,
        y: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> Vec<Vec<Ring::Elem>> {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());

        //https://math.stackexchange.com/questions/2554408/how-to-find-the-basis-of-a-quotient-space
        debug_assert!(self.contains(y, x));

        let n = self.module().rank();
        // form matrix of all vectors from [other self]
        // row reduce and get pivots - take cols from orig to form basies of the quotient space

        let mat_structure = MatrixStructure::<Ring, _>::new(self.ring());
        let row_span = Matrix::join_rows(n, vec![&x.row_basis, &y.row_basis]);
        let (_h, _u, _u_det, pivs) = mat_structure.col_hermite_algorithm(row_span.clone());

        let mut extension_basis = vec![];
        for r in pivs {
            //dont take vectors which form a basis of lat_a
            if r >= x.rank() {
                extension_basis.push(row_span.get_row(r));
            }
        }

        debug_assert_eq!(x.rank() + extension_basis.len(), y.rank());

        extension_basis
    }

    pub fn coset(
        &self,
        x: &FinitelyFreeSubmodule<Ring::Elem>,
        offset: &Module::Elem,
    ) -> FinitelyFreeSubmoduleCoset<Ring::Elem> {
        self.module()
            .cosets()
            .from_offset_and_submodule(offset, x.clone())
    }

    pub fn into_coset(
        &self,
        x: FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmoduleCoset<Ring::Elem> {
        self.module()
            .cosets()
            .from_offset_and_submodule(&self.module().zero(), x)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> EqSignature for FinitelyFreeSubmodulesStructure<Set, Ring, Module, ModuleB>
{
    fn equal(
        &self,
        x: &FinitelyFreeSubmodule<Ring::Elem>,
        y: &FinitelyFreeSubmodule<Ring::Elem>,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        MatrixStructure::<Ring, _>::new(self.ring()).equal(&x.row_basis, &y.row_basis)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_finitely_free_submodule_kernel() {
        let submodules = Integer::structure()
            .into_free_module(EnumeratedFiniteSetStructure::new(3))
            .into_submodules();

        let a = submodules.span(vec![&vec![1.into(), 1.into(), (-1).into()]]);

        let b = submodules.kernel(vec![
            &vec![1.into(), 1.into(), 2.into(), 2.into()],
            &vec![2.into(), 2.into(), 1.into(), 1.into()],
            &vec![3.into(), 3.into(), 3.into(), 3.into()],
        ]);

        assert!(submodules.equal(&a, &b));
    }

    #[test]
    fn test_finitely_free_submodule_reduction() {
        let a = Matrix::from_rows(vec![vec![1, 2, 4, 5], vec![1, 2, 3, 4]]);
        a.pprint();
        let a_reduced = Integer::structure()
            .free_module(EnumeratedFiniteSetStructure::new(4))
            .submodules()
            .matrix_row_span(a)
            .into_row_basis_matrix();
        let a_expected = Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]);
        a_reduced.pprint();
        a_expected.pprint();
        assert_eq!(a_reduced, a_expected);
    }

    #[test]
    fn test_finitely_free_submodule_unreduced_equal() {
        let modules = Integer::structure().into_free_module(EnumeratedFiniteSetStructure::new(4));
        assert!(
            modules.submodules().equal(
                &modules
                    .submodules()
                    .matrix_row_span(Matrix::from_rows(vec![vec![1, 2, 3, 4], vec![0, 0, 1, 1]])),
                &modules
                    .submodules()
                    .matrix_row_span(Matrix::from_rows(vec![vec![1, 2, 0, 1], vec![0, 0, 1, 1]]))
            )
        );
    }

    #[test]
    fn test_finitely_free_submodule_intersect() {
        let modules = Integer::structure().into_free_module(EnumeratedFiniteSetStructure::new(4));

        let a = modules.submodules().matrix_row_span(Matrix::from_rows(vec![
            vec![2, 0, 0, 0],
            vec![0, 2, 0, 0],
            vec![0, 0, 2, 0],
            vec![0, 0, 0, 2],
        ]));

        let b = modules.submodules().matrix_row_span(Matrix::from_rows(vec![
            vec![3, 3, 3, 3],
            vec![3, 3, -3, -3],
        ]));

        let c = modules
            .submodules()
            .matrix_row_span(Matrix::from_rows(vec![vec![6, 6, 0, 0], vec![0, 0, 6, 6]]));

        let s = modules.submodules().intersect(a, b);

        s.clone().into_row_basis_matrix().pprint();

        assert!(modules.submodules().equal(&c, &s));
    }

    #[test]
    fn test_finitely_free_submodule_element_reduction() {
        let modules = Integer::structure().into_free_module(EnumeratedFiniteSetStructure::new(4));

        let a = modules.submodules().matrix_row_span(Matrix::from_rows(vec![
            vec![3, 2, 0, 3],
            vec![0, 14, 3, 1],
            vec![0, 0, 0, 10],
        ]));
        a.clone().into_row_basis_matrix().pprint();

        let element = vec![20, 20, 20, 20]
            .into_iter()
            .map(Integer::from)
            .collect::<Vec<_>>();
        println!("element = {:?}", element);

        let (offset, reduced_element) = modules.submodules().reduce_element(&a, &element);

        println!("offset = {:?}", offset);
        println!("reduced_element = {:?}", reduced_element);
        debug_assert_eq!(
            reduced_element,
            vec![2, 8, 20, 2]
                .into_iter()
                .map(Integer::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_finitely_free_submodule_extension_basis() {
        let modules = Rational::structure().into_free_module(EnumeratedFiniteSetStructure::new(3));

        let a = Matrix::<Rational>::from_rows(vec![vec![1, 0, 0], vec![1, 0, 0], vec![-1, 0, 0]]);

        let b = Matrix::<Rational>::from_rows(vec![vec![1, 1, 0], vec![1, 1, 0], vec![1, -1, 0]]);

        println!("a");
        a.pprint();
        println!("b");
        b.pprint();

        let ext = modules
            .submodules()
            .extension_basis(&a.col_span(), &b.col_span());

        println!("ext");
        for v in ext {
            println!("{:?}", v);
        }
    }
}
