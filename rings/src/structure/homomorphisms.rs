use super::*;
use crate::matrix::{Matrix, MatrixStructure};
use crate::polynomial::Polynomial;
use algebraeon_structures::*;
use std::borrow::Borrow;
use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::rc::Rc;

mod range_module {
    use std::rc::Rc;

    use super::*;
    #[derive(Debug, Clone)]
    pub struct RingHomomorphismRangeModuleStructure<
        Domain: RingSignature,
        Range: RingSignature,
        Hom: RingHomomorphism<Domain, Range>,
    > {
        _domain: PhantomData<Domain>,
        _range: PhantomData<Range>,
        hom: Rc<Hom>,
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn new(hom: Rc<Hom>) -> Rc<Self> {
            Self {
                _domain: PhantomData,
                _range: PhantomData,
                hom,
            }
            .into()
        }

        pub fn module(&self) -> Rc<Range> {
            self.hom.range()
        }

        pub fn homomorphism(&self) -> &Rc<Hom> {
            &self.hom
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        PartialEq for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn eq(&self, other: &Self) -> bool {
            std::ptr::eq(self.hom.as_ref(), other.hom.as_ref())
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>> Eq
        for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        Signature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        RinglikeSpecializationSignature
        for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        SetSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        type Elem = Range::Elem;

        fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
            self.hom.range().validate_element(x)
        }
    }

    impl<Domain: RingSignature, Range: RingEqSignature, Hom: RingHomomorphism<Domain, Range>>
        EqSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
            self.hom.range().equal(a, b)
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        ZeroSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn zero(self: &Rc<Self>) -> Self::Elem {
            self.hom.range().zero()
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        AdditionSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn add(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
            self.hom.range().add(a, b)
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        CancellativeAdditionSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn try_sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
            Some(self.sub(a, b))
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        TryNegateSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn try_neg(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
            Some(self.neg(a))
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        AdditiveMonoidSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        AdditiveGroupSignature for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn neg(self: &Rc<Self>, a: &Self::Elem) -> Self::Elem {
            self.hom.range().neg(a)
        }
    }

    impl<Domain: RingSignature, Range: RingSignature, Hom: RingHomomorphism<Domain, Range>>
        SemiModuleSignature<Domain> for RingHomomorphismRangeModuleStructure<Domain, Range, Hom>
    {
        fn ring(self: &Rc<Self>) -> Rc<Domain> {
            self.hom.domain()
        }

        fn scalar_mul(self: &Rc<Self>, a: &Self::Elem, x: &Domain::Elem) -> Self::Elem {
            self.hom.range().mul(&self.hom.image(x), a)
        }
    }

    pub trait RingHomomorphism<Domain: RingSignature, Range: RingSignature>:
        FunctionMorphism<Domain, Range>
    {
        fn range_module_structure(
            self: &Rc<Self>,
        ) -> Rc<RingHomomorphismRangeModuleStructure<Domain, Range, Self>> {
            RingHomomorphismRangeModuleStructure::new(self.clone())
        }
    }

    impl<
        A: RingSignature,
        B: RingSignature,
        C: RingSignature,
        AB: RingHomomorphism<A, B>,
        BC: RingHomomorphism<B, C>,
    > RingHomomorphism<A, C> for CompositionMorphism<A, B, C, AB, BC>
    {
    }
}
pub use range_module::*;

mod principal_subring_inclusion {
    use super::*;
    use std::rc::Rc;

    /// The unique ring homomorphism Z -> R of the integers into any ring R
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct PrincipalIntegerMap<Ring: RingSignature> {
        ring: Rc<Ring>,
    }

    impl<Ring: RingSignature> PrincipalIntegerMap<Ring> {
        pub fn new(ring: Rc<Ring>) -> Rc<Self> {
            Self { ring }.into()
        }
    }

    impl<Ring: RingSignature> Morphism<IntegerCanonicalStructure, Ring> for PrincipalIntegerMap<Ring> {
        fn domain(self: &Rc<Self>) -> Rc<IntegerCanonicalStructure> {
            Integer::structure()
        }

        fn range(self: &Rc<Self>) -> Rc<Ring> {
            self.ring.clone()
        }
    }

    impl<Ring: RingSignature> FunctionMorphism<IntegerCanonicalStructure, Ring>
        for PrincipalIntegerMap<Ring>
    {
        fn image(self: &Rc<Self>, x: &Integer) -> <Ring as SetSignature>::Elem {
            self.range().from_int(x)
        }
    }

    impl<Ring: CharZeroRingSignature> InjectiveFunctionMorphism<IntegerCanonicalStructure, Ring>
        for PrincipalIntegerMap<Ring>
    {
        fn try_preimage(self: &Rc<Self>, x: &<Ring as SetSignature>::Elem) -> Option<Integer> {
            self.range().try_to_int(x)
        }
    }

    impl<Ring: RingSignature> RingHomomorphism<IntegerCanonicalStructure, Ring>
        for PrincipalIntegerMap<Ring>
    {
    }

    /// The unique field embedding Q -> K of the rationals into any field of characteristic zero
    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct PrincipalRationalMap<Field: CharZeroFieldSignature> {
        field: Rc<Field>,
    }

    impl<Field: CharZeroFieldSignature> PrincipalRationalMap<Field> {
        pub fn new(field: Rc<Field>) -> Rc<Self> {
            Self { field }.into()
        }
    }

    impl<Field: CharZeroFieldSignature> Morphism<RationalCanonicalStructure, Field>
        for PrincipalRationalMap<Field>
    {
        fn domain(self: &Rc<Self>) -> Rc<RationalCanonicalStructure> {
            Rational::structure()
        }

        fn range(self: &Rc<Self>) -> Rc<Field> {
            self.field.clone()
        }
    }

    impl<Field: CharZeroFieldSignature> FunctionMorphism<RationalCanonicalStructure, Field>
        for PrincipalRationalMap<Field>
    {
        fn image(self: &Rc<Self>, x: &Rational) -> <Field as SetSignature>::Elem {
            self.range().try_from_rat(x).unwrap()
        }
    }

    impl<Field: CharZeroFieldSignature> InjectiveFunctionMorphism<RationalCanonicalStructure, Field>
        for PrincipalRationalMap<Field>
    {
        fn try_preimage(self: &Rc<Self>, x: &<Field as SetSignature>::Elem) -> Option<Rational> {
            self.range().try_to_rat(x)
        }
    }

    impl<Field: CharZeroFieldSignature> RingHomomorphism<RationalCanonicalStructure, Field>
        for PrincipalRationalMap<Field>
    {
    }
}
pub use principal_subring_inclusion::*;

/// The inclusion of an integral domain into its field of fractions
pub trait FieldOfFractionsInclusion<Ring: RingSignature, Field: FieldSignature>:
    RingHomomorphism<Ring, Field> + InjectiveFunctionMorphism<Ring, Field>
{
    fn numerator_and_denominator(self: &Rc<Self>, a: &Field::Elem) -> (Ring::Elem, Ring::Elem);
    fn numerator(self: &Rc<Self>, a: &Field::Elem) -> Ring::Elem {
        self.numerator_and_denominator(a).0
    }
    fn denominator(self: &Rc<Self>, a: &Field::Elem) -> Ring::Elem {
        self.numerator_and_denominator(a).1
    }
}

/// An injective homomorphism A -> B of integral domains where there is a way to get all roots in B of a polynomial over A
pub trait IntegralDomainExtensionAllPolynomialRoots<
    A: IntegralDomainSignature,
    B: IntegralDomainSignature,
>: RingHomomorphism<A, B> + InjectiveFunctionMorphism<A, B>
{
    /// Return all roots of the polynomial in B with duplicate elements according to multiplicity
    fn all_roots(self: &Rc<Self>, polynomial: &Polynomial<A::Elem>) -> Vec<B::Elem>;
}

/// A ring extension Z -> R such that R is a finitely free Z-module
pub trait FiniteRankFreeRingExtension<Basis: FiniteSetSignature, Z: RingSignature, R: RingSignature>:
    RingHomomorphism<Z, R> + InjectiveFunctionMorphism<Z, R>
{
    // things inherited from the finitely free domain-module structure on the range
    fn degree(self: &Rc<Self>) -> usize;
    fn to_vec(self: &Rc<Self>, a: &R::Elem) -> Vec<Z::Elem>;
    fn from_vec(self: &Rc<Self>, v: Vec<impl Borrow<Z::Elem>>) -> R::Elem;
    fn to_col(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem>;
    fn from_col(self: &Rc<Self>, v: Matrix<Z::Elem>) -> R::Elem;
    fn to_row(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem>;
    fn from_row(self: &Rc<Self>, v: Matrix<Z::Elem>) -> R::Elem;

    /// matrix representing column vector multiplication by `a` on the left
    fn col_multiplication_matrix(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem>;

    /// matrix representing row vector multiplication by `a` on the right
    fn row_multiplication_matrix(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem>;
}

impl<
    Basis: FiniteSetSignature,
    Z: RingSignature,
    R: RingSignature,
    Hom: RingHomomorphism<Z, R> + InjectiveFunctionMorphism<Z, R>,
> FiniteRankFreeRingExtension<Basis, Z, R> for Hom
where
    for<'h> RingHomomorphismRangeModuleStructure<Z, R, Self>:
        FinitelyFreeModuleSignature<Basis, Z, Elem = R::Elem>,
{
    fn degree(self: &Rc<Self>) -> usize {
        self.range_module_structure().rank()
    }
    fn to_vec(self: &Rc<Self>, a: &R::Elem) -> Vec<Z::Elem> {
        self.range_module_structure().to_vec(a)
    }
    fn from_vec(self: &Rc<Self>, v: Vec<impl Borrow<Z::Elem>>) -> R::Elem {
        self.range_module_structure().from_vec(v)
    }
    fn to_col(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem> {
        self.range_module_structure().to_col(a)
    }
    fn from_col(self: &Rc<Self>, v: Matrix<Z::Elem>) -> R::Elem {
        self.range_module_structure().from_col(v)
    }
    fn to_row(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem> {
        self.range_module_structure().to_row(a)
    }
    fn from_row(self: &Rc<Self>, v: Matrix<Z::Elem>) -> R::Elem {
        self.range_module_structure().from_row(v)
    }

    fn col_multiplication_matrix(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem> {
        let basis = self.range_module_structure().basis_vecs();
        Matrix::from_cols(
            (0..self.degree())
                .map(|i| {
                    self.range_module_structure()
                        .to_vec(&self.range().mul(a, &basis[i]))
                })
                .collect(),
        )
    }

    fn row_multiplication_matrix(self: &Rc<Self>, a: &R::Elem) -> Matrix<Z::Elem> {
        self.col_multiplication_matrix(a).transpose()
    }
}

/// A finite dimensional field extension F -> K
pub trait FiniteDimensionalFieldExtension<
    Basis: FiniteSetSignature,
    F: FieldSignature,
    K: FieldSignature,
>:
    RingHomomorphism<F, K> + InjectiveFunctionMorphism<F, K> + FiniteRankFreeRingExtension<Basis, F, K>
{
    fn norm(self: &Rc<Self>, a: &K::Elem) -> F::Elem;

    fn trace(self: &Rc<Self>, a: &K::Elem) -> F::Elem;

    /// The monic minimal polynomial of a
    fn min_poly(self: &Rc<Self>, a: &K::Elem) -> Polynomial<F::Elem>;

    fn trace_form_matrix(self: &Rc<Self>, elems: &[K::Elem]) -> Matrix<F::Elem> {
        let n = self.degree();
        assert_eq!(n, elems.len());
        Matrix::construct(n, n, |r, c| {
            self.trace(&self.range().mul(&elems[r], &elems[c]))
        })
    }

    fn discriminant(self: &Rc<Self>, elems: &[K::Elem]) -> F::Elem {
        MatrixStructure::new(self.domain().clone())
            .det(self.trace_form_matrix(elems))
            .unwrap()
    }
}

impl<
    Basis: FiniteSetSignature,
    F: FieldSignature,
    K: FieldSignature,
    Hom: RingHomomorphism<F, K> + InjectiveFunctionMorphism<F, K>,
> FiniteDimensionalFieldExtension<Basis, F, K> for Hom
where
    for<'h> RingHomomorphismRangeModuleStructure<F, K, Self>:
        FinitelyFreeModuleSignature<Basis, F, Elem = K::Elem>,
{
    fn norm(self: &Rc<Self>, a: &K::Elem) -> F::Elem {
        MatrixStructure::new(self.domain().clone())
            .det(self.col_multiplication_matrix(a))
            .unwrap()
    }

    fn trace(self: &Rc<Self>, a: &K::Elem) -> F::Elem {
        MatrixStructure::new(self.domain().clone())
            .trace(&self.col_multiplication_matrix(a))
            .unwrap()
    }

    fn min_poly(self: &Rc<Self>, a: &K::Elem) -> Polynomial<F::Elem> {
        MatrixStructure::new(self.domain().clone())
            .minimal_polynomial(self.col_multiplication_matrix(a))
            .unwrap()
    }
}

/// Represent all ring homomorphisms from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingHomomorphisms<Domain: RingSignature, Range: RingSignature> {
    domain: Rc<Domain>,
    range: Rc<Range>,
}

impl<Domain: RingSignature, Range: RingSignature> RingHomomorphisms<Domain, Range> {
    pub fn new(domain: Rc<Domain>, range: Rc<Range>) -> Self {
        Self { domain, range }
    }
}

impl<Domain: RingSignature, Range: RingSignature> Signature for RingHomomorphisms<Domain, Range> {}

impl<Domain: FreeRingSignature, Range: RingSignature> SetSignature
    for RingHomomorphisms<Domain, Range>
{
    type Elem = HashMap<Domain::Generator, Range::Elem>;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if self.domain.free_generators() != x.keys().cloned().collect::<HashSet<_>>() {
            return Err("missing key".to_string());
        }

        for v in x.values() {
            self.range.validate_element(v)?;
        }

        Ok(())
    }
}
