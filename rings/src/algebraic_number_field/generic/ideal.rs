use crate::algebraic_number_field::{
    FullRankSublatticeWithBasisSignature, RingOfIntegersIntegralExtensionWithIdeals,
    RingOfIntegersToAlgebraicNumberFieldInclusion,
};
use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, OrderWithBasis,
    },
    integer::ideal::IntegerIdealsStructure,
    matrix::Matrix,
    module::{
        finitely_free_affine::FinitelyFreeSubmoduleAffineSubset,
        finitely_free_submodule::{FinitelyFreeSubmodule, FinitelyFreeSubmoduleStructure},
    },
    polynomial::ToPolynomialSignature,
    structure::*,
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Natural, traits::Abs};
use algebraeon_sets::{
    combinatorics::num_partitions_part_pool,
    structure::{
        BorrowedMorphism, BorrowedStructure, EqSignature, Function, InjectiveFunction, MetaType,
        SetSignature, Signature,
    },
};
use itertools::Itertools;
use std::{borrow::Cow, marker::PhantomData};

#[derive(Debug, Clone)]
pub enum OrderIdeal {
    Zero,
    NonZero(FinitelyFreeSubmodule<Integer>),
}

impl OrderIdeal {
    /// A basis of this ideal as a Z-module.
    pub fn basis(&self) -> Option<Vec<Vec<Integer>>> {
        match self {
            OrderIdeal::Zero => None,
            OrderIdeal::NonZero(sublattice) => Some(sublattice.basis()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrderIdealsStructure<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> {
    _k: PhantomData<K>,
    _kb: PhantomData<KB>,
    order: OB,
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    pub fn order(&self) -> &OrderWithBasis<K, KB, MAXIMAL> {
        self.order.borrow()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    RingToIdealsSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    type Ideals<SelfB: BorrowedStructure<Self>> = OrderIdealsStructure<K, KB, MAXIMAL, SelfB>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        OrderIdealsStructure {
            _k: PhantomData,
            _kb: PhantomData,
            order: self,
        }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        OrderIdealsStructure {
            _k: PhantomData,
            _kb: PhantomData,
            order: self,
        }
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> Signature for OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    fn does_sublattice_define_an_ideal(&self, sublattice: &FinitelyFreeSubmodule<Integer>) -> bool {
        #[cfg(debug_assertions)]
        self.order()
            .free_lattice_restructure()
            .into_submodules()
            .is_element(sublattice)
            .unwrap();

        // check it's an ideal
        for ideal_basis_elem in sublattice.basis() {
            for ring_basis_elem in (0..self.order().n()).map(|i| {
                self.order
                    .borrow()
                    .free_lattice_restructure()
                    .basis_element(i)
            }) {
                let x = self.order().mul(&ideal_basis_elem, &ring_basis_elem);
                if !self
                    .order
                    .borrow()
                    .free_lattice_restructure()
                    .submodules()
                    .contains_element(sublattice, &x)
                {
                    return false;
                }
            }
        }
        // All ideals should have full rank in an order
        debug_assert_eq!(sublattice.rank(), self.order().n());
        true
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> SetSignature for OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    type Set = OrderIdeal;

    fn is_element(&self, ideal: &Self::Set) -> Result<(), String> {
        match ideal {
            OrderIdeal::Zero => Ok(()),
            OrderIdeal::NonZero(sublattice) => {
                // check it's a sublattice
                self.order
                    .borrow()
                    .free_lattice_restructure()
                    .submodules()
                    .is_element(sublattice)?;
                // check it's an ideal
                if !self.does_sublattice_define_an_ideal(sublattice) {
                    return Err("sublattice is not an ideal".to_string());
                }
                Ok(())
            }
        }
    }
}

mod sublattices_to_ideals {
    use super::*;
    use algebraeon_sets::structure::Morphism;

    #[derive(Debug, Clone)]
    pub struct SubmoduleToIdeals<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
        IB: BorrowedStructure<OrderIdealsStructure<K, KB, MAXIMAL, OB>>,
        IntB: BorrowedStructure<IntegerCanonicalStructure>,
        SB: BorrowedStructure<FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>>,
    > {
        _k: PhantomData<K>,
        _kb: PhantomData<KB>,
        _ob: PhantomData<OB>,
        _intb: PhantomData<IntB>,
        ideals: IB,
        sublattices: SB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
        IB: BorrowedStructure<OrderIdealsStructure<K, KB, MAXIMAL, OB>>,
        IntB: BorrowedStructure<IntegerCanonicalStructure>,
        SB: BorrowedStructure<FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>>,
    > SubmoduleToIdeals<K, KB, MAXIMAL, OB, IB, IntB, SB>
    {
        pub fn new(ideals: IB, sublattices: SB) -> Self {
            Self {
                _k: PhantomData,
                _kb: PhantomData,
                _ob: PhantomData,
                _intb: PhantomData,
                ideals,
                sublattices,
            }
        }

        fn ideals(&self) -> &OrderIdealsStructure<K, KB, MAXIMAL, OB> {
            self.ideals.borrow()
        }

        fn sublattices(&self) -> &FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB> {
            self.sublattices.borrow()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
        IB: BorrowedStructure<OrderIdealsStructure<K, KB, MAXIMAL, OB>>,
        IntB: BorrowedStructure<IntegerCanonicalStructure>,
        SB: BorrowedStructure<FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>>,
    >
        Morphism<
            OrderIdealsStructure<K, KB, MAXIMAL, OB>,
            FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>,
        > for SubmoduleToIdeals<K, KB, MAXIMAL, OB, IB, IntB, SB>
    {
        fn domain(&self) -> &OrderIdealsStructure<K, KB, MAXIMAL, OB> {
            self.ideals()
        }

        fn range(&self) -> &FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB> {
            self.sublattices()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
        IB: BorrowedStructure<OrderIdealsStructure<K, KB, MAXIMAL, OB>>,
        IntB: BorrowedStructure<IntegerCanonicalStructure>,
        SB: BorrowedStructure<FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>>,
    >
        Function<
            OrderIdealsStructure<K, KB, MAXIMAL, OB>,
            FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>,
        > for SubmoduleToIdeals<K, KB, MAXIMAL, OB, IB, IntB, SB>
    {
        fn image(&self, x: &OrderIdeal) -> FinitelyFreeSubmodule<Integer> {
            #[cfg(debug_assertions)]
            self.ideals().is_element(x).unwrap();
            match x {
                OrderIdeal::Zero => self.sublattices().zero_submodule(),
                OrderIdeal::NonZero(sublattice) => sublattice.clone(),
            }
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
        IB: BorrowedStructure<OrderIdealsStructure<K, KB, MAXIMAL, OB>>,
        IntB: BorrowedStructure<IntegerCanonicalStructure>,
        SB: BorrowedStructure<FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>>,
    >
        InjectiveFunction<
            OrderIdealsStructure<K, KB, MAXIMAL, OB>,
            FinitelyFreeSubmoduleStructure<IntegerCanonicalStructure, IntB>,
        > for SubmoduleToIdeals<K, KB, MAXIMAL, OB, IB, IntB, SB>
    {
        fn try_preimage(&self, y: &FinitelyFreeSubmodule<Integer>) -> Option<OrderIdeal> {
            #[cfg(debug_assertions)]
            self.sublattices().is_element(y).unwrap();
            if y.rank() == 0 {
                Some(OrderIdeal::Zero)
            } else if self.ideals().does_sublattice_define_an_ideal(y) {
                let ideal = OrderIdeal::NonZero(y.clone());
                #[cfg(debug_assertions)]
                self.ideals().is_element(&ideal).unwrap();
                Some(ideal)
            } else {
                None
            }
        }
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    pub fn outbound_sublattices_inclusion<'a>(
        &'a self,
    ) -> sublattices_to_ideals::SubmoduleToIdeals<
        K,
        KB,
        MAXIMAL,
        OB,
        &'a Self,
        &'static IntegerCanonicalStructure,
        FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        >,
    > {
        sublattices_to_ideals::SubmoduleToIdeals::new(
            self,
            self.order().free_lattice_restructure().into_submodules(),
        )
    }

    pub fn into_outbound_sublattices_inclusion(
        self,
    ) -> sublattices_to_ideals::SubmoduleToIdeals<
        K,
        KB,
        MAXIMAL,
        OB,
        Self,
        &'static IntegerCanonicalStructure,
        FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        >,
    > {
        let sublattices = self.order().free_lattice_restructure().into_submodules();
        sublattices_to_ideals::SubmoduleToIdeals::new(self, sublattices)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> IdealsSignature<OrderWithBasis<K, KB, MAXIMAL>, OB> for OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    fn ring(&self) -> &OrderWithBasis<K, KB, MAXIMAL> {
        self.order()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    /// Construct an ideal from a Z-linear span
    pub fn from_integer_span(&self, span: Vec<Vec<Integer>>) -> OrderIdeal {
        for elem in &span {
            debug_assert!(self.order().is_element(elem).is_ok());
        }
        let n = self.order().n();
        let ideal = OrderIdeal::NonZero(
            Matrix::join_cols(
                n,
                span.into_iter()
                    .map(|elem| Matrix::from_cols(vec![elem]))
                    .collect(),
            )
            .col_span(),
        );
        #[cfg(debug_assertions)]
        self.is_element(&ideal).unwrap();
        ideal
    }

    /// The cardinality of the quotient ring, or 0 if the ideal is 0
    pub fn norm(&self, ideal: &OrderIdeal) -> Natural {
        debug_assert!(self.is_element(ideal).is_ok());
        match ideal {
            OrderIdeal::Zero => Natural::ZERO,
            OrderIdeal::NonZero(sublattice) => {
                let n = self.order().n();
                let cols = sublattice.basis();
                #[cfg(debug_assertions)]
                for col in &cols {
                    assert_eq!(col.len(), n);
                }
                let mat = Matrix::construct(n, n, |i, j| cols[i][j].clone());
                mat.det().unwrap().abs()
            }
        }
    }

    /// The sub z-module of points defining this ideal
    pub fn to_sublattice<'a>(&self, i: &'a OrderIdeal) -> Cow<'a, FinitelyFreeSubmodule<Integer>> {
        match i {
            OrderIdeal::Zero => Cow::Owned(
                self.order()
                    .free_lattice_restructure()
                    .submodules()
                    .zero_submodule(),
            ),
            OrderIdeal::NonZero(i) => Cow::Borrowed(i),
        }
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> IdealsArithmeticSignature<OrderWithBasis<K, KB, MAXIMAL>, OB>
    for OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    fn principal_ideal(&self, a: &Vec<Integer>) -> Self::Set {
        if self.order().is_zero(a) {
            Self::Set::Zero
        } else {
            let n = self.order().n();
            let ideal = self.from_integer_span(
                (0..n)
                    .map(|i| {
                        self.order()
                            .outbound_order_to_anf_inclusion()
                            .try_preimage(&self.order().anf().mul(
                                self.order().basis_vector(i),
                                &self.order().outbound_order_to_anf_inclusion().image(a),
                            ))
                            .unwrap()
                    })
                    .collect(),
            );
            debug_assert!(self.is_element(&ideal).is_ok());
            ideal
        }
    }

    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (OrderIdeal::Zero, OrderIdeal::Zero) => true,
            (OrderIdeal::NonZero(a_sublattice), OrderIdeal::NonZero(b_sublattice)) => self
                .ring()
                .free_lattice_restructure()
                .submodules()
                .equal(a_sublattice, b_sublattice),
            _ => false,
        }
    }

    fn contains_ideal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (_, OrderIdeal::Zero) => true,
            (OrderIdeal::Zero, OrderIdeal::NonZero { .. }) => {
                debug_assert_ne!(self.order().n(), 0);
                false
            }
            (OrderIdeal::NonZero(a_sublattice), OrderIdeal::NonZero(b_sublattice)) => self
                .ring()
                .free_lattice_restructure()
                .submodules()
                .contains(a_sublattice, b_sublattice),
        }
    }

    fn contains_element(&self, a: &Self::Set, x: &Vec<Integer>) -> bool {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.order().is_element(x).is_ok());
        match a {
            OrderIdeal::Zero => self.order().is_zero(x),
            OrderIdeal::NonZero(sublattice) => self
                .ring()
                .free_lattice_restructure()
                .submodules()
                .contains_element(sublattice, x),
        }
    }

    fn intersect(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (OrderIdeal::NonZero(a_sublattice), OrderIdeal::NonZero(b_sublattice)) => {
                Self::Set::NonZero(
                    self.order()
                        .free_lattice_restructure()
                        .submodules()
                        .intersect(a_sublattice.clone(), b_sublattice.clone()),
                )
            }
            _ => Self::Set::Zero,
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (OrderIdeal::Zero, OrderIdeal::Zero) => OrderIdeal::Zero,
            (OrderIdeal::Zero, OrderIdeal::NonZero { .. }) => b.clone(),
            (OrderIdeal::NonZero { .. }, OrderIdeal::Zero) => a.clone(),
            (OrderIdeal::NonZero(a_sublattice), OrderIdeal::NonZero(b_sublattice)) => {
                Self::Set::NonZero(
                    self.order()
                        .free_lattice_restructure()
                        .submodules()
                        .add(a_sublattice.clone(), b_sublattice.clone()),
                )
            }
        }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        match (a, b) {
            (OrderIdeal::NonZero(a_sublattice), OrderIdeal::NonZero(b_sublattice)) => {
                let n = self.order().n();
                let a_basis = a_sublattice.basis();
                let b_basis = b_sublattice.basis();
                debug_assert_eq!(a_basis.len(), n);
                debug_assert_eq!(b_basis.len(), n);

                let mut span = vec![];
                for i in 0..n {
                    for j in 0..n {
                        span.push(self.order().mul(&a_basis[i], &b_basis[j]));
                    }
                }
                self.from_integer_span(span)
            }
            _ => Self::Set::Zero,
        }
    }

    fn quotient(&self, i: &Self::Set, j: &Self::Set) -> Self::Set {
        let j = self.to_sublattice(j);
        let j = j.as_ref();
        self.quotient_ideal_by_sublattice(i, j)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
    OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
> OrderIdealsStructure<K, KB, MAXIMAL, OB>
{
    fn quotient_ideal_by_sublattice(
        &self,
        i: &<Self as SetSignature>::Set,
        s: &FinitelyFreeSubmodule<Integer>,
    ) -> <Self as SetSignature>::Set {
        let i = self.to_sublattice(i);
        let i = i.as_ref();

        let quotient_sublattice = self.order().quotient_sublattice(i, s);

        let quotient_ideal = if quotient_sublattice.rank() == 0 {
            OrderIdeal::Zero
        } else {
            OrderIdeal::NonZero(quotient_sublattice)
        };
        #[cfg(debug_assertions)]
        self.is_element(&quotient_ideal).unwrap();
        quotient_ideal
    }
}

// only maximal orders are Dedekind
impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<OrderWithBasis<K, KB, true>>,
> DedekindDomainIdealsSignature<OrderWithBasis<K, KB, true>, OB>
    for OrderIdealsStructure<K, KB, true, OB>
{
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<OrderWithBasis<K, KB, true>>,
> FactorableIdealsSignature<OrderWithBasis<K, KB, true>, OB>
    for OrderIdealsStructure<K, KB, true, OB>
where
    for<'a> RingOfIntegersIntegralExtensionWithIdeals<
        K,
        OrderWithBasis<K, KB, true>,
        &'a OrderWithBasis<K, KB, true>,
        RingOfIntegersToAlgebraicNumberFieldInclusion<
            K,
            OrderWithBasis<K, KB, true>,
            &'a OrderWithBasis<K, KB, true>,
        >,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        &'a OrderWithBasis<K, KB, true>,
        OrderIdealsStructure<K, KB, true, &'a OrderWithBasis<K, KB, true>>,
    >: IntegralClosureExtension<Z = IntegerCanonicalStructure, R = OrderWithBasis<K, KB, true>>
        + DedekindDomainExtension<
            IntegerCanonicalStructure,
            &'a OrderWithBasis<K, KB, true>,
            IdealsZ = IntegerIdealsStructure<IntegerCanonicalStructure>,
            IdealsR = OrderIdealsStructure<K, KB, true, &'a OrderWithBasis<K, KB, true>>,
        >,
{
    fn factor_ideal(
        &self,
        ideal: &Self::Set,
    ) -> Option<DedekindDomainIdealFactorization<Self::Set>> {
        Some(
            self.order()
                .outbound_roi_to_anf_inclusion()
                .zq_extension()
                .into_with_ideals()
                .factor_ideal(ideal)?
                .into_full_factorization(),
        )
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<OrderWithBasis<K, KB, true>>,
> OrderIdealsStructure<K, KB, true, OB>
where
    for<'a> RingOfIntegersIntegralExtensionWithIdeals<
        K,
        OrderWithBasis<K, KB, true>,
        &'a OrderWithBasis<K, KB, true>,
        RingOfIntegersToAlgebraicNumberFieldInclusion<
            K,
            OrderWithBasis<K, KB, true>,
            &'a OrderWithBasis<K, KB, true>,
        >,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        &'a OrderWithBasis<K, KB, true>,
        OrderIdealsStructure<K, KB, true, &'a OrderWithBasis<K, KB, true>>,
    >: IntegralClosureExtension<Z = IntegerCanonicalStructure, R = OrderWithBasis<K, KB, true>>
        + DedekindDomainExtension<
            IntegerCanonicalStructure,
            &'a OrderWithBasis<K, KB, true>,
            IdealsZ = IntegerIdealsStructure<IntegerCanonicalStructure>,
            IdealsR = OrderIdealsStructure<K, KB, true, &'a OrderWithBasis<K, KB, true>>,
        >,
{
    /// The order of the multiplicative group of the quotient modulo the ideal.
    pub fn euler_phi(&self, ideal: &OrderIdeal) -> Option<Natural> {
        Some(
            self.factorizations()
                .into_powers(self.factor_ideal(ideal)?)
                .iter()
                .map(|(prime_ideal, exponent)| {
                    let norm = self.norm(prime_ideal.ideal());
                    let e_minus_1 = exponent - Natural::ONE;
                    (&norm - Natural::ONE) * norm.pow(&e_minus_1)
                })
                .fold(Natural::ONE, |acc, x| acc * x),
        )
    }

    /// generate all ideals of norm equal to n
    pub fn all_ideals_norm_eq<'a>(
        &'a self,
        n: &Natural,
    ) -> Box<dyn 'a + Iterator<Item = OrderIdeal>> {
        match Integer::ideals().factor_ideal(n) {
            Some(n) => {
                let roi_to_anf =
                    RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(
                        self.order().clone(),
                    );
                let sq = roi_to_anf.zq_extension().into_with_ideals();
                Box::new(
                    Integer::structure()
                        .ideals()
                        .factorizations()
                        .into_powers(n)
                        .into_iter()
                        .map(|(p, k)| {
                            let k: usize = k.try_into().unwrap();
                            let primes_over_p = sq.factor_prime_ideal(p).into_factors();
                            num_partitions_part_pool(
                                k,
                                primes_over_p
                                    .iter()
                                    .map(|f| f.residue_class_degree)
                                    .collect(),
                            )
                            .map(|idxs| {
                                self.product(
                                    idxs.into_iter()
                                        .map(|i| primes_over_p[i].prime_ideal.ideal().clone())
                                        .collect(),
                                )
                            })
                            .collect::<Vec<OrderIdeal>>()
                        })
                        .multi_cartesian_product()
                        .map(|ideals| self.product(ideals)),
                )
            }
            None => Box::new(vec![self.zero_ideal()].into_iter()),
        }
    }

    /// generate all non-zero ideals of norm at most n
    pub fn all_nonzero_ideals_norm_le<'a>(
        &'a self,
        n: &'a Natural,
    ) -> Box<dyn 'a + Iterator<Item = OrderIdeal>> {
        Box::new(
            (1usize..)
                .map(Natural::from)
                .take_while(|m| m <= n)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }

    /// generate all ideals
    pub fn all_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = OrderIdeal>> {
        Box::new(
            (0usize..)
                .map(Natural::from)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }

    /// generate all non-zero ideals
    pub fn all_nonzero_ideals<'a>(&'a self) -> Box<dyn 'a + Iterator<Item = OrderIdeal>> {
        Box::new(
            (1usize..)
                .map(Natural::from)
                .flat_map(|m| self.all_ideals_norm_eq(&m)),
        )
    }

    /// given an ideal I and element a find an element b such that I = (a, b)
    pub fn ideal_other_generator(&self, g: &Vec<Integer>, ideal: &OrderIdeal) -> Vec<Integer> {
        debug_assert!(self.contains_element(ideal, g));
        debug_assert!(!self.order().is_zero(g));
        // prod_i p^{e_i}
        let ideal_factored = self.factor_ideal(ideal).unwrap();
        // prod_i p^{f_i} * prod_j q^{g_j}
        let g_factored = self.factor_ideal(&self.principal_ideal(g)).unwrap();
        // want b not in any q and in all p^{e_i} and not in any p^{e_i+1}

        // this is all b not in any q and in all p^{e_i}
        let b_set = self
            .ring()
            .free_lattice_restructure()
            .affine_subsets()
            .intersect_list(
                self.factorizations()
                    .to_powers(&ideal_factored)
                    .into_iter()
                    .map(|(p, k)| match self.nat_pow(p.ideal(), k) {
                        OrderIdeal::Zero => unreachable!(),
                        OrderIdeal::NonZero(pk_sublattice) => {
                            FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                                self.order()
                                    .free_lattice_restructure()
                                    .cosets()
                                    .from_submodule(pk_sublattice),
                            )
                        }
                    })
                    .chain(
                        self.factorizations()
                            .into_prime_support(g_factored)
                            .into_iter()
                            .filter(|prime_ideal| {
                                !self
                                    .factorizations()
                                    .to_prime_support(&ideal_factored)
                                    .into_iter()
                                    .any(|p| self.equal(p.ideal(), prime_ideal.ideal()))
                            })
                            .map(|q| match q.into_ideal() {
                                OrderIdeal::Zero => unreachable!(),
                                OrderIdeal::NonZero(q_sublattice) => {
                                    FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                                        self.order()
                                            .free_lattice_restructure()
                                            .cosets()
                                            .from_offset_and_submodule(
                                                &self.order().one(),
                                                q_sublattice,
                                            ),
                                    )
                                }
                            }),
                    )
                    .collect(),
            );

        //need to filter out the b in some p^{e_i+1}
        let rm_b_set = self
            .ring()
            .free_lattice_restructure()
            .affine_subsets()
            .intersect_list(
                self.factorizations()
                    .to_powers(&ideal_factored)
                    .into_iter()
                    .map(
                        |(p, k)| match self.nat_pow(p.ideal(), &(k + Natural::ONE)) {
                            OrderIdeal::Zero => unreachable!(),
                            OrderIdeal::NonZero(pk_sublattice) => {
                                FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                                    self.order()
                                        .free_lattice_restructure()
                                        .cosets()
                                        .from_submodule(pk_sublattice),
                                )
                            }
                        },
                    )
                    .collect(),
            );

        //if all basis elements of b_set were contain in rm_b_set then we'd have b_set contained in rm_b_set
        //but this is not the case, so some basis of b_set is not in rm_b_set

        self.order()
            .free_lattice_restructure()
            .affine_subsets()
            .affine_basis(&b_set)
            .into_iter()
            .find(|b| {
                !self
                    .ring()
                    .free_lattice_restructure()
                    .affine_subsets()
                    .contains_element(&rm_b_set, b)
            })
            .unwrap()
    }

    /// return two elements which generate the ideal
    pub fn ideal_two_generators(&self, ideal: &OrderIdeal) -> (Vec<Integer>, Vec<Integer>) {
        let (a, b) = match ideal {
            OrderIdeal::Zero => (self.order().zero(), self.order().zero()),
            OrderIdeal::NonZero(sublattice) => {
                let a = sublattice.basis().into_iter().next().unwrap();
                let b = self.ideal_other_generator(&a, ideal);
                (a, b)
            }
        };
        debug_assert!(self.equal(ideal, &self.generated_ideal(vec![a.clone(), b.clone()])));
        (a, b)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    OB: BorrowedStructure<OrderWithBasis<K, KB, true>>,
    OtoK: BorrowedMorphism<
            OrderWithBasis<K, KB, true>,
            K,
            RingOfIntegersToAlgebraicNumberFieldInclusion<K, OrderWithBasis<K, KB, true>, OB>,
        >,
    OIB: BorrowedStructure<OrderWithBasis<K, KB, true>>,
> DedekindDomainExtension<IntegerCanonicalStructure, OIB>
    for RingOfIntegersIntegralExtensionWithIdeals<
        K,
        OrderWithBasis<K, KB, true>,
        OB,
        OtoK,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        OIB,
        OrderIdealsStructure<K, KB, true, OIB>,
    >
{
    type IdealsZ = IntegerIdealsStructure<IntegerCanonicalStructure>;
    type IdealsR = OrderIdealsStructure<K, KB, true, OIB>;

    fn z_ideals(&self) -> &Self::IdealsZ {
        self.z_ideals()
    }

    fn r_ideals(&self) -> &Self::IdealsR {
        self.r_ideals()
    }

    fn ideal_norm(&self, ideal: &<Self::IdealsR as SetSignature>::Set) -> Natural {
        self.r_ideals().norm(ideal)
    }

    fn factor_prime_ideal(
        &self,
        prime_ideal: DedekindDomainPrimeIdeal<Natural>,
    ) -> DedekindExtensionIdealFactorsAbovePrime<Natural, <Self::IdealsR as SetSignature>::Set>
    {
        // https://en.wikipedia.org/wiki/Dedekind%E2%80%93Kummer_theorem
        let p = Integer::ideals().ideal_generator(prime_ideal.ideal());
        let anf = self.k_field();
        let roi = self.r_ring();
        let roi_ideals = self.r_ideals();
        let mod_p = Integer::structure().into_quotient_field_unchecked(p.clone());
        let poly_mod_p = mod_p.polynomials();
        let poly_roi = roi.polynomials();

        // alpha generates the algebraic number field but it is not necessarily an algebraic integer
        let alpha = anf.generator();
        // beta generates the algebraic number field and belongs to the ring of integers
        let beta = self.integral_scalar_multiple_r(&alpha);
        // factor the minimal polynomial of beta over the integers modulo p
        let beta_min_poly = self.min_poly_r_over_z(&beta);
        let beta_min_poly_factored = poly_mod_p.factor(&beta_min_poly).unwrap();
        // there is one prime ideal factor for each irreducible factor of beta's minimal polynomial modulo p
        // the prime ideal corresponding to an irreducible factor g(x) is generated by (p, g(beta))
        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
            prime_ideal,
            poly_mod_p
                .factorizations()
                .into_powers(beta_min_poly_factored)
                .into_iter()
                .map(|(g, power)| {
                    debug_assert!(g.is_monic());
                    let prime_ideal = roi_ideals.generated_ideal(vec![
                        self.z_to_r().image(&p),
                        poly_roi.evaluate(&g.apply_map(|c| self.z_to_r().image(c)), &beta),
                    ]);
                    // norm(I) = p^deg(g)
                    debug_assert_eq!(
                        roi_ideals.norm(&prime_ideal),
                        p.clone().abs().nat_pow(&g.degree().unwrap().into())
                    );
                    DedekindExtensionIdealFactorsAbovePrimeFactor {
                        prime_ideal: DedekindDomainPrimeIdeal::from_ideal_unchecked(prime_ideal),
                        residue_class_degree: g.degree().unwrap(),
                        power,
                    }
                })
                .collect(),
        )
    }

    fn factor_ideal(
        &self,
        ideal: &OrderIdeal,
    ) -> Option<DedekindExtensionIdealFactorization<Natural, OrderIdeal>> {
        let norm = self.ideal_norm(ideal);
        let norm_prime_factors = Integer::ideals().factor_ideal(&norm)?;
        Some(
            DedekindExtensionIdealFactorization::from_ideal_factors_above_primes(
                Integer::ideals()
                    .factorizations()
                    .into_prime_support(norm_prime_factors)
                    .into_iter()
                    .map(|prime| {
                        DedekindExtensionIdealFactorsAbovePrime::from_powers_unchecked(
                            prime.clone(),
                            self.factor_prime_ideal(prime.clone())
                                .into_factors()
                                .into_iter()
                                .filter_map(|factor_above_prime| {
                                    let k = self.r_ideals().largest_prime_ideal_factor_power(
                                        &factor_above_prime.prime_ideal,
                                        ideal,
                                    );
                                    if k == Natural::ZERO {
                                        None
                                    } else {
                                        Some(DedekindExtensionIdealFactorsAbovePrimeFactor {
                                            prime_ideal: factor_above_prime.prime_ideal,
                                            residue_class_degree: factor_above_prime
                                                .residue_class_degree,
                                            power: k,
                                        })
                                    }
                                })
                                .collect(),
                        )
                    })
                    .collect(),
            ),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{parsing::parse_rational_polynomial, polynomial::*};
    use algebraeon_nzq::*;

    #[test]
    fn ring_of_integers_ideals() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        let a = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE, Rational::ZERO]);
        let b = Polynomial::<Rational>::from_coeffs(vec![Rational::ZERO, Rational::ONE]);

        // Q[sqrt(2)]
        let anf = (x.pow(2) - 2)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = OrderWithBasis::new_maximal(anf.clone(), vec![a.clone(), b.clone()]).unwrap();
        let roi_ideals = roi.ideals();

        {
            // 1 + sqrt(2)
            let alpha = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&(&x + 1).into_verbose())
                .unwrap();

            // (a + b sqrt(2)) * (1 + sqrt(2)) = a(1 + sqrt(2)) + b(2 + sqrt(2))
            assert!(roi_ideals.equal(
                &roi_ideals.principal_ideal(&alpha),
                &roi_ideals.from_integer_span(vec![
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(1 + &x).into_verbose())
                        .unwrap(),
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(2 + &x).into_verbose())
                        .unwrap()
                ])
            ));
        }

        {
            // 6
            let alpha = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&(6 * x.pow(0)).into_verbose())
                .unwrap();
            // 15
            let beta = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&(15 * x.pow(0)).into_verbose())
                .unwrap();

            let alpha_ideal = roi_ideals.principal_ideal(&alpha);
            let beta_ideal = roi_ideals.principal_ideal(&beta);

            let alpha_beta_add = roi_ideals.add(&alpha_ideal, &beta_ideal);
            let alpha_beta_intersect = roi_ideals.intersect(&alpha_ideal, &beta_ideal);
            let alpha_beta_mul = roi_ideals.mul(&alpha_ideal, &beta_ideal);

            // sum is 3
            assert!(roi_ideals.equal(
                &alpha_beta_add,
                &roi_ideals.from_integer_span(vec![
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(3 * x.pow(0)).into_verbose())
                        .unwrap(),
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(3 * x.pow(1)).into_verbose())
                        .unwrap()
                ])
            ));

            // intersection is 30
            assert!(roi_ideals.equal(
                &alpha_beta_intersect,
                &roi_ideals.from_integer_span(vec![
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(30 * x.pow(0)).into_verbose())
                        .unwrap(),
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(30 * x.pow(1)).into_verbose())
                        .unwrap()
                ])
            ));

            // product is 90
            assert!(roi_ideals.equal(
                &alpha_beta_mul,
                &roi_ideals.from_integer_span(vec![
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(90 * x.pow(0)).into_verbose())
                        .unwrap(),
                    roi.outbound_order_to_anf_inclusion()
                        .try_preimage(&(90 * x.pow(1)).into_verbose())
                        .unwrap()
                ])
            ));
        }
    }

    #[test]
    fn test_sublattices_to_ideals() {
        let anf = parse_rational_polynomial("x^2+3", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();
        let order = anf
            .order(vec![
                parse_rational_polynomial("1", "x").unwrap(),
                parse_rational_polynomial("x", "x").unwrap(),
            ])
            .unwrap();
        let ideals = order.ideals();
        let sublattices = order.free_lattice_restructure().into_submodules();

        // f : {ideals of the order} -> {z-sublattices of the order}
        let f = ideals.outbound_sublattices_inclusion();

        assert!(
            f.try_preimage(&sublattices.matrix_row_span(Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1)],
            ])))
            .is_some()
        );

        assert!(f.try_preimage(&sublattices.zero_submodule()).is_some());

        assert!(
            f.try_preimage(
                &sublattices
                    .matrix_row_span(Matrix::from_row(vec![Integer::from(1), Integer::from(0)]))
            )
            .is_none()
        );

        // S = {2*a + b*sqrt(3) : a,b in ZZ} is not an ideal because sqrt(3) in R times sqrt(3) in S is 3 which is not in S
        assert!(
            f.try_preimage(&sublattices.matrix_row_span(Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1)],
            ])))
            .is_none()
        );
    }

    #[test]
    fn test_count_all_ideals_norm_eq() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.into_ring_of_integers();
        let roi_ideals = roi.ideals();

        assert_eq!(
            roi_ideals
                .all_ideals_norm_eq(&Natural::from(5040 as u32))
                .collect::<Vec<_>>()
                .len(),
            0
        );
        assert_eq!(
            roi_ideals
                .all_ideals_norm_eq(&Natural::from(5040 * 7 as u32))
                .collect::<Vec<_>>()
                .len(),
            2
        );
    }

    #[test]
    fn test_euler_phi_of_principal_ideal() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Construct the number field Q(i), which has ring of integers Z[i]
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.ring_of_integers();
        let roi_ideals = roi.ideals();

        // Consider the ideal (5)
        let ideal = roi_ideals.principal_ideal(&roi.from_int(5));

        let phi = roi_ideals.euler_phi(&ideal).unwrap();
        assert_eq!(phi, Natural::from(16u32));
    }

    #[test]
    fn test_is_square_ideal() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Work in the ring of integers Z[i]
        let anf = (x.pow(2) + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi = anf.ring_of_integers();
        let roi_ideals = roi.ideals();
        let one_plus_i = roi
            .outbound_order_to_anf_inclusion()
            .try_preimage(&(&x + 1).into_verbose())
            .unwrap();
        let gaussian_prime = roi_ideals.principal_ideal(&one_plus_i);

        let zero_sqrt = roi_ideals.sqrt_if_square(&roi_ideals.zero_ideal()).unwrap();
        assert!(roi_ideals.equal(&zero_sqrt, &roi_ideals.zero_ideal()));
        assert!(roi_ideals.is_square(&roi_ideals.zero_ideal()));
        assert!(!roi_ideals.is_squarefree(&roi_ideals.zero_ideal()));

        // (2) = (1 + i)^2 in Z[i]
        let two_ideal = roi_ideals.principal_ideal(&roi.from_int(2));
        assert!(roi_ideals.is_square(&two_ideal));
        let sqrt_two = roi_ideals.sqrt_if_square(&two_ideal).unwrap();
        assert!(roi_ideals.equal(&sqrt_two, &gaussian_prime));
        assert!(!roi_ideals.is_squarefree(&two_ideal));

        // (3) stays prime in Z[i] so it cannot be a square.
        let three_ideal = roi_ideals.principal_ideal(&roi.from_int(3));
        assert!(!roi_ideals.is_square(&three_ideal));
        assert!(roi_ideals.sqrt_if_square(&three_ideal).is_none());
        assert!(roi_ideals.is_squarefree(&three_ideal));

        // Squares built from a non-trivial prime ideal should be detected as well.
        assert!(!roi_ideals.is_square(&gaussian_prime));
        assert!(roi_ideals.is_squarefree(&gaussian_prime));
        let gaussian_prime_square = roi_ideals.mul(&gaussian_prime, &gaussian_prime);
        assert!(roi_ideals.is_square(&gaussian_prime_square));
        let sqrt_gaussian_square = roi_ideals.sqrt_if_square(&gaussian_prime_square).unwrap();
        assert!(roi_ideals.equal(&sqrt_gaussian_square, &gaussian_prime));
        assert!(!roi_ideals.is_squarefree(&gaussian_prime_square));
    }
}
