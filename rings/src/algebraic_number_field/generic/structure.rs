use crate::{
    algebraic_number_field::{OrderWithBasis, RingOfIntegersIntegralExtension},
    linear::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    structure::{
        AdditiveGroupSignature, CharZeroFieldSignature, CharZeroRingSignature,
        DedekindDomainSignature, FieldOfFractionsInclusion, FiniteDimensionalFieldExtension,
        FiniteRankFreeRingExtension, MetaGreatestCommonDivisorSignature, RingHomomorphism,
    },
};
use algebraeon_nzq::{
    Integer, IntegerCanonicalStructure, Rational, RationalCanonicalStructure, traits::Fraction,
};
use algebraeon_sets::structure::{
    BorrowedStructure, FiniteSetSignature, Function, InjectiveFunction, Morphism, SetSignature,
};
use std::marker::PhantomData;

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn inbound_finite_dimensional_rational_extension(&self) -> Self::RationalInclusion<&Self>;
    fn into_inbound_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

    /// The dimension of this algebraic number field as a vector space over the rational numbers
    fn n(&self) -> usize {
        self.inbound_finite_dimensional_rational_extension()
            .degree()
    }

    /// An element which generates this algebraic number field when adjoined to the rational numbers
    /// Such an element always exists by the primitive element theorem
    fn generator(&self) -> Self::Set;

    /// Determine whether an element is integral over the integers i.e. is it a root of a monic integer polynomial
    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    /// The discriminant of this algebraic number field i.e. the discriminant of its ring of integers
    /// Implementations should not compute this by constructing the ring of integers, as the constructor for a maximal OrderWithBasis calls this function to validate its input
    fn discriminant(&self) -> Integer;

    /// A list of self.n() elements which generate the ring of integers as a Z-module
    fn integral_basis(&self) -> Vec<Self::Set>;

    fn ring_of_integers(&self) -> OrderWithBasis<Self, &Self, true> {
        OrderWithBasis::new_maximal_unchecked(self, self.integral_basis())
    }
    fn into_ring_of_integers(self) -> OrderWithBasis<Self, Self, true> {
        let basis = self.integral_basis();
        OrderWithBasis::new_maximal_unchecked(self, basis)
    }

    fn order(&self, basis: Vec<Self::Set>) -> Result<OrderWithBasis<Self, &Self, false>, String> {
        OrderWithBasis::new(self, basis)
    }
    fn into_order(
        self,
        basis: Vec<Self::Set>,
    ) -> Result<OrderWithBasis<Self, Self, false>, String> {
        OrderWithBasis::new(self, basis)
    }

    /// The LCM of the denominators of the coefficients of the minimal polynomial of a.
    ///
    /// It may well be >1 even when the element a is an algebraic integer.
    fn min_poly_denominator_lcm(&self, a: &Self::Set) -> Integer {
        Integer::lcm_list(
            self.inbound_finite_dimensional_rational_extension()
                .min_poly(a)
                .coeffs()
                .map(|c| Integer::from(c.denominator()))
                .collect(),
        )
    }

    /// A scalar multiple of $a$ which is an algebraic integer.
    ///
    /// It need not return $a$ itself when $a$ is already an algebraic integer.
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingSignature<K: AlgebraicNumberFieldSignature>:
    DedekindDomainSignature + CharZeroRingSignature
{
    fn n(&self) -> usize {
        self.anf().n()
    }

    fn anf(&self) -> &K;

    /// A list of self.n() elements which generate this ring as a Z-module
    fn integral_basis(&self) -> Vec<Self::Set>;

    fn to_anf(&self, x: &Self::Set) -> K::Set;

    fn try_from_anf(&self, y: &K::Set) -> Option<Self::Set>;

    fn order(&self) -> OrderWithBasis<K, &K, true> {
        OrderWithBasis::new_maximal_unchecked(
            self.anf(),
            self.integral_basis()
                .into_iter()
                .map(|v| self.outbound_roi_to_anf_inclusion().image(&v))
                .collect(),
        )
    }

    fn into_outbound_roi_to_anf_inclusion(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }

    fn outbound_roi_to_anf_inclusion(
        &self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, &Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }

    fn inbound_order_inclusion<
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    >(
        &self,
        order: OB,
    ) -> order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion<
        K,
        Self,
        &Self,
        KOB,
        MAXIMAL,
        OB,
    > {
        order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion::new(self, order)
    }

    fn inbound_order_isomorphism<
        KOB: BorrowedStructure<K>,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, true>>,
    >(
        &self,
        order: OB,
    ) -> order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion<
        K,
        Self,
        &Self,
        KOB,
        true,
        OB,
    > {
        order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion::new(self, order)
    }
}

mod ring_of_integers_to_algebraic_number_field_inclusion {
    use crate::structure::IntegralClosureExtension;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > {
        _roi: PhantomData<R>,
        roi: RB,
        _anf: PhantomData<K>,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        pub fn from_ring_of_integers(roi: RB) -> Self {
            Self {
                _roi: PhantomData,
                _anf: PhantomData,
                roi,
            }
        }

        pub fn roi(&self) -> &R {
            self.roi.borrow()
        }

        pub fn anf(&self) -> &K {
            self.roi().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > Morphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        fn domain(&self) -> &R {
            self.roi()
        }

        fn range(&self) -> &K {
            self.anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > Function<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        fn image(&self, x: &<R as SetSignature>::Set) -> <K as SetSignature>::Set {
            self.roi().to_anf(x)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > InjectiveFunction<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<<R as SetSignature>::Set> {
            self.roi().try_from_anf(y)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > RingHomomorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > FieldOfFractionsInclusion<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        fn numerator_and_denominator(&self, a: &<K>::Set) -> (<R>::Set, <R>::Set) {
            self.zq_extension()
                .r_to_k_field_of_fractions()
                .numerator_and_denominator(a)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
    > RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
    {
        pub fn zq_extension(
            &self,
        ) -> RingOfIntegersIntegralExtension<
            K,
            R,
            &R,
            RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, &R>,
        > {
            RingOfIntegersIntegralExtension::new_integer_extension(
                RingOfIntegersToAlgebraicNumberFieldInclusion {
                    _roi: PhantomData,
                    _anf: PhantomData,
                    roi: self.domain(),
                },
            )
        }
    }
}
pub(crate) use ring_of_integers_to_algebraic_number_field_inclusion::RingOfIntegersToAlgebraicNumberFieldInclusion;

mod order_to_ring_of_integers_inclusion {
    use algebraeon_sets::structure::BijectiveFunction;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct OrderToRingOfIntegersInclusion<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > {
        _k: PhantomData<K>,
        _r: PhantomData<R>,
        roi: RB,
        _kob: PhantomData<KOB>,
        order: OB,
        order_basis_in_roi: Vec<R::Set>,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > OrderToRingOfIntegersInclusion<K, R, RB, KOB, MAXIMAL, OB>
    {
        pub fn new(roi: RB, order: OB) -> Self {
            let order_basis_in_roi = order
                .borrow()
                .basis()
                .iter()
                .map(|bv| {
                    roi.borrow()
                        .outbound_roi_to_anf_inclusion()
                        .try_preimage(bv)
                        .unwrap()
                })
                .collect();
            Self {
                _k: PhantomData,
                _r: PhantomData,
                roi,
                _kob: PhantomData,
                order,
                order_basis_in_roi,
            }
        }

        pub fn n(&self) -> usize {
            self.anf().n()
        }

        pub fn anf(&self) -> &K {
            debug_assert_eq!(self.roi().anf(), self.order().anf());
            self.order().anf()
        }

        pub fn roi(&self) -> &R {
            self.roi.borrow()
        }

        pub fn order(&self) -> &OrderWithBasis<K, KOB, MAXIMAL> {
            self.order.borrow()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > Morphism<OrderWithBasis<K, KOB, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, RB, KOB, MAXIMAL, OB>
    {
        fn domain(&self) -> &OrderWithBasis<K, KOB, MAXIMAL> {
            self.order()
        }

        fn range(&self) -> &R {
            self.roi()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > Function<OrderWithBasis<K, KOB, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, RB, KOB, MAXIMAL, OB>
    {
        fn image(&self, x: &Vec<Integer>) -> <R as SetSignature>::Set {
            self.roi().sum(
                (0..self.n())
                    .map(|i| {
                        self.roi()
                            .mul(&self.roi().from_int(&x[i]), &self.order_basis_in_roi[i])
                    })
                    .collect(),
            )
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > RingHomomorphism<OrderWithBasis<K, KOB, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, RB, KOB, MAXIMAL, OB>
    {
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > InjectiveFunction<OrderWithBasis<K, KOB, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, RB, KOB, MAXIMAL, OB>
    {
        fn try_preimage(&self, y: &<R as SetSignature>::Set) -> Option<Vec<Integer>> {
            self.order()
                .outbound_order_to_anf_inclusion()
                .try_preimage(&self.roi().outbound_roi_to_anf_inclusion().image(y))
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        RB: BorrowedStructure<R>,
        KOB: BorrowedStructure<K>,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, true>>,
    > BijectiveFunction<OrderWithBasis<K, KOB, true>, R>
        for OrderToRingOfIntegersInclusion<K, R, RB, KOB, true, OB>
    {
    }
}

mod anf_inclusion {
    use super::*;
    use crate::{
        matrix::Matrix,
        structure::{
            FieldOfFractionsInclusion, IntegralClosureExtension, MetaCharZeroRingSignature,
        },
    };

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
    > {
        _k: PhantomData<K>,
        _lat: PhantomData<IntegerModule>,
        integer_submodule: IntegerModuleB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
    >
        AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
        >
    {
        pub fn new(integer_submodule: IntegerModuleB) -> Self {
            Self {
                _k: PhantomData,
                _lat: PhantomData,
                integer_submodule,
            }
        }

        pub fn integer_submodule(&self) -> &IntegerModule {
            self.integer_submodule.borrow()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
    > Morphism<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
        >
    {
        fn domain(&self) -> &IntegerModule {
            self.integer_submodule()
        }

        fn range(&self) -> &K {
            self.integer_submodule().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
    > Function<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
        >
    {
        fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
            debug_assert!(self.integer_submodule().validate_element(x).is_ok());
            let k = self.integer_submodule().anf();
            let n = k.n();
            debug_assert_eq!(n, x.len());
            k.sum(
                (0..n)
                    .map(|i| k.mul(&k.from_int(&x[i]), self.integer_submodule().basis_vector(i)))
                    .collect(),
            )
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
    > InjectiveFunction<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
        >
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
            let k = self.integer_submodule().anf();
            let n = k.n();
            debug_assert!(k.validate_element(y).is_ok());
            let mat = Matrix::join_cols(
                n,
                (0..n)
                    .map(|i| {
                        k.inbound_finite_dimensional_rational_extension()
                            .to_col(self.integer_submodule().basis_vector(i))
                    })
                    .collect(),
            );
            let y = k.inbound_finite_dimensional_rational_extension().to_vec(y);
            let x_rat = mat.col_solve(&y)?;
            let mut x_int = Vec::with_capacity(n);
            for c_rat in x_rat {
                x_int.push(c_rat.try_to_int()?);
            }
            Some(x_int)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > RingHomomorphism<OrderWithBasis<K, KOB, MAXIMAL>, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, KOB, MAXIMAL>,
            OB,
        >
    {
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    > FieldOfFractionsInclusion<OrderWithBasis<K, KOB, MAXIMAL>, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, KOB, MAXIMAL>,
            OB,
        >
    {
        fn numerator_and_denominator(
            &self,
            a: &<K>::Set,
        ) -> (
            <OrderWithBasis<K, KOB, MAXIMAL> as SetSignature>::Set,
            <OrderWithBasis<K, KOB, MAXIMAL> as SetSignature>::Set,
        ) {
            self.zq_extension()
                .r_to_k_field_of_fractions()
                .numerator_and_denominator(a)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KOB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
    >
        AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, KOB, MAXIMAL>,
            OB,
        >
    {
        pub fn zq_extension(
            &self,
        ) -> order_integral_extension::OrderIntegralExtension<
            K,
            KOB,
            MAXIMAL,
            OB,
            &AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                K,
                OrderWithBasis<K, KOB, MAXIMAL>,
                OB,
            >,
        > {
            order_integral_extension::OrderIntegralExtension::new_integer_extension(self)
        }
    }

    mod order_integral_extension {
        use std::borrow::Cow;

        use crate::{
            algebraic_number_field::OrderIdealsStructure,
            integer::ideal::IntegerIdealsStructure,
            structure::{
                IdealsSignature, IntegralClosureExtension, PrincipalIntegerMap, RingSignature,
                RingToIdealsSignature,
            },
        };

        use super::*;
        use algebraeon_sets::structure::{BorrowedMorphism, MetaType};

        /// Q -> K
        /// ↑    ↑
        /// Z -> R
        ///
        /// Where Q is the rationals, Z is the integers, K is an algebraic number field, R is its ring of integers
        ///
        #[derive(Debug, Clone)]
        pub struct OrderIntegralExtension<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
        > {
            _k: PhantomData<K>,
            _kob: PhantomData<KOB>,
            _rb: PhantomData<RB>,
            r_to_k: RtoK,
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
        > OrderIntegralExtension<K, KOB, MAXIMAL, RB, RtoK>
        {
            pub fn new_integer_extension(r_to_k: RtoK) -> Self {
                Self {
                    _k: PhantomData,
                    _kob: PhantomData,
                    _rb: PhantomData,
                    r_to_k,
                }
            }

            pub fn with_ideals(
                &self,
            ) -> OrderIntegralExtensionWithIdeals<
                K,
                KOB,
                MAXIMAL,
                RB,
                &AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                    K,
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    RB,
                >,
                IntegerIdealsStructure<IntegerCanonicalStructure>,
                &OrderWithBasis<K, KOB, MAXIMAL>,
                OrderIdealsStructure<K, KOB, MAXIMAL, &OrderWithBasis<K, KOB, MAXIMAL>>,
            > {
                let ideals_z = Integer::structure().into_ideals();
                let ideals_r = self.r_to_k.borrow().domain().ideals();
                OrderIntegralExtensionWithIdeals::new_integer_extension(
                    self.r_to_k.borrow(),
                    ideals_z,
                    ideals_r,
                )
            }

            pub fn into_with_ideals(
                self,
            ) -> OrderIntegralExtensionWithIdeals<
                K,
                KOB,
                MAXIMAL,
                RB,
                RtoK,
                IntegerIdealsStructure<IntegerCanonicalStructure>,
                OrderWithBasis<K, KOB, MAXIMAL>,
                OrderIdealsStructure<K, KOB, MAXIMAL, OrderWithBasis<K, KOB, MAXIMAL>>,
            > {
                let ideals_z = Integer::structure().into_ideals();
                let ideals_r = self.r_to_k().domain().clone().into_ideals();
                OrderIntegralExtensionWithIdeals::new_integer_extension(
                    self.r_to_k,
                    ideals_z,
                    ideals_r,
                )
            }
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
        > IntegralClosureExtension for OrderIntegralExtension<K, KOB, MAXIMAL, RB, RtoK>
        {
            type QKBasis = K::Basis;
            type Z = IntegerCanonicalStructure;
            type Q = RationalCanonicalStructure;
            type R = OrderWithBasis<K, KOB, MAXIMAL>;
            type K = K;
            type ZQ<BZ: BorrowedStructure<Self::Z>, BQ: BorrowedStructure<Self::Q>> =
                PrincipalIntegerMap<RationalCanonicalStructure, RationalCanonicalStructure>;
            type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>> =
                PrincipalIntegerMap<OrderWithBasis<K, KOB, MAXIMAL>, BR>;
            type QK<BQ: BorrowedStructure<Self::Q>, BK: BorrowedStructure<Self::K>> =
                K::RationalInclusion<BK>;
            type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>> =
                AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                    K,
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    RB,
                >;

            fn z_ring(&self) -> &Self::Z {
                Integer::structure_ref()
            }
            fn r_ring(&self) -> &Self::R {
                self.r_to_k.borrow().domain()
            }
            fn q_field(&self) -> &Self::Q {
                Rational::structure_ref()
            }
            fn k_field(&self) -> &Self::K {
                self.r_to_k.borrow().range()
            }

            fn z_to_q<'a>(&'a self) -> Cow<'a, Self::ZQ<&'a Self::Z, &'a Self::Q>> {
                Cow::Owned(Rational::structure().into_inbound_principal_integer_map())
            }
            fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>> {
                Cow::Owned(self.r_ring().inbound_principal_integer_map())
            }
            fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>> {
                Cow::Owned(
                    self.k_field()
                        .inbound_finite_dimensional_rational_extension(),
                )
            }
            fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>> {
                Cow::Borrowed(self.r_to_k.borrow())
            }

            fn integralize_multiplier(&self, alpha: &<Self::K as SetSignature>::Set) -> Integer {
                if self.k_field().is_algebraic_integer(alpha) {
                    Integer::ONE
                } else {
                    self.k_field().min_poly_denominator_lcm(alpha)
                }
            }
        }

        #[derive(Debug, Clone)]
        pub struct OrderIntegralExtensionWithIdeals<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
            RIB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            IdealsR: IdealsSignature<OrderWithBasis<K, KOB, MAXIMAL>, RIB>,
        > {
            _k: PhantomData<K>,
            _kob: PhantomData<KOB>,
            _rb: PhantomData<RB>,
            r_to_k: RtoK,
            ideals_z: IdealsZ,
            _rib: PhantomData<RIB>,
            ideals_r: IdealsR,
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
            RIB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            IdealsR: IdealsSignature<OrderWithBasis<K, KOB, MAXIMAL>, RIB>,
        > OrderIntegralExtensionWithIdeals<K, KOB, MAXIMAL, RB, RtoK, IdealsZ, RIB, IdealsR>
        {
            pub fn new_integer_extension(
                r_to_k: RtoK,
                ideals_z: IdealsZ,
                ideals_r: IdealsR,
            ) -> Self {
                Self {
                    _k: PhantomData,
                    _kob: PhantomData,
                    _rb: PhantomData,
                    r_to_k,
                    ideals_z,
                    _rib: PhantomData,
                    ideals_r,
                }
            }

            pub fn z_ideals(&self) -> &IdealsZ {
                &self.ideals_z
            }

            pub fn r_ideals(&self) -> &IdealsR {
                &self.ideals_r
            }
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            KOB: BorrowedStructure<K>,
            const MAXIMAL: bool,
            RB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            RtoK: BorrowedMorphism<
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    K,
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, KOB, MAXIMAL>,
                        RB,
                    >,
                >,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure, IntegerCanonicalStructure>,
            RIB: BorrowedStructure<OrderWithBasis<K, KOB, MAXIMAL>>,
            IdealsR: IdealsSignature<OrderWithBasis<K, KOB, MAXIMAL>, RIB>,
        > IntegralClosureExtension
            for OrderIntegralExtensionWithIdeals<K, KOB, MAXIMAL, RB, RtoK, IdealsZ, RIB, IdealsR>
        {
            type QKBasis = K::Basis;
            type Z = IntegerCanonicalStructure;
            type Q = RationalCanonicalStructure;
            type R = OrderWithBasis<K, KOB, MAXIMAL>;
            type K = K;
            type ZQ<BZ: BorrowedStructure<Self::Z>, BQ: BorrowedStructure<Self::Q>> =
                PrincipalIntegerMap<RationalCanonicalStructure, RationalCanonicalStructure>;
            type ZR<BZ: BorrowedStructure<Self::Z>, BR: BorrowedStructure<Self::R>> =
                PrincipalIntegerMap<OrderWithBasis<K, KOB, MAXIMAL>, BR>;
            type QK<BQ: BorrowedStructure<Self::Q>, BK: BorrowedStructure<Self::K>> =
                K::RationalInclusion<BK>;
            type RK<BR: BorrowedStructure<Self::R>, BK: BorrowedStructure<Self::K>> =
                AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                    K,
                    OrderWithBasis<K, KOB, MAXIMAL>,
                    RB,
                >;

            fn z_ring(&self) -> &Self::Z {
                Integer::structure_ref()
            }
            fn r_ring(&self) -> &Self::R {
                self.r_to_k.borrow().domain()
            }
            fn q_field(&self) -> &Self::Q {
                Rational::structure_ref()
            }
            fn k_field(&self) -> &Self::K {
                self.r_to_k.borrow().range()
            }

            fn z_to_q<'a>(&'a self) -> Cow<'a, Self::ZQ<&'a Self::Z, &'a Self::Q>> {
                Cow::Owned(Rational::structure().into_inbound_principal_integer_map())
            }
            fn z_to_r<'a>(&'a self) -> Cow<'a, Self::ZR<&'a Self::Z, &'a Self::R>> {
                Cow::Owned(self.r_ring().inbound_principal_integer_map())
            }
            fn q_to_k<'a>(&'a self) -> Cow<'a, Self::QK<&'a Self::Q, &'a Self::K>> {
                Cow::Owned(
                    self.k_field()
                        .inbound_finite_dimensional_rational_extension(),
                )
            }
            fn r_to_k<'a>(&'a self) -> Cow<'a, Self::RK<&'a Self::R, &'a Self::K>> {
                Cow::Borrowed(self.r_to_k.borrow())
            }

            fn integralize_multiplier(&self, alpha: &<Self::K as SetSignature>::Set) -> Integer {
                if self.k_field().is_algebraic_integer(alpha) {
                    Integer::ONE
                } else {
                    self.k_field().min_poly_denominator_lcm(alpha)
                }
            }
        }
    }
}

mod integer_submodule_inclusion {
    use super::*;
    use crate::{
        linear::finitely_free_submodule::{FinitelyFreeSubmodule, FinitelyFreeSubmoduleStructure},
        structure::RingSignature,
    };
    use algebraeon_sets::structure::{BorrowedMorphism, BorrowedStructure};
    use std::marker::PhantomData;

    #[derive(Debug, Clone)]
    pub struct IntegerSubmoduleInclusion<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    > {
        _k: PhantomData<K>,
        _lat: PhantomData<IntegerModule>,
        _sublat: PhantomData<IntegerSubmodule>,
        integer_module: IntegerModuleB,
        integer_submodule: IntegerSubmoduleB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    >
        IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >
    {
        pub fn integer_module(&self) -> &IntegerModule {
            self.integer_module.borrow()
        }

        pub fn integer_submodule(&self) -> &IntegerSubmodule {
            self.integer_submodule.borrow()
        }

        fn check(&self) -> Result<(), String> {
            if self
                .integer_module()
                .contains_integer_submodule(self.integer_submodule())
            {
                Ok(())
            } else {
                Err("Integer module does not contain submodule".to_string())
            }
        }

        fn new_impl(integer_module: IntegerModuleB, integer_submodule: IntegerSubmoduleB) -> Self {
            Self {
                _k: PhantomData,
                _lat: PhantomData,
                _sublat: PhantomData,
                integer_module,
                integer_submodule,
            }
        }

        pub fn new(
            integer_module: IntegerModuleB,
            integer_submodule: IntegerSubmoduleB,
        ) -> Option<Self> {
            let s = Self::new_impl(integer_module, integer_submodule);
            if s.check().is_ok() { Some(s) } else { None }
        }

        pub fn new_unchecked(
            integer_module: IntegerModuleB,
            integer_submodule: IntegerSubmoduleB,
        ) -> Self {
            let s = Self::new_impl(integer_module, integer_submodule);
            #[cfg(debug_assertions)]
            s.check().unwrap();
            s
        }

        pub fn integer_submodules_inclusion(
            &self,
        ) -> IntegerSubmoduleIntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
            &Self,
        > {
            IntegerSubmoduleIntegerSubmoduleInclusion::new(self)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    > Morphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >
    {
        fn domain(&self) -> &IntegerSubmodule {
            self.integer_submodule()
        }

        fn range(&self) -> &IntegerModule {
            self.integer_module()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K> + RingSignature,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K> + RingSignature,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    > RingHomomorphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >
    {
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    > Function<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >
    {
        fn image(
            &self,
            x: &<IntegerSubmodule as SetSignature>::Set,
        ) -> <IntegerModule as SetSignature>::Set {
            self.integer_module()
                .outbound_order_to_anf_inclusion()
                .try_preimage(
                    &self
                        .integer_submodule()
                        .outbound_order_to_anf_inclusion()
                        .image(x),
                )
                .unwrap()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    > InjectiveFunction<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >
    {
        fn try_preimage(
            &self,
            y: &<IntegerModule as SetSignature>::Set,
        ) -> Option<<IntegerSubmodule as SetSignature>::Set> {
            self.integer_submodule()
                .outbound_order_to_anf_inclusion()
                .try_preimage(
                    &self
                        .integer_module()
                        .outbound_order_to_anf_inclusion()
                        .image(y),
                )
        }
    }

    #[derive(Debug, Clone)]
    pub struct IntegerSubmoduleIntegerSubmoduleInclusion<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
        IntegerSubmoduleInclusionB: BorrowedMorphism<
                IntegerSubmodule,
                IntegerModule,
                IntegerSubmoduleInclusion<
                    K,
                    IntegerModule,
                    IntegerModuleB,
                    IntegerSubmodule,
                    IntegerSubmoduleB,
                >,
            >,
    > {
        _k: PhantomData<K>,
        _lat: PhantomData<IntegerModule>,
        _latb: PhantomData<IntegerModuleB>,
        _sublat: PhantomData<IntegerSubmodule>,
        _sublatb: PhantomData<IntegerSubmoduleB>,
        integer_submodule_to_module: IntegerSubmoduleInclusionB,
        integer_module_to_submodule: FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        >,
        integer_submodule_to_submodules: FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        >,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
        IntegerSubmoduleInclusionB: BorrowedMorphism<
                IntegerSubmodule,
                IntegerModule,
                IntegerSubmoduleInclusion<
                    K,
                    IntegerModule,
                    IntegerModuleB,
                    IntegerSubmodule,
                    IntegerSubmoduleB,
                >,
            >,
    >
        IntegerSubmoduleIntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
            IntegerSubmoduleInclusionB,
        >
    {
        pub fn new(integer_submodule_to_module: IntegerSubmoduleInclusionB) -> Self {
            Self {
                _k: PhantomData,
                _lat: PhantomData,
                _latb: PhantomData,
                _sublat: PhantomData,
                _sublatb: PhantomData,
                integer_module_to_submodule: integer_submodule_to_module
                    .borrow()
                    .integer_module()
                    .free_integer_submodule_restructure()
                    .into_submodules(),
                integer_submodule_to_submodules: integer_submodule_to_module
                    .borrow()
                    .integer_submodule()
                    .free_integer_submodule_restructure()
                    .into_submodules(),
                integer_submodule_to_module,
            }
        }

        pub fn integer_submodule_to_module(
            &self,
        ) -> &IntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
        > {
            self.integer_submodule_to_module.borrow()
        }

        pub fn integer_module_to_submodule(
            &self,
        ) -> &FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        > {
            &self.integer_module_to_submodule
        }

        pub fn integer_submodule_to_submodules(
            &self,
        ) -> &FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        > {
            &self.integer_submodule_to_submodules
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
        IntegerSubmoduleInclusionB: BorrowedMorphism<
                IntegerSubmodule,
                IntegerModule,
                IntegerSubmoduleInclusion<
                    K,
                    IntegerModule,
                    IntegerModuleB,
                    IntegerSubmodule,
                    IntegerSubmoduleB,
                >,
            >,
    >
        Morphism<
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
        >
        for IntegerSubmoduleIntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
            IntegerSubmoduleInclusionB,
        >
    {
        fn domain(
            &self,
        ) -> &FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        > {
            &self.integer_submodule_to_submodules
        }

        fn range(
            &self,
        ) -> &FinitelyFreeSubmoduleStructure<
            IntegerCanonicalStructure,
            &'static IntegerCanonicalStructure,
        > {
            &self.integer_module_to_submodule
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
        IntegerSubmoduleInclusionB: BorrowedMorphism<
                IntegerSubmodule,
                IntegerModule,
                IntegerSubmoduleInclusion<
                    K,
                    IntegerModule,
                    IntegerModuleB,
                    IntegerSubmodule,
                    IntegerSubmoduleB,
                >,
            >,
    >
        Function<
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
        >
        for IntegerSubmoduleIntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
            IntegerSubmoduleInclusionB,
        >
    {
        fn image(&self, x: &FinitelyFreeSubmodule<Integer>) -> FinitelyFreeSubmodule<Integer> {
            self.integer_module_to_submodule().span(
                x.basis()
                    .into_iter()
                    .map(|bv| self.integer_submodule_to_module().image(&bv))
                    .collect(),
            )
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerModuleB: BorrowedStructure<IntegerModule>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
        IntegerSubmoduleInclusionB: BorrowedMorphism<
                IntegerSubmodule,
                IntegerModule,
                IntegerSubmoduleInclusion<
                    K,
                    IntegerModule,
                    IntegerModuleB,
                    IntegerSubmodule,
                    IntegerSubmoduleB,
                >,
            >,
    >
        InjectiveFunction<
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
            FinitelyFreeSubmoduleStructure<
                IntegerCanonicalStructure,
                &'static IntegerCanonicalStructure,
            >,
        >
        for IntegerSubmoduleIntegerSubmoduleInclusion<
            K,
            IntegerModule,
            IntegerModuleB,
            IntegerSubmodule,
            IntegerSubmoduleB,
            IntegerSubmoduleInclusionB,
        >
    {
        fn try_preimage(
            &self,
            y: &FinitelyFreeSubmodule<Integer>,
        ) -> Option<FinitelyFreeSubmodule<Integer>> {
            Some(
                self.integer_submodule_to_submodules().span(
                    y.basis()
                        .into_iter()
                        .map(|bv| self.integer_submodule_to_module().try_preimage(&bv))
                        .collect::<Option<Vec<_>>>()?,
                ),
            )
        }
    }
}

pub trait FullRankIntegerSubmoduleWithBasisSignature<K: AlgebraicNumberFieldSignature>:
    AdditiveGroupSignature<Set = Vec<Integer>>
{
    fn anf(&self) -> &K;

    fn basis(&self) -> &Vec<K::Set>;

    fn basis_vector(&self, i: usize) -> &K::Set {
        debug_assert!(i < self.n());
        self.basis().get(i).unwrap()
    }

    fn n(&self) -> usize {
        debug_assert_eq!(self.anf().n(), self.basis().len());
        self.anf().n()
    }

    fn free_integer_submodule_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, &'static IntegerCanonicalStructure>
    {
        Integer::structure_ref().free_module(self.n())
    }

    fn contains_element(&self, p: &K::Set) -> bool {
        self.outbound_order_to_anf_inclusion()
            .try_preimage(p)
            .is_some()
    }

    fn contains_integer_submodule<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >(
        &self,
        integer_submodule: &IntegerSubmodule,
    ) -> bool {
        integer_submodule
            .basis()
            .iter()
            .all(|sublat_basis_vector| self.contains_element(sublat_basis_vector))
    }

    fn discriminant(&self) -> Rational {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(self.basis())
    }

    fn into_outbound_order_to_anf_inclusion(
        self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, Self, Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion::new(self)
    }

    fn outbound_order_to_anf_inclusion(
        &self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, Self, &Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion::new(self)
    }

    fn into_inbound_integer_submodule_inclusion<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    >(
        self,
        integer_submodule: IntegerSubmoduleB,
    ) -> Option<
        integer_submodule_inclusion::IntegerSubmoduleInclusion<
            K,
            Self,
            Self,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >,
    > {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new(self, integer_submodule)
    }

    fn inbound_integer_submodule_inclusion<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    >(
        &self,
        integer_submodule: IntegerSubmoduleB,
    ) -> Option<
        integer_submodule_inclusion::IntegerSubmoduleInclusion<
            K,
            Self,
            &Self,
            IntegerSubmodule,
            IntegerSubmoduleB,
        >,
    > {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new(self, integer_submodule)
    }

    fn into_inbound_integer_submodule_inclusion_unchecked<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    >(
        self,
        integer_submodule: IntegerSubmoduleB,
    ) -> integer_submodule_inclusion::IntegerSubmoduleInclusion<
        K,
        Self,
        Self,
        IntegerSubmodule,
        IntegerSubmoduleB,
    > {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new_unchecked(
            self,
            integer_submodule,
        )
    }

    fn inbound_integer_submodule_inclusion_unchecked<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmoduleB: BorrowedStructure<IntegerSubmodule>,
    >(
        &self,
        integer_submodule: IntegerSubmoduleB,
    ) -> integer_submodule_inclusion::IntegerSubmoduleInclusion<
        K,
        Self,
        &Self,
        IntegerSubmodule,
        IntegerSubmoduleB,
    > {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new_unchecked(
            self,
            integer_submodule,
        )
    }
}
