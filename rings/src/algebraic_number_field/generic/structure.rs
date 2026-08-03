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
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use std::marker::PhantomData;
use std::sync::Arc;

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RationalInclusion: FiniteDimensionalFieldExtension<Self::Basis, RationalCanonicalStructure, Self>;

    fn inbound_finite_dimensional_rational_extension(
        self: &Arc<Self>,
    ) -> Arc<Self::RationalInclusion>;

    /// The dimension of this algebraic number field as a vector space over the rational numbers
    fn n(self: &Arc<Self>) -> usize {
        self.inbound_finite_dimensional_rational_extension()
            .degree()
    }

    /// An element which generates this algebraic number field when adjoined to the rational numbers
    /// Such an element always exists by the primitive element theorem
    fn generator(self: &Arc<Self>) -> Self::Elem;

    /// Determine whether an element is integral over the integers i.e. is it a root of a monic integer polynomial
    fn is_algebraic_integer(self: &Arc<Self>, a: &Self::Elem) -> bool;

    /// The discriminant of this algebraic number field i.e. the discriminant of its ring of integers
    /// Implementations should not compute this by constructing the ring of integers, as the constructor for a maximal OrderWithBasis calls this function to validate its input
    fn discriminant(self: &Arc<Self>) -> Integer;

    /// A list of self.n() elements which generate the ring of integers as a Z-module
    fn integral_basis(self: &Arc<Self>) -> Vec<Self::Elem>;

    fn ring_of_integers(self: &Arc<Self>) -> Arc<OrderWithBasis<Self, true>> {
        OrderWithBasis::new_maximal_unchecked(self.clone(), self.integral_basis())
    }

    fn order(
        self: &Arc<Self>,
        basis: Vec<Self::Elem>,
    ) -> Result<Arc<OrderWithBasis<Self, false>>, String> {
        OrderWithBasis::new(self.clone(), basis)
    }

    /// The LCM of the denominators of the coefficients of the minimal polynomial of a.
    ///
    /// It may well be >1 even when the element a is an algebraic integer.
    fn min_poly_denominator_lcm(self: &Arc<Self>, a: &Self::Elem) -> Integer {
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
    fn integral_multiple(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingSignature<K: AlgebraicNumberFieldSignature>:
    DedekindDomainSignature + CharZeroRingSignature
{
    fn n(self: &Arc<Self>) -> usize {
        self.anf().n()
    }

    fn anf(self: &Arc<Self>) -> Arc<K>;

    /// A list of self.n() elements which generate this ring as a Z-module
    fn integral_basis(self: &Arc<Self>) -> Vec<Self::Elem>;

    fn to_anf(self: &Arc<Self>, x: &Self::Elem) -> K::Elem;

    fn try_from_anf(self: &Arc<Self>, y: &K::Elem) -> Option<Self::Elem>;

    fn order(self: &Arc<Self>) -> Arc<OrderWithBasis<K, true>> {
        OrderWithBasis::new_maximal_unchecked(
            self.anf(),
            self.integral_basis()
                .into_iter()
                .map(|v| self.outbound_roi_to_anf_inclusion().image(&v))
                .collect(),
        )
    }

    fn outbound_roi_to_anf_inclusion(
        self: &Arc<Self>,
    ) -> Arc<RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self>> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self.clone())
    }

    fn inbound_order_inclusion<const MAXIMAL: bool>(
        self: &Arc<Self>,
        order: Arc<OrderWithBasis<K, MAXIMAL>>,
    ) -> Arc<order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion<K, Self, MAXIMAL>>
    {
        order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion::new(
            self.clone(),
            order,
        )
    }

    fn inbound_order_isomorphism(
        self: &Arc<Self>,
        order: Arc<OrderWithBasis<K, true>>,
    ) -> Arc<order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion<K, Self, true>>
    {
        order_to_ring_of_integers_inclusion::OrderToRingOfIntegersInclusion::new(
            self.clone(),
            order,
        )
    }
}

mod ring_of_integers_to_algebraic_number_field_inclusion {
    use crate::structure::IntegralClosureExtension;

    use super::*;

    #[derive(Debug, Clone)]
    pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
    > {
        _anf: PhantomData<K>,
        roi: Arc<R>,
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        pub fn from_ring_of_integers(roi: Arc<R>) -> Arc<Self> {
            Self {
                _anf: PhantomData,
                roi,
            }
            .into()
        }

        pub fn roi(&self) -> Arc<R> {
            self.roi.clone()
        }

        pub fn anf(&self) -> Arc<K> {
            self.roi().anf()
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>> Morphism<R, K>
        for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        fn domain(self: &Arc<Self>) -> Arc<R> {
            self.roi()
        }

        fn range(self: &Arc<Self>) -> Arc<K> {
            self.anf()
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        FunctionMorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        fn image(self: &Arc<Self>, x: &<R as SetSignature>::Elem) -> <K as SetSignature>::Elem {
            self.roi().to_anf(x)
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        InjectiveFunctionMorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        fn try_preimage(
            self: &Arc<Self>,
            y: &<K as SetSignature>::Elem,
        ) -> Option<<R as SetSignature>::Elem> {
            self.roi().try_from_anf(y)
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        RingHomomorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        FieldOfFractionsInclusion<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        fn numerator_and_denominator(self: &Arc<Self>, a: &<K>::Elem) -> (<R>::Elem, <R>::Elem) {
            self.zq_extension()
                .r_to_k_field_of_fractions()
                .numerator_and_denominator(a)
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R>
    {
        pub fn zq_extension(self: &Arc<Self>) -> Arc<RingOfIntegersIntegralExtension<K, R>> {
            RingOfIntegersIntegralExtension::new_integer_extension(
                RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self.domain()),
            )
        }
    }
}
pub(crate) use ring_of_integers_to_algebraic_number_field_inclusion::RingOfIntegersToAlgebraicNumberFieldInclusion;

mod order_to_ring_of_integers_inclusion {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct OrderToRingOfIntegersInclusion<
        K: AlgebraicNumberFieldSignature,
        R: AlgebraicIntegerRingSignature<K>,
        const MAXIMAL: bool,
    > {
        roi: Arc<R>,
        order: Arc<OrderWithBasis<K, MAXIMAL>>,
        order_basis_in_roi: Vec<R::Elem>,
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>, const MAXIMAL: bool>
        OrderToRingOfIntegersInclusion<K, R, MAXIMAL>
    {
        pub fn new(roi: Arc<R>, order: Arc<OrderWithBasis<K, MAXIMAL>>) -> Arc<Self> {
            let order_basis_in_roi = order
                .basis()
                .iter()
                .map(|bv| {
                    roi.outbound_roi_to_anf_inclusion()
                        .try_preimage(bv)
                        .unwrap()
                })
                .collect();
            Self {
                roi,
                order,
                order_basis_in_roi,
            }
            .into()
        }

        pub fn n(&self) -> usize {
            self.anf().n()
        }

        pub fn anf(&self) -> Arc<K> {
            debug_assert_eq!(self.roi().anf(), self.order().anf());
            self.order().anf()
        }

        pub fn roi(&self) -> &Arc<R> {
            &self.roi
        }

        pub fn order(&self) -> &Arc<OrderWithBasis<K, MAXIMAL>> {
            &self.order
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>, const MAXIMAL: bool>
        Morphism<OrderWithBasis<K, MAXIMAL>, R> for OrderToRingOfIntegersInclusion<K, R, MAXIMAL>
    {
        fn domain(self: &Arc<Self>) -> Arc<OrderWithBasis<K, MAXIMAL>> {
            self.order().clone()
        }

        fn range(self: &Arc<Self>) -> Arc<R> {
            self.roi().clone()
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>, const MAXIMAL: bool>
        FunctionMorphism<OrderWithBasis<K, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, MAXIMAL>
    {
        fn image(self: &Arc<Self>, x: &Vec<Integer>) -> <R as SetSignature>::Elem {
            self.roi().sum(
                &(0..self.n())
                    .map(|i| {
                        self.roi()
                            .mul(&self.roi().from_int(&x[i]), &self.order_basis_in_roi[i])
                    })
                    .collect::<Vec<_>>(),
            )
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>, const MAXIMAL: bool>
        RingHomomorphism<OrderWithBasis<K, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, MAXIMAL>
    {
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>, const MAXIMAL: bool>
        InjectiveFunctionMorphism<OrderWithBasis<K, MAXIMAL>, R>
        for OrderToRingOfIntegersInclusion<K, R, MAXIMAL>
    {
        fn try_preimage(self: &Arc<Self>, y: &<R as SetSignature>::Elem) -> Option<Vec<Integer>> {
            self.order()
                .outbound_order_to_anf_inclusion()
                .try_preimage(&self.roi().outbound_roi_to_anf_inclusion().image(y))
        }
    }

    impl<K: AlgebraicNumberFieldSignature, R: AlgebraicIntegerRingSignature<K>>
        BijectiveFunctionMorphism<OrderWithBasis<K, true>, R>
        for OrderToRingOfIntegersInclusion<K, R, true>
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
    > {
        _k: PhantomData<K>,
        integer_submodule: Arc<IntegerModule>,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, IntegerModule>
    {
        pub fn new(integer_submodule: Arc<IntegerModule>) -> Arc<Self> {
            Self {
                _k: PhantomData,
                integer_submodule,
            }
            .into()
        }

        pub fn integer_submodule(&self) -> &Arc<IntegerModule> {
            &self.integer_submodule
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > Morphism<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, IntegerModule>
    {
        fn domain(self: &Arc<Self>) -> Arc<IntegerModule> {
            self.integer_submodule().clone()
        }

        fn range(self: &Arc<Self>) -> Arc<K> {
            self.integer_submodule().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > FunctionMorphism<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, IntegerModule>
    {
        fn image(self: &Arc<Self>, x: &Vec<Integer>) -> <K as SetSignature>::Elem {
            debug_assert!(self.integer_submodule().validate_element(x).is_ok());
            let k = self.integer_submodule().anf();
            let n = k.n();
            debug_assert_eq!(n, x.len());
            k.sum(
                &(0..n)
                    .map(|i| k.mul(&k.from_int(&x[i]), self.integer_submodule().basis_vector(i)))
                    .collect::<Vec<_>>(),
            )
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > InjectiveFunctionMorphism<IntegerModule, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, IntegerModule>
    {
        fn try_preimage(self: &Arc<Self>, y: &<K as SetSignature>::Elem) -> Option<Vec<Integer>> {
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

    impl<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool>
        RingHomomorphism<OrderWithBasis<K, MAXIMAL>, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, MAXIMAL>,
        >
    {
    }

    impl<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool>
        FieldOfFractionsInclusion<OrderWithBasis<K, MAXIMAL>, K>
        for AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, MAXIMAL>,
        >
    {
        fn numerator_and_denominator(
            self: &Arc<Self>,
            a: &<K>::Elem,
        ) -> (
            <OrderWithBasis<K, MAXIMAL> as SetSignature>::Elem,
            <OrderWithBasis<K, MAXIMAL> as SetSignature>::Elem,
        ) {
            self.zq_extension()
                .r_to_k_field_of_fractions()
                .numerator_and_denominator(a)
        }
    }

    impl<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool>
        AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
            K,
            OrderWithBasis<K, MAXIMAL>,
        >
    {
        pub fn zq_extension(
            self: &Arc<Self>,
        ) -> Arc<order_integral_extension::OrderIntegralExtension<K, MAXIMAL>> {
            order_integral_extension::OrderIntegralExtension::new_integer_extension(self.clone())
        }
    }

    mod order_integral_extension {
        use super::*;
        use crate::{
            algebraic_number_field::OrderIdealsStructure,
            num_theory::integer_ideal::IntegerIdealsStructure,
            structure::{
                IdealsSignature, IntegralClosureExtension, PrincipalIntegerMap, RingSignature,
                RingToIdealsSignature,
            },
        };

        /// Q -> K
        /// ↑    ↑
        /// Z -> R
        ///
        /// Where Q is the rationals, Z is the integers, K is an algebraic number field, R is its ring of integers
        ///
        #[derive(Debug, Clone)]
        pub struct OrderIntegralExtension<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool> {
            r_to_k: Arc<
                AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                    K,
                    OrderWithBasis<K, MAXIMAL>,
                >,
            >,
        }

        impl<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool> OrderIntegralExtension<K, MAXIMAL> {
            pub fn new_integer_extension(
                r_to_k: Arc<
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, MAXIMAL>,
                    >,
                >,
            ) -> Arc<Self> {
                Self { r_to_k }.into()
            }

            pub fn with_ideals(
                &self,
            ) -> Arc<
                OrderIntegralExtensionWithIdeals<
                    K,
                    MAXIMAL,
                    IntegerIdealsStructure,
                    OrderIdealsStructure<K, MAXIMAL>,
                >,
            > {
                let ideals_z = Integer::structure().ideals();
                let ideals_r = self.r_to_k.domain().ideals();
                OrderIntegralExtensionWithIdeals::new_integer_extension(
                    self.r_to_k.clone(),
                    ideals_z,
                    ideals_r,
                )
            }
        }

        impl<K: AlgebraicNumberFieldSignature, const MAXIMAL: bool> IntegralClosureExtension
            for OrderIntegralExtension<K, MAXIMAL>
        {
            type QKBasis = K::Basis;
            type Z = IntegerCanonicalStructure;
            type Q = RationalCanonicalStructure;
            type R = OrderWithBasis<K, MAXIMAL>;
            type K = K;
            type ZQ = PrincipalIntegerMap<RationalCanonicalStructure>;
            type ZR = PrincipalIntegerMap<OrderWithBasis<K, MAXIMAL>>;
            type QK = K::RationalInclusion;
            type RK = AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                K,
                OrderWithBasis<K, MAXIMAL>,
            >;

            fn z_ring(self: &Arc<Self>) -> Arc<Self::Z> {
                Integer::structure()
            }
            fn r_ring(self: &Arc<Self>) -> Arc<Self::R> {
                self.r_to_k.domain()
            }
            fn q_field(self: &Arc<Self>) -> Arc<Self::Q> {
                Rational::structure()
            }
            fn k_field(self: &Arc<Self>) -> Arc<Self::K> {
                self.r_to_k.range()
            }

            fn z_to_q(self: &Arc<Self>) -> Arc<Self::ZQ> {
                Rational::structure().inbound_principal_integer_map()
            }
            fn z_to_r(self: &Arc<Self>) -> Arc<Self::ZR> {
                self.r_ring().inbound_principal_integer_map()
            }
            fn q_to_k(self: &Arc<Self>) -> Arc<Self::QK> {
                self.k_field()
                    .inbound_finite_dimensional_rational_extension()
            }
            fn r_to_k(self: &Arc<Self>) -> Arc<Self::RK> {
                self.r_to_k.clone()
            }

            fn integralize_multiplier(
                self: &Arc<Self>,
                alpha: &<Self::K as SetSignature>::Elem,
            ) -> Integer {
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
            const MAXIMAL: bool,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure>,
            IdealsR: IdealsSignature<OrderWithBasis<K, MAXIMAL>>,
        > {
            r_to_k: Arc<
                AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                    K,
                    OrderWithBasis<K, MAXIMAL>,
                >,
            >,
            ideals_z: Arc<IdealsZ>,
            ideals_r: Arc<IdealsR>,
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            const MAXIMAL: bool,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure>,
            IdealsR: IdealsSignature<OrderWithBasis<K, MAXIMAL>>,
        > OrderIntegralExtensionWithIdeals<K, MAXIMAL, IdealsZ, IdealsR>
        {
            pub fn new_integer_extension(
                r_to_k: Arc<
                    AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                        K,
                        OrderWithBasis<K, MAXIMAL>,
                    >,
                >,
                ideals_z: Arc<IdealsZ>,
                ideals_r: Arc<IdealsR>,
            ) -> Arc<Self> {
                Self {
                    r_to_k,
                    ideals_z,
                    ideals_r,
                }
                .into()
            }

            pub fn z_ideals(self: &Arc<Self>) -> &Arc<IdealsZ> {
                &self.ideals_z
            }

            pub fn r_ideals(self: &Arc<Self>) -> &Arc<IdealsR> {
                &self.ideals_r
            }
        }

        impl<
            K: AlgebraicNumberFieldSignature,
            const MAXIMAL: bool,
            IdealsZ: IdealsSignature<IntegerCanonicalStructure>,
            IdealsR: IdealsSignature<OrderWithBasis<K, MAXIMAL>>,
        > IntegralClosureExtension
            for OrderIntegralExtensionWithIdeals<K, MAXIMAL, IdealsZ, IdealsR>
        {
            type QKBasis = K::Basis;
            type Z = IntegerCanonicalStructure;
            type Q = RationalCanonicalStructure;
            type R = OrderWithBasis<K, MAXIMAL>;
            type K = K;
            type ZQ = PrincipalIntegerMap<RationalCanonicalStructure>;
            type ZR = PrincipalIntegerMap<OrderWithBasis<K, MAXIMAL>>;
            type QK = K::RationalInclusion;
            type RK = AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<
                K,
                OrderWithBasis<K, MAXIMAL>,
            >;

            fn z_ring(self: &Arc<Self>) -> Arc<Self::Z> {
                Integer::structure()
            }
            fn r_ring(self: &Arc<Self>) -> Arc<Self::R> {
                self.r_to_k.domain()
            }
            fn q_field(self: &Arc<Self>) -> Arc<Self::Q> {
                Rational::structure()
            }
            fn k_field(self: &Arc<Self>) -> Arc<Self::K> {
                self.r_to_k.range()
            }

            fn z_to_q(self: &Arc<Self>) -> Arc<Self::ZQ> {
                Rational::structure().inbound_principal_integer_map()
            }
            fn z_to_r(self: &Arc<Self>) -> Arc<Self::ZR> {
                self.r_ring().inbound_principal_integer_map()
            }
            fn q_to_k(self: &Arc<Self>) -> Arc<Self::QK> {
                self.k_field()
                    .inbound_finite_dimensional_rational_extension()
            }
            fn r_to_k(self: &Arc<Self>) -> Arc<Self::RK> {
                self.r_to_k.clone()
            }

            fn integralize_multiplier(
                self: &Arc<Self>,
                alpha: &<Self::K as SetSignature>::Elem,
            ) -> Integer {
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
        linear::finitely_free_submodules::{
            FinitelyFreeSubmodule, FinitelyFreeSubmodulesStructure,
        },
        structure::{FinitelyFreeModuleSignature, RingSignature},
    };
    use std::marker::PhantomData;

    #[derive(Debug, Clone)]
    pub struct IntegerSubmoduleInclusion<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > {
        _k: PhantomData<K>,
        integer_module: Arc<IntegerModule>,
        integer_submodule: Arc<IntegerSubmodule>,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        pub fn integer_module(&self) -> &Arc<IntegerModule> {
            &self.integer_module
        }

        pub fn integer_submodule(&self) -> &Arc<IntegerSubmodule> {
            &self.integer_submodule
        }

        fn check(self: &Arc<Self>) -> Result<(), String> {
            if self
                .integer_module()
                .contains_integer_submodule(self.integer_submodule())
            {
                Ok(())
            } else {
                Err("Integer module does not contain submodule".to_string())
            }
        }

        fn new_impl(
            integer_module: Arc<IntegerModule>,
            integer_submodule: Arc<IntegerSubmodule>,
        ) -> Arc<Self> {
            Self {
                _k: PhantomData,
                integer_module,
                integer_submodule,
            }
            .into()
        }

        pub fn new(
            integer_module: Arc<IntegerModule>,
            integer_submodule: Arc<IntegerSubmodule>,
        ) -> Option<Arc<Self>> {
            let s = Self::new_impl(integer_module, integer_submodule);
            if s.check().is_ok() { Some(s) } else { None }
        }

        pub fn new_unchecked(
            integer_module: Arc<IntegerModule>,
            integer_submodule: Arc<IntegerSubmodule>,
        ) -> Arc<Self> {
            let s = Self::new_impl(integer_module, integer_submodule);
            #[cfg(debug_assertions)]
            s.check().unwrap();
            s
        }

        pub fn integer_submodules_inclusion(
            self: &Arc<Self>,
        ) -> Arc<IntegerSubmoduleIntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>>
        {
            IntegerSubmoduleIntegerSubmoduleInclusion::new(self.clone())
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > Morphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn domain(self: &Arc<Self>) -> Arc<IntegerSubmodule> {
            self.integer_submodule().clone()
        }

        fn range(self: &Arc<Self>) -> Arc<IntegerModule> {
            self.integer_module().clone()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K> + RingSignature,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K> + RingSignature,
    > RingHomomorphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > FunctionMorphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn image(
            self: &Arc<Self>,
            x: &<IntegerSubmodule as SetSignature>::Elem,
        ) -> <IntegerModule as SetSignature>::Elem {
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
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > InjectiveFunctionMorphism<IntegerSubmodule, IntegerModule>
        for IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn try_preimage(
            self: &Arc<Self>,
            y: &<IntegerModule as SetSignature>::Elem,
        ) -> Option<<IntegerSubmodule as SetSignature>::Elem> {
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
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > {
        _k: PhantomData<K>,
        integer_submodule_to_module:
            Arc<IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>>,
        integer_module_to_submodule: Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        >,
        integer_submodule_to_submodules: Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        >,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    > IntegerSubmoduleIntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        pub fn new(
            integer_submodule_to_module: Arc<
                IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>,
            >,
        ) -> Arc<Self> {
            Self {
                _k: PhantomData,
                integer_module_to_submodule: integer_submodule_to_module
                    .integer_module()
                    .free_integer_submodule_restructure()
                    .submodules(),
                integer_submodule_to_submodules: integer_submodule_to_module
                    .integer_submodule()
                    .free_integer_submodule_restructure()
                    .submodules(),
                integer_submodule_to_module,
            }
            .into()
        }

        pub fn integer_submodule_to_module(
            &self,
        ) -> &Arc<IntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>> {
            &self.integer_submodule_to_module
        }

        pub fn integer_module_to_submodule(
            self: &Arc<Self>,
        ) -> &Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > {
            &self.integer_module_to_submodule
        }

        pub fn integer_submodule_to_submodules(
            self: &Arc<Self>,
        ) -> &Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > {
            &self.integer_submodule_to_submodules
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >
        Morphism<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > for IntegerSubmoduleIntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn domain(
            self: &Arc<Self>,
        ) -> Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > {
            self.integer_submodule_to_submodules.clone()
        }

        fn range(
            self: &Arc<Self>,
        ) -> Arc<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > {
            self.integer_module_to_submodule.clone()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        IntegerModule: FullRankIntegerSubmoduleWithBasisSignature<K>,
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >
        FunctionMorphism<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > for IntegerSubmoduleIntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn image(
            self: &Arc<Self>,
            x: &FinitelyFreeSubmodule<Integer>,
        ) -> FinitelyFreeSubmodule<Integer> {
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
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >
        InjectiveFunctionMorphism<
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
            FinitelyFreeSubmodulesStructure<
                EnumeratedFiniteSetStructure,
                IntegerCanonicalStructure,
                FinitelyFreeModuleStructure<
                    EnumeratedFiniteSetStructure,
                    IntegerCanonicalStructure,
                >,
            >,
        > for IntegerSubmoduleIntegerSubmoduleInclusion<K, IntegerModule, IntegerSubmodule>
    {
        fn try_preimage(
            self: &Arc<Self>,
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
    AdditiveGroupSignature<Elem = Vec<Integer>>
{
    fn anf(self: &Arc<Self>) -> Arc<K>;

    fn basis(self: &Arc<Self>) -> &Vec<K::Elem>;

    fn basis_vector(self: &Arc<Self>, i: usize) -> &K::Elem {
        debug_assert!(i < self.n());
        self.basis().get(i).unwrap()
    }

    fn n(self: &Arc<Self>) -> usize {
        debug_assert_eq!(self.anf().n(), self.basis().len());
        self.anf().n()
    }

    fn free_integer_submodule_restructure(
        self: &Arc<Self>,
    ) -> Arc<FinitelyFreeModuleStructure<EnumeratedFiniteSetStructure, IntegerCanonicalStructure>>
    {
        Integer::structure().free_module(EnumeratedFiniteSetStructure::new(self.n()))
    }

    fn contains_element(self: &Arc<Self>, p: &K::Elem) -> bool {
        self.outbound_order_to_anf_inclusion()
            .try_preimage(p)
            .is_some()
    }

    fn contains_integer_submodule<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >(
        self: &Arc<Self>,
        integer_submodule: &Arc<IntegerSubmodule>,
    ) -> bool {
        integer_submodule
            .basis()
            .iter()
            .all(|sublat_basis_vector| self.contains_element(sublat_basis_vector))
    }

    fn discriminant(self: &Arc<Self>) -> Rational {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(self.basis())
    }

    fn outbound_order_to_anf_inclusion(
        self: &Arc<Self>,
    ) -> Arc<anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion<K, Self>>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankIntegerSubmoduleWithBasisInclusion::new(
            self.clone(),
        )
    }

    fn inbound_integer_submodule_inclusion<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >(
        self: &Arc<Self>,
        integer_submodule: Arc<IntegerSubmodule>,
    ) -> Option<
        Arc<integer_submodule_inclusion::IntegerSubmoduleInclusion<K, Self, IntegerSubmodule>>,
    > {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new(self.clone(), integer_submodule)
    }

    fn inbound_integer_submodule_inclusion_unchecked<
        IntegerSubmodule: FullRankIntegerSubmoduleWithBasisSignature<K>,
    >(
        self: &Arc<Self>,
        integer_submodule: Arc<IntegerSubmodule>,
    ) -> Arc<integer_submodule_inclusion::IntegerSubmoduleInclusion<K, Self, IntegerSubmodule>>
    {
        integer_submodule_inclusion::IntegerSubmoduleInclusion::new_unchecked(
            self.clone(),
            integer_submodule,
        )
    }
}
