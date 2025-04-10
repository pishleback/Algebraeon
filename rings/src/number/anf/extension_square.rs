use super::number_field::*;
use super::ring_of_integers::*;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
struct RingOfIntegersSquare<
    ZR: RingHomomorphism<IntegerCanonicalStructure, RingOfIntegersWithIntegralBasisStructure>
        + InjectiveFunction<IntegerCanonicalStructure, RingOfIntegersWithIntegralBasisStructure>,
    QK: FiniteDimensionalFieldExtension<RationalCanonicalStructure, AlgebraicNumberFieldStructure>,
    ZQ: FieldOfFractionsInclusion<IntegerCanonicalStructure, RationalCanonicalStructure>,
    RK: RingHomomorphism<RingOfIntegersWithIntegralBasisStructure, AlgebraicNumberFieldStructure>
        + InjectiveFunction<RingOfIntegersWithIntegralBasisStructure, AlgebraicNumberFieldStructure>,
> {
    z_to_r: ZR,
    q_to_k: QK,
    z_to_q: ZQ,
    r_to_k: RK,
}

impl AlgebraicNumberFieldStructure {
    fn ring_of_integers_commuting_square(
        &self,
    ) -> RingOfIntegersSquare<
        PrincipalSubringInclusion<RingOfIntegersWithIntegralBasisStructure>,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldStructure>,
        PrincipalSubringInclusion<RationalCanonicalStructure>,
        RingOfIntegersToAlgebraicNumberFieldInclusion,
    > {
        let anf = self;
        let roi = self.ring_of_integers();
        RingOfIntegersSquare {
            z_to_r: PrincipalSubringInclusion::new(roi.clone()),
            q_to_k: PrincipalRationalSubfieldInclusion::new(anf.clone()),
            z_to_q: PrincipalSubringInclusion::new(Rational::structure()),
            r_to_k: RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(roi),
        }
    }
}

// impl<
//     Z: IntegralDomainStructure,
//     R: IntegralDomainStructure,
//     Q: FieldStructure,
//     K: FieldStructure,
//     ZR: RingHomomorphism<Z, R> + InjectiveFunction<Z, R>,
//     QK: FiniteDimensionalFieldExtension<Q, K>,
//     ZQ: FieldOfFractionsInclusion<Z, Q>,
//     RK: RingHomomorphism<R, K> + InjectiveFunction<R, K>,
// > RingOfIntegersSquare<Z, R, Q, K, ZR, QK, ZQ, RK>
// {
//     pub fn new(z_to_r: ZR, q_to_k: QK, z_to_q: ZQ, r_to_k: RK) -> Self {
//         assert_eq!(z_to_r.domain(), z_to_q.domain());
//         assert_eq!(q_to_k.domain(), z_to_q.range());
//         assert_eq!(z_to_r.range(), r_to_k.domain());
//         assert_eq!(q_to_k.range(), r_to_k.range());
//         Self {
//             z: PhantomData,
//             r: PhantomData,
//             q: PhantomData,
//             k: PhantomData,
//             z_to_r,
//             q_to_k,
//             z_to_q,
//             r_to_k,
//         }
//     }
// }
