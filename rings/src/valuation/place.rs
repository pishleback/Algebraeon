use algebraeon_nzq::Integer;
use algebraeon_sets::structure::{Function, SetSignature};

use super::AdditiveValuation;
use crate::{
    algebraic_number_field::{
        embedded_anf::EmbeddedAnf,
        polynomial_quotient_number_field::AlgebraicNumberFieldPolynomialQuotientStructure,
    },
    isolated_algebraic::ComplexAlgebraic,
    valuation::AdditiveValueGroup,
};

/// A place for an `AlgebraicNumberFieldStructure` can be
/// a `Valuation` to the integers
/// or an embedding specified by `EmbeddedAnf` which fixes
/// where the `generator` goes in the real or complex numbers.
/// The `Valuation` variant really represents it's equivalence class
/// of valuations
pub enum Place<V>
where
    V: AdditiveValuation<
            ValuationGamma = Integer,
            DomainFieldSignature = AlgebraicNumberFieldPolynomialQuotientStructure,
        >,
{
    Embedding(Box<EmbeddedAnf>),
    Valuation(V),
}

impl<V> Place<V>
where
    V: AdditiveValuation<
            ValuationGamma = Integer,
            DomainFieldSignature = AlgebraicNumberFieldPolynomialQuotientStructure,
        >,
{
    pub fn is_finite(&self) -> bool {
        matches!(self, Self::Valuation(_))
    }

    pub fn is_archemedian(&self) -> bool {
        matches!(self, Self::Embedding(_))
    }

    pub fn is_real(&self) -> bool {
        matches!(self, Self::Embedding(e) if e.is_real())
    }

    /// If we are at an Archimedean place, use the embedding to embed this `element`
    /// as a `ComplexAlgebraic`.
    fn embed_archimedian(
        &self,
        element: &<AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set,
    ) -> Result<ComplexAlgebraic, ()> {
        match self {
            Place::Embedding(embedded_anf) => Ok(embedded_anf.embed(element)),
            Place::Valuation(_) => Err(()),
        }
    }

    /// If we are at an non-Archimedean place and the element is in the corresponding `ValuationRing`
    /// give the quotient in `ResidueField`.
    /// If it is not in the corresponding `ValuationRing`, give the valuation which should be less than `0`.
    fn try_to_residue(
        &self,
        element: &<AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set,
    ) -> Result<
        Result<<V::ResidueField as SetSignature>::Set, AdditiveValueGroup<V::ValuationGamma>>,
        (),
    > {
        match self {
            Place::Embedding(_) => Err(()),
            Place::Valuation(v) => {
                let in_valuation_ring = match v.try_to_valuation_ring(element) {
                    Ok(z) => z,
                    Err(w) => {
                        return Ok(Err(w));
                    }
                };
                Ok(Ok(v.valuation_quotient().image(&in_valuation_ring)))
            }
        }
    }
}
