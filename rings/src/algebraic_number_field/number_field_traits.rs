use crate::structure::{FieldSignature, FiniteDimensionalFieldExtension};
use algebraeon_nzq::RationalCanonicalStructure;

pub trait AlgebraicNumberField<K: FieldSignature>:
    FiniteDimensionalFieldExtension<RationalCanonicalStructure, K>
{
}
