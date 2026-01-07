use std::marker::PhantomData;

use crate::structure::{AdditiveMonoidSignature, MultiplicativeMonoidWithZeroSignature};
use algebraeon_sets::structure::{EqSignature, OrdSignature};

pub struct FactoringStructure<
    Unfactored: MultiplicativeMonoidWithZeroSignature + EqSignature,
    Power: AdditiveMonoidSignature + OrdSignature,
> {
    _unfactored: PhantomData<Unfactored>,
    _power: PhantomData<Power>,
}
