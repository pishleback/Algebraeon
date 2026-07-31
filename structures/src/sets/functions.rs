use crate::*;
use algebraeon_macros::{ref_meta, signature_meta_trait};

#[signature_meta_trait]
pub trait FunctionsSignature<Domain: SetSignature, Range: SetSignature>: SetSignature {
    fn function(&self, f: impl Fn(&Domain::Elem) -> Range::Elem) -> Option<Self::Elem>;
    fn image<'a>(&self, f: &'a Self::Elem, x: &Domain::Elem) -> &'a Range::Elem;
}

#[signature_meta_trait]
pub trait FunctionsDomainPermutationActionSignature<Domain: SetSignature, Range: SetSignature>:
    FunctionsSignature<Domain, Range>
{
    type DomainPerms: PermutationsSignature<Domain>;
    type DomainPermsRef<'a>: PermutationsSignature<Domain>
    where
        Self: 'a;

    fn into_domain_permutation_action(
        self,
    ) -> impl RightGroupActionSignature<Self, Self::DomainPerms>;

    #[ref_meta]
    fn domain_permutation_action<'a>(
        &'a self,
    ) -> impl RightGroupActionSignature<Self, Self::DomainPermsRef<'a>>;
}

#[signature_meta_trait]
pub trait FunctionsRangePermutationActionSignature<Domain: SetSignature, Range: SetSignature>:
    FunctionsSignature<Domain, Range>
{
    type RangePerms: PermutationsSignature<Range>;
    type RangePermsRef<'a>: PermutationsSignature<Range>
    where
        Self: 'a;

    fn into_range_permutation_action(self)
    -> impl LeftGroupActionSignature<Self::RangePerms, Self>;

    #[ref_meta]
    fn range_permutation_action<'a>(
        &'a self,
    ) -> impl LeftGroupActionSignature<Self::RangePermsRef<'a>, Self>;
}
