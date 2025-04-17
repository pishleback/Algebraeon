
use crate::structure::*;
use algebraeon_sets::structure::*;

pub trait ModuleSignature<Ring: RingSignature>: SetSignature {
    fn ring(&self) ->&Ring;
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn neg(&self, a: &Self::Set) -> Self::Set;
    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set;
}

pub trait FreeModuleSignature<Ring: RingSignature>: ModuleSignature<Ring> {}

pub trait FiniteRankModuleSignature<Ring: RingSignature>: FreeModuleSignature<Ring> {
    fn basis(&self) -> Vec<Self::Set>;
    fn rank(&self) -> usize {
        self.basis().len()
    }
}

pub trait ModuleHomomorphism<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: Function<Domain, Range>
{
}

// mod free_modules {
//     use crate::structure::*;

//     pub struct FreeModuleOverSetStructure<Ring: RingSignature> {
//         ring: Ring,
//     }

//     pub struct FreeModuleFiniteNumberedBasisStructure<Ring: RingSignature> {
//         ring: Ring,
//     }
// }

/*

// the free module of fixed finite rank over a ring
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteRankFreeModuleStructure<Ring: RingStructure> {
    ring: Ring,
    rank: usize,
}
pub type FiniteDimensionalVectorSpace<Field: FieldStructure> = FiniteRankFreeModuleStructure<Field>;

impl<Ring: RingStructure> Structure for FiniteRankFreeModuleStructure<Ring> {}

impl<Ring: RingStructure> SetStructure for FiniteRankFreeModuleStructure<Ring> {
    type Set = Vec<Ring::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.rank == x.len()
    }
}

impl<Ring: RingStructure> ModuleStructure<Ring> for FiniteRankFreeModuleStructure<Ring> {
    fn ring(&self) -> Ring {
        self.ring.clone()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        (0..self.rank).map(|i| self.ring.add(a[i], b[i])).collect()
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(a));
        a.iter().map(|x| self.ring.neg(x)).collect()
    }

    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set {
        debug_assert!(self.ring.is_element(x));
        debug_assert!(self.is_element(&a));
        a.iter().map(|y| self.ring.mul(x, y)).collect()
    }
}

impl<Ring: RingStructure> FreeFiniteRankModuleStructure<Ring>
    for FiniteRankFreeModuleStructure<Ring>
{
    fn rank(&self) -> usize {
        self.rank
    }
}

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct LinearTransformation<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: Ring,
    domain: Domain,
    range: Range,
    matrix: Matrix<Ring::Set>, // v -> Mv
}

impl<
    Ring: BezoutDomainStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> LinearTransformation<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
    pub fn new(ring: Ring, domain: Domain, range: Range, matrix: Matrix<Ring::Set>) -> Self {
        debug_assert_eq!(ring, domain.ring());
        debug_assert_eq!(ring, range.ring());
        debug_assert_eq!(domain.rank(), matrix.cols());
        debug_assert_eq!(range.rank(), matrix.rows());
        let rank = MatrixStructure::new(ring.clone()).rank(matrix.clone());
        if INJECTIVE {
            debug_assert_eq!(rank, domain.rank());
        }
        if SURJECTIVE {
            debug_assert_eq!(rank, range.rank());
        }
        Self {
            ring,
            domain,
            range,
            matrix,
        }
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> Morphism<Domain, Range> for LinearTransformation<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
    fn domain(&self) -> &Domain {
        &self.domain
    }

    fn range(&self) -> &Range {
        &self.range
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> Function<Domain, Range> for LinearTransformation<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
    fn image(&self, x: &Domain::Set) -> Range::Set {
        compile_error!("grrr")
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const SURJECTIVE: bool,
> InjectiveFunction<Domain, Range> for LinearTransformation<Ring, Domain, Range, true, SURJECTIVE>
{
    fn try_preimage(&self, x: &Range::Set) -> Option<<Domain as SetStructure>::Set> {
        compile_error!("grrr")
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> BijectiveFunction<Domain, Range> for LinearTransformation<Ring, Domain, Range, true, true>
{
    fn preimage(&self, x: &Range::Set) -> Domain::Set {
        compile_error!("grrr")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearTransformationStructure<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> {
    ring: Ring,
    domain: Domain,
    range: Range,
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> Structure for LinearTransformationStructure<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> SetStructure for LinearTransformationStructure<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
    type Set = LinearTransformation<Ring, Domain, Range, INJECTIVE, SURJECTIVE>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.ring == x.ring && self.domain == x.domain && self.range == x.range
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
    const INJECTIVE: bool,
    const SURJECTIVE: bool,
> MorphismsStructure<Domain, Range>
    for LinearTransformationStructure<Ring, Domain, Range, INJECTIVE, SURJECTIVE>
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tmp() {
        println!("hi");
    }
}
*/
