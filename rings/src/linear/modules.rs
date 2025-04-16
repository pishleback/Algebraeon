use super::matrix::Matrix;
use crate::structure::*;
use algebraeon_sets::structure::*;

// modules and vector spaces
pub trait ModuleStructure<Ring: RingStructure>: SetStructure {
    fn ring(&self) -> Ring;
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn neg(&self, a: &Self::Set) -> Self::Set;
    fn scalar_mul(&self, x: &Ring::Set, a: Self::Set) -> Self::Set;
}

pub trait VectorSpaceStructure<Field: FieldStructure>: ModuleStructure<Field> {
    fn field(&self) -> Field;
}
impl<Field: FieldStructure, FMod: ModuleStructure<Field>> VectorSpaceStructure<Field> for FMod {
    fn field(&self) -> Field {
        self.ring()
    }
}

// modules of finite rank and vector spaces of finite dimension
pub trait FreeFiniteRankModuleStructure<Ring: RingStructure>: ModuleStructure<Ring> {
    fn rank(&self) -> usize;
}

pub trait FiniteDimensionalVectorSpaceStructure<Field: FieldStructure>:
    FreeFiniteRankModuleStructure<Field>
{
    fn dimension(&self) -> usize;
}
impl<Field: FieldStructure, FMod: FreeFiniteRankModuleStructure<Field>>
    FiniteDimensionalVectorSpaceStructure<Field> for FMod
{
    fn dimension(&self) -> usize {
        self.rank()
    }
}

// modules of finite rank and vector spaces of finite dimension with a prefered basis
pub trait FiniteRankModuleWithBasisStructure<Ring: RingStructure>:
    FreeFiniteRankModuleStructure<Ring>
{
    fn basis(&self) -> Vec<Self::Set>;
    fn coordinates(&self, a: &Self::Set) -> Vec<Ring::Set>;
    fn into_coordinates(&self, a: Self::Set) -> Vec<Ring::Set>;
    fn coordinate(&self, a: &Self::Set, i: usize) -> Ring::Set {
        self.coordinates(a).into_iter().nth(i).unwrap()
    }
    fn into_coordinate(&self, a: Self::Set, i: usize) -> Ring::Set {
        self.into_coordinates(a).into_iter().nth(i).unwrap()
    }
}

pub trait FiniteDimensionalVectorSpaceWithBasisStructure<Field: FieldStructure>:
    FiniteRankModuleWithBasisStructure<Field>
{
}
impl<Field: FieldStructure, FMod: FiniteRankModuleWithBasisStructure<Field>>
    FiniteDimensionalVectorSpaceWithBasisStructure<Field> for FMod
{
}

// linear maps of finite rank free modules with a basis
#[derive(Debug, Clone)]
pub struct LinearTransformation<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> {
    ring: Ring,
    domain: Domain,
    range: Range,
    matrix: Matrix<Ring::Set>,
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> LinearTransformation<Ring, Domain, Range>
{
    pub fn new(ring: Ring, domain: Domain, range: Range, matrix: Matrix<Ring::Set>) -> Self {
        debug_assert_eq!(ring, domain.ring());
        debug_assert_eq!(ring, range.ring());
        debug_assert_eq!(domain.rank(), matrix.cols());
        debug_assert_eq!(range.rank(), matrix.rows());
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
> Morphism<Domain, Range> for LinearTransformation<Ring, Domain, Range>
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
> Function<Domain, Range> for LinearTransformation<Ring, Domain, Range>
{
    fn image(&self, x: &Domain::Set) -> Range::Set {
        todo!()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearTransformationStructure<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> {
    ring: Ring,
    domain: Domain,
    range: Range,
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> Structure for LinearTransformationStructure<Ring, Domain, Range>
{
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> SetStructure for LinearTransformationStructure<Ring, Domain, Range>
{
    type Set = LinearTransformation<Ring, Domain, Range>;

    fn is_element(&self, x: &Self::Set) -> bool {
        true
    }
}

impl<
    Ring: RingStructure,
    Domain: FiniteRankModuleWithBasisStructure<Ring>,
    Range: FiniteRankModuleWithBasisStructure<Ring>,
> MorphismsStructure<Domain, Range> for LinearTransformationStructure<Ring, Domain, Range>
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
