use std::{
    collections::{HashMap, HashSet},
    marker::PhantomData,
};

use super::*;

#[derive(Debug, Clone)]
struct SimplicialComplexSimplex {
    points: Vec<usize>, //non-empty ordered distinct
}

impl SimplicialComplexSimplex {
    pub fn new(mut points: Vec<usize>) -> Self {
        debug_assert!(!points.is_empty());
        points.sort_unstable();
        for i in 0..(points.len() - 1) {
            if points[i] == points[i + 1] {
                panic!("Points must be distinct");
            }
        }
        Self { points }
    }
}

#[derive(Debug, Clone)]
pub struct SimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
> {
    ambient_space: SP,
    points: Vec<Vector<FS, SP>>,
    simplexes: Vec<SimplicialComplexSimplex>,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialComplex<FS, SP>
where
    FS::Set: Hash,
{
    fn new_impl(ambient_space: SP, simplexes: Vec<Simplex<FS, SP>>) -> Self {
        let mut idx_to_point: Vec<Vector<FS, SP>> = vec![];
        let mut existing_points: HashMap<Vector<FS, SP>, usize> = HashMap::new();

        let mut simplexes_as_idxs = vec![];
        for simplex in &simplexes {
            let mut simplex_points = vec![];
            for point in simplex.points() {
                simplex_points.push(match existing_points.get(&point) {
                    Some(idx) => *idx,
                    None => {
                        let idx = idx_to_point.len();
                        idx_to_point.push(point.clone());
                        existing_points.insert(point.clone(), idx);
                        idx
                    }
                });
            }
            simplexes_as_idxs.push(SimplicialComplexSimplex::new(simplex_points));
        }
        Self {
            ambient_space,
            points: idx_to_point,
            simplexes: simplexes_as_idxs,
        }
    }

    pub fn new_unchecked(ambient_space: SP, simplexes: Vec<Simplex<FS, SP>>) -> Self {
        #[cfg(debug_assertions)]
        Self::new(ambient_space.clone(), simplexes.clone()).unwrap();
        Self::new_impl(ambient_space, simplexes)
    }

    pub fn new(ambient_space: SP, simplexes: Vec<Simplex<FS, SP>>) -> Result<Self, &'static str> {
        for simplex in &simplexes {
            assert_eq!(simplex.ambient_space().borrow(), ambient_space.borrow());
            if simplex.n() == 0 {
                return Err("Simplicial complex musn't contain the null simplex");
            }
        }
        let mut all_simplexes: HashSet<_> = simplexes.iter().collect();
        if simplexes.len() != all_simplexes.len() {
            return Err("Simplicial complex simplexes must be distinct");
        }
        for simplex in &simplexes {
            for bdry_simplex in simplex.sub_simplices_not_null() {
                if !all_simplexes.contains(&bdry_simplex) {
                    return Err("Simplicial complex must be closed under taking facets");
                }
            }
        }
        Ok(Self::new_impl(ambient_space, simplexes))
    }

    fn idx_to_simplex(&self, idx: usize) -> Simplex<FS, SP> {
        Simplex::new(
            self.ambient_space(),
            self.simplexes[idx]
                .points
                .iter()
                .map(|pt| self.points[*pt].clone())
                .collect(),
        )
        .unwrap()
    }

    pub fn simplexes(&self) -> Vec<Simplex<FS, SP>> {
        (0..self.simplexes.len())
            .map(|idx| self.idx_to_simplex(idx))
            .collect()
    }

    pub fn ambient_space(&self) -> SP {
        self.ambient_space.clone()
    }
}

#[derive(Debug, Clone)]
pub struct PartialSubSimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    SC: Borrow<SimplicialComplex<FS, SP>> + Clone,
> {
    ordered_field: PhantomData<FS>,
    ambient_space: PhantomData<SP>,
    simplicial_complex: SC,
    subset: HashSet<usize>,
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        SC: Borrow<SimplicialComplex<FS, SP>> + Clone,
    > PartialSubSimplicialComplex<FS, SP, SC>
where
    FS::Set: Hash,
{
    pub fn new_unchecked(simplicial_complex: SC, subset: HashSet<usize>) -> Self {
        #[cfg(debug_assertions)]
        Self::new(simplicial_complex.clone(), subset.clone()).unwrap();
        Self {
            ordered_field: PhantomData,
            ambient_space: PhantomData,
            simplicial_complex,
            subset,
        }
    }

    pub fn new(simplicial_complex: SC, subset: HashSet<usize>) -> Result<Self, &'static str> {
        for idx in &subset {
            if *idx > simplicial_complex.borrow().simplexes.len() {
                return Err("Simplicial complex subset simplex index out of range");
            }
        }
        Ok(Self {
            ordered_field: PhantomData,
            ambient_space: PhantomData,
            simplicial_complex,
            subset,
        })
    }

    pub fn simplexes(&self) -> HashSet<Simplex<FS, SP>> {
        self.subset
            .iter()
            .map(|idx| self.simplicial_complex.borrow().idx_to_simplex(*idx))
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct FullSubSimplicialComplex<
    FS: OrderedRingStructure + FieldStructure,
    SP: Borrow<AffineSpace<FS>> + Clone,
    SC: Borrow<SimplicialComplex<FS, SP>> + Clone,
> {
    ordered_field: PhantomData<FS>,
    ambient_space: PhantomData<SP>,
    simplicial_complex: SC,
    subset: HashSet<usize>,
}

impl<
        FS: OrderedRingStructure + FieldStructure,
        SP: Borrow<AffineSpace<FS>> + Clone,
        SC: Borrow<SimplicialComplex<FS, SP>> + Clone,
    > FullSubSimplicialComplex<FS, SP, SC>
where
    FS::Set: Hash,
{
    pub fn new_unchecked(simplicial_complex: SC, subset: HashSet<usize>) -> Self {
        #[cfg(debug_assertions)]
        Self::new(simplicial_complex.clone(), subset.clone()).unwrap();
        Self {
            ordered_field: PhantomData,
            ambient_space: PhantomData,
            simplicial_complex,
            subset,
        }
    }

    pub fn new(simplicial_complex: SC, subset: HashSet<usize>) -> Result<Self, &'static str> {
        let ambient_space = simplicial_complex.borrow().ambient_space();
        let partial = PartialSubSimplicialComplex::new(simplicial_complex, subset)?;
        let simplexes_subset = partial.simplexes();
        for spx in &simplexes_subset {
            for sub_spx in spx.sub_simplices_not_null() {
                if !(simplexes_subset.contains(&sub_spx)) {
                    return Err("Full sub simplicial complex must be closed under taking facets");
                }
            }
        }
        Ok(Self {
            ordered_field: PhantomData,
            ambient_space: PhantomData,
            simplicial_complex: partial.simplicial_complex,
            subset: partial.subset,
        })
    }

    pub fn simplexes(&self) -> HashSet<Simplex<FS, SP>> {
        self.subset
            .iter()
            .map(|idx| self.simplicial_complex.borrow().idx_to_simplex(*idx))
            .collect()
    }
}
