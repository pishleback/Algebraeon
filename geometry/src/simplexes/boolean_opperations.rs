use std::collections::HashSet;

use super::*;

pub struct VennResult<W, X> {
    pub left: W,
    pub middle: X,
    pub right: W,
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone> Simplex<FS, SP>
where
    FS::Set: Hash,
{
    pub fn venn(
        &self,
        other: &Self,
    ) -> VennResult<SimplicialDisjointUnion<FS, SP>, SimplicialDisjointUnion<FS, SP>> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        let overlap = ConvexHull::intersect(
            &ConvexHull::from_simplex(self.clone()),
            &ConvexHull::from_simplex(other.clone()),
        );

        let mut self_ext = overlap.clone();
        for pt in self.points() {
            self_ext.extend_by_point(pt.clone());
        }
        let self_parts = self_ext
            .as_simplicial_complex()
            .interior
            .simplexes()
            .into_iter()
            .map(|s| s.clone())
            .collect::<HashSet<_>>();

        let mut other_ext = overlap.clone();
        for pt in other.points() {
            other_ext.extend_by_point(pt.clone());
        }
        let other_parts = other_ext
            .as_simplicial_complex()
            .interior
            .simplexes()
            .into_iter()
            .map(|s| s.clone())
            .collect::<HashSet<_>>();

        let mut left = self_parts.clone();
        let mut middle = HashSet::new();
        let mut right = HashSet::new();
        for other_part in other_parts {
            if left.contains(&other_part) {
                left.remove(&other_part);
                middle.insert(other_part);
            } else {
                right.insert(other_part);
            }
        }

        if middle.len() == 0 {
            VennResult {
                left: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::from([self.clone()]),
                ),
                middle: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::new(),
                ),
                right: SimplicialDisjointUnion::new_unchecked(
                    ambient_space.clone(),
                    HashSet::from([other.clone()]),
                ),
            }
        } else {
            VennResult {
                left: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), left),
                middle: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), middle),
                right: SimplicialDisjointUnion::new_unchecked(ambient_space.clone(), right),
            }
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialDisjointUnion<FS, SP>
where
    FS::Set: Hash,
{
    pub fn subtract_raw(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        Self::new_unchecked(ambient_space.clone(), {
            let mut simplexes = HashSet::new();
            for self_spx in self.simplexes() {
                let mut self_leftover = HashSet::from([self_spx.clone()]);
                for other_spx in other.simplexes() {
                    self_leftover = self_leftover
                        .into_iter()
                        .map(|self_leftover_spx| {
                            Simplex::venn(&self_leftover_spx, other_spx)
                                .left
                                .into_simplexes()
                        })
                        .flatten()
                        .collect();
                }
                for spx in self_leftover {
                    simplexes.insert(spx);
                }
            }
            simplexes
        })
    }

    pub fn intersection_raw(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        Self::new_unchecked(ambient_space.clone(), {
            let mut simplexes = HashSet::new();
            for self_spx in self.simplexes() {
                for other_spx in other.simplexes() {
                    for spx in Simplex::venn(self_spx, other_spx).middle.into_simplexes() {
                        simplexes.insert(spx);
                    }
                }
            }
            simplexes
        })
    }

    pub fn union_raw(&self, other: &Self) -> SimplicialDisjointUnion<FS, SP> {
        let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();
        let mut simplexes = HashSet::new();
        for spx in Self::subtract_raw(other, self).into_simplexes() {
            simplexes.insert(spx);
        }
        for spx in self.simplexes() {
            simplexes.insert(spx.clone());
        }
        return Self::new_unchecked(ambient_space, simplexes);
    }

    pub fn venn_raw(
        &self,
        other: &Self,
    ) -> VennResult<SimplicialDisjointUnion<FS, SP>, SimplicialDisjointUnion<FS, SP>> {
        // let ambient_space = common_space(self.ambient_space(), other.ambient_space()).unwrap();

        VennResult {
            left: Self::subtract_raw(self, other),
            middle: Self::intersection_raw(self, other),
            right: Self::subtract_raw(other, self),
        }
    }
}

impl<FS: OrderedRingStructure + FieldStructure, SP: Borrow<AffineSpace<FS>> + Clone>
    SimplicialComplex<FS, SP>
where
    FS::Set: Hash,
{
    pub fn union(&self, other: &Self) -> Self {
        self.union_raw(other).simplify()
    }

    pub fn union_raw(&self, other: &Self) -> Self {
        SimplicialDisjointUnion::union_raw(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .closure_as_simplicial_complex()
    }

    pub fn intersection_raw(&self, other: &Self) -> Self {
        SimplicialDisjointUnion::intersection_raw(&self.into(), &other.into())
            .refine_to_partial_simplicial_complex()
            .closure_as_simplicial_complex()
    }

    pub fn intersection(&self, other: &Self) -> Self {
        self.intersection_raw(other).simplify()
    }
}
