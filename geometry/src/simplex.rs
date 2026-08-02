use super::*;
use crate::{
    affine_subspace::EmbeddedAffineSubspace, ambient_space::AffineSpace,
    oriented_simplex::OrientedSimplex, vector::Vector,
};
use itertools::Itertools;

#[derive(Clone)]
pub struct Simplex<FS: OrderedRingSignature + FieldSignature> {
    ambient_space: AffineSpace<FS>,
    // points must be ordered w.r.t the ordering on vectors
    // points must be non-degenerate
    // points must belong to the ambient_space
    points: Vec<Vector<FS>>,
}

impl<FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for Simplex<FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simplex")
            .field("points", &self.points)
            .finish()
    }
}

impl<FS: OrderedRingSignature + FieldSignature> PartialEq for Simplex<FS> {
    fn eq(&self, other: &Self) -> bool {
        self.ambient_space == other.ambient_space && self.points == other.points
    }
}

impl<FS: OrderedRingSignature + FieldSignature> Eq for Simplex<FS> {}

impl<FS: OrderedRingSignature + FieldSignature> Hash for Simplex<FS>
where
    FS::Elem: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.points.hash(state);
    }
}

impl<FS: OrderedRingSignature + FieldSignature> Simplex<FS>
where
    AffineSpace<FS>: Clone,
{
    fn new(
        ambient_space: &AffineSpace<FS>,
        mut points: Vec<Vector<FS>>,
    ) -> Result<Self, &'static str> {
        for point in &points {
            assert_eq!(ambient_space, point.ambient_space());
        }
        points.sort_unstable();
        if ambient_space.are_points_affine_independent(points.iter().collect()) {
            Ok(Self {
                ambient_space: ambient_space.clone(),
                points,
            })
        } else {
            Err("Can't make a simplex using degenerate points")
        }
    }
}

impl<FS: OrderedRingSignature + FieldSignature> AffineSpace<FS>
where
    AffineSpace<FS>: Clone,
{
    pub fn simplex(&self, points: Vec<Vector<FS>>) -> Result<Simplex<FS>, &'static str> {
        Simplex::new(self, points)
    }
}

impl<FS: OrderedRingSignature + FieldSignature> Simplex<FS>
where
    AffineSpace<FS>: Clone,
{
    pub fn ambient_space(&self) -> &AffineSpace<FS> {
        &self.ambient_space
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn points(&self) -> &Vec<Vector<FS>> {
        &self.points
    }

    pub fn into_points(self) -> Vec<Vector<FS>> {
        self.points
    }

    pub fn point(&self, i: usize) -> &Vector<FS> {
        &self.points[i]
    }

    pub fn skeleton(&self, skel_n: isize) -> Vec<Self> {
        if skel_n < 0 {
            vec![]
        } else {
            let skel_n = skel_n as usize;
            let mut parts = vec![];
            for subset in (0..self.points.len()).combinations(skel_n) {
                let part = Self::new(
                    &self.ambient_space,
                    subset.into_iter().map(|i| self.points[i].clone()).collect(),
                )
                .unwrap();
                parts.push(part);
            }
            parts
        }
    }

    pub fn vertices(&self) -> Vec<Self> {
        self.skeleton(1)
    }

    pub fn edges(&self) -> Vec<Self> {
        self.skeleton(2)
    }

    pub fn faces(&self) -> Vec<Self> {
        self.skeleton(3)
    }

    pub fn ridges(&self) -> Vec<Self> {
        self.skeleton(self.points.len() as isize - 2)
    }

    pub fn facets(&self) -> Vec<Self> {
        self.skeleton(self.points.len() as isize - 1)
    }

    pub fn facet(&self, k: usize) -> Self {
        assert!(k <= self.points.len());
        Self::new(&self.ambient_space, {
            let mut facet_points = self.points.clone();
            facet_points.remove(k);
            facet_points
        })
        .unwrap()
    }

    pub fn sub_simplices(&self) -> Vec<Self> {
        self.points()
            .clone()
            .into_iter()
            .powerset()
            .map(|sub_points| Self::new(&self.ambient_space, sub_points).unwrap())
            .collect()
    }

    pub fn sub_simplices_not_null(&self) -> Vec<Self> {
        self.sub_simplices()
            .into_iter()
            .filter(|spx| spx.n() != 0)
            .collect()
    }

    pub fn proper_sub_simplices_not_null(&self) -> Vec<Self> {
        self.sub_simplices()
            .into_iter()
            .filter(|spx| spx.n() != 0 && spx.n() != self.n())
            .collect()
    }

    pub fn sub_simplex(&self, pts: Vec<usize>) -> Self {
        Self::new(
            self.ambient_space(),
            pts.iter().map(|idx| self.points[*idx].clone()).collect(),
        )
        .unwrap()
    }

    pub fn oriented_facet(&self, k: usize) -> OrientedSimplex<FS> {
        //return the oriented facet of self with negative side on the outside and positive side on the inside
        assert!(k <= self.points.len());
        let mut facet_points = self.points.clone();
        let other_pt = facet_points.remove(k);
        OrientedSimplex::new_with_positive_point(&self.ambient_space, facet_points, &other_pt)
            .unwrap()
    }

    pub fn oriented_facets(&self) -> Vec<OrientedSimplex<FS>> {
        assert_eq!(self.ambient_space.borrow().affine_dimension(), self.n());
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }

    pub fn into_affine_span(self) -> (EmbeddedAffineSubspace<FS>, Vec<Vector<FS>>) {
        EmbeddedAffineSubspace::new_affine_independent_span(
            &self.ambient_space,
            self.clone().into_points(),
        )
        .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_structures::{MetaType, Rational};

    #[test]
    fn make_simplex() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = space.vector([1, 1]);
        let v2 = space.vector([1, 0]);
        let v3 = space.vector([0, 1]);
        let s = Simplex::new(&space, vec![v1, v2, v3]);
        assert!(s.is_ok());

        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = space.vector([0, 0]);
        let v2 = space.vector([1, 0]);
        let v3 = space.vector([2, 0]);
        let s = Simplex::new(&space, vec![v1, v2, v3]);
        assert!(s.is_err());
    }

    #[test]
    fn simplex_skeleton() {
        let space = AffineSpace::new_linear(Rational::structure(), 2);
        let v1 = space.vector([1, 1]);
        let v2 = space.vector([1, 0]);
        let v3 = space.vector([0, 1]);
        let s = Simplex::new(&space, vec![v1, v2, v3]).unwrap();

        assert_eq!(s.skeleton(-2).len(), 0);
        assert_eq!(s.skeleton(-1).len(), 0);
        assert_eq!(s.skeleton(0).len(), 1);
        assert_eq!(s.vertices().len(), 3);
        assert_eq!(s.edges().len(), 3);
        assert_eq!(s.faces().len(), 1);
        assert_eq!(s.skeleton(4).len(), 0);
    }
}
