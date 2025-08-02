use super::*;
use crate::{ambient_space::AffineSpace, coordinates::Vector, oriented_simplex::OrientedSimplex};
use itertools::Itertools;

#[derive(Clone)]
pub struct Simplex<'f, FS: OrderedRingSignature + FieldSignature> {
    ambient_space: AffineSpace<'f, FS>,
    // points must be ordered w.r.t the ordering on vectors
    // points must be non-degenerate
    // points must belong to the ambient_space
    points: Vec<Vector<'f, FS>>,
}

impl<'f, FS: OrderedRingSignature + FieldSignature> std::fmt::Debug for Simplex<'f, FS> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Simplex")
            .field("points", &self.points)
            .finish()
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> PartialEq for Simplex<'f, FS> {
    fn eq(&self, other: &Self) -> bool {
        self.ambient_space == other.ambient_space && self.points == other.points
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Eq for Simplex<'f, FS> {}

impl<'f, FS: OrderedRingSignature + FieldSignature> Hash for Simplex<'f, FS>
where
    FS::Set: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.points.hash(state);
    }
}

impl<'f, FS: OrderedRingSignature + FieldSignature> Simplex<'f, FS>
where
    AffineSpace<'f, FS>: Clone,
{
    pub fn new(
        ambient_space: AffineSpace<'f, FS>,
        mut points: Vec<Vector<'f, FS>>,
    ) -> Result<Self, &'static str> {
        for point in &points {
            assert_eq!(ambient_space.borrow(), point.ambient_space());
        }
        points.sort_unstable();
        if ambient_space
            .borrow()
            .are_points_affine_independent(points.iter().collect())
        {
            Ok(Self {
                ambient_space,
                points,
            })
        } else {
            Err("Can't make a simplex using degenerate points")
        }
    }

    pub fn ambient_space(&self) -> &AffineSpace<'f, FS> {
        &self.ambient_space
    }

    pub fn n(&self) -> usize {
        self.points.len()
    }

    pub fn points(&self) -> &Vec<Vector<'f, FS>> {
        &self.points
    }

    pub fn into_points(self) -> Vec<Vector<'f, FS>> {
        self.points
    }

    pub fn point(&self, i: usize) -> &Vector<'f, FS> {
        &self.points[i]
    }

    pub fn skeleton(&self, skel_n: isize) -> Vec<Self> {
        if skel_n < 0 {
            vec![]
        } else {
            #[allow(clippy::cast_sign_loss)]
            let skel_n = skel_n as usize;
            let mut parts = vec![];
            for subset in (0..self.points.len()).combinations(skel_n) {
                let part = Self::new(
                    self.ambient_space.clone(),
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
        #[allow(clippy::cast_possible_wrap)]
        self.skeleton(self.points.len() as isize - 2)
    }

    pub fn facets(&self) -> Vec<Self> {
        #[allow(clippy::cast_possible_wrap)]
        self.skeleton(self.points.len() as isize - 1)
    }

    pub fn facet(&self, k: usize) -> Self {
        assert!(k <= self.points.len());
        Self::new(self.ambient_space.clone(), {
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
            .map(|sub_points| Self::new(self.ambient_space.clone(), sub_points).unwrap())
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

    #[allow(clippy::needless_pass_by_value)]
    pub fn sub_simplex(&self, pts: Vec<usize>) -> Self {
        Self::new(
            self.ambient_space().clone(),
            pts.iter().map(|idx| self.points[*idx].clone()).collect(),
        )
        .unwrap()
    }

    pub fn oriented_facet(&self, k: usize) -> OrientedSimplex<'f, FS> {
        //return the oriented facet of self with negative side on the outside and positive side on the inside
        assert!(k <= self.points.len());
        let mut facet_points = self.points.clone();
        let other_pt = facet_points.remove(k);
        OrientedSimplex::new_with_positive_point(
            self.ambient_space.clone(),
            facet_points,
            &other_pt,
        )
        .unwrap()
    }

    pub fn oriented_facets(&self) -> Vec<OrientedSimplex<'f, FS>> {
        assert_eq!(self.ambient_space.borrow().affine_dimension(), self.n());
        (0..self.n()).map(|k| self.oriented_facet(k)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;

    #[test]
    fn make_simplex() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(1)]);
        let v2 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(1)]);
        let s = Simplex::new(space.clone(), vec![v1, v2, v3]);
        assert!(s.is_ok());

        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(0)]);
        let v2 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(space.clone(), vec![Rational::from(2), Rational::from(0)]);
        let s = Simplex::new(space.clone(), vec![v1, v2, v3]);
        assert!(s.is_err());
    }

    #[test]
    fn simplex_skeleton() {
        let space = AffineSpace::new_linear(Rational::structure_ref(), 2);
        let v1 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(1)]);
        let v2 = Vector::new(space.clone(), vec![Rational::from(1), Rational::from(0)]);
        let v3 = Vector::new(space.clone(), vec![Rational::from(0), Rational::from(1)]);
        let s = Simplex::new(space.clone(), vec![v1, v2, v3]).unwrap();

        assert_eq!(s.skeleton(-2).len(), 0);
        assert_eq!(s.skeleton(-1).len(), 0);
        assert_eq!(s.skeleton(0).len(), 1);
        assert_eq!(s.vertices().len(), 3);
        assert_eq!(s.edges().len(), 3);
        assert_eq!(s.faces().len(), 1);
        assert_eq!(s.skeleton(4).len(), 0);
    }
}
