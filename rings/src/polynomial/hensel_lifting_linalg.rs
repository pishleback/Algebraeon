use crate::{
    matrix::{Matrix, MatrixStructure},
    polynomial::{Polynomial, PolynomialStructure, ToPolynomialSignature},
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeMultiplicationSignature, EuclideanDomainSignature, FactoringMonoidSignature,
        FactoringStructure, FieldSignature, GreatestCommonDivisorSignature,
        MultiplicationSignature, MultiplicativeMonoidSignature, NonZeroFactored,
        QuotientRingGetPrincipalIdealSignature, RingToQuotientFieldSignature,
        RingToQuotientRingSignature, SemiModuleSignature, TryReciprocalSignature, ZeroSignature,
    },
};
use algebraeon_nzq::{Natural, NaturalCanonicalStructure};
use algebraeon_sets::structure::{BorrowedStructure, EqSignature};

#[derive(Debug, Clone)]
pub struct HenselFactorization<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingModP: QuotientRingGetPrincipalIdealSignature<Ring> + FieldSignature,
> {
    // A factorization h = h.lc() * f1 * f2 * ... * fn mod p^k
    ring_mod_p: RingModP,
    polys: PolynomialStructure<Ring, Ring>,
    polys_mod_p: PolynomialStructure<RingModP, RingModP>,
    mats_mod_p: MatrixStructure<RingModP, RingModP>,
    p: Ring::Elem, // An irreducible element of Ring
    k: Natural,
    pk: Ring::Elem,            // p^k
    h: Polynomial<Ring::Elem>, // A polynomial over Ring
    h_deg: usize,
    fs_count: usize,
    fs: Vec<Polynomial<Ring::Elem>>, // monic Polynomials over Ring mod p^k
    fs_deg: Vec<usize>,
    t: Polynomial<RingModP::Elem>, // (h - h.lc() * f1 * f2 * ... * fn) / p^k mod p
    // A matrix used to compute lifts. We cache it here since it is unchanged during lifts.
    lifting_mat_inv_mod_p: Matrix<RingModP::Elem>,
    // The inverse of the leading coefficient of h modulo p
    lc_inv_mod_p: RingModP::Elem,
}

impl<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingModP: QuotientRingGetPrincipalIdealSignature<Ring> + FieldSignature,
> HenselFactorization<Ring, RingModP>
{
    pub fn ring(&self) -> &Ring {
        self.ring_mod_p.pre_quotient_set()
    }

    pub fn ring_mod_p(&self) -> &RingModP {
        &self.ring_mod_p
    }

    pub fn modulus(&self) -> Ring::Elem {
        self.ring().nat_pow(&self.p, &self.k)
    }

    fn lc(&self) -> &Ring::Elem {
        self.ring().polynomials().leading_coeff(&self.h).unwrap()
    }

    fn compute_lifting_mat_inv_mod_p(
        fs_count: usize,
        fs: &[Polynomial<Ring::Elem>],
        fs_deg: &[usize],
        h_deg: usize,
        ring_mod_p: &RingModP,
        polys_mod_p: &PolynomialStructure<RingModP, RingModP>,
        mats_mod_p: &MatrixStructure<RingModP, RingModP>,
    ) -> Matrix<RingModP::Elem> {
        debug_assert_eq!(fs_count, fs.len());
        debug_assert_eq!(fs_count, fs_deg.len());
        debug_assert_eq!(fs_deg.iter().cloned().sum::<usize>(), h_deg);

        println!("make prods");

        let fs_mod_p = fs
            .iter()
            .map(|f| f.apply_map(|c| ring_mod_p.project_ref(c)))
            .collect::<Vec<_>>();
        let fs_prod = polys_mod_p.product(&fs_mod_p);
        let fs_prod_excluding_each_mod_p = (0..fs_count)
            .map(|i| polys_mod_p.try_divide(&fs_prod, &fs_mod_p[i]).unwrap())
            .collect::<Vec<_>>();

        println!("make {}x{} mat", h_deg, h_deg);

        let mat = Matrix::<RingModP::Elem>::from_cols({
            let mut cols = vec![];
            for idx in 0..fs_count {
                let d = fs_deg[idx];
                for i in 0..d {
                    let mut col = vec![];
                    for _ in 0..i {
                        col.push(ring_mod_p.zero());
                    }
                    for j in 0..(h_deg - d + 1) {
                        col.push(
                            polys_mod_p
                                .coeff(&fs_prod_excluding_each_mod_p[idx], j)
                                .as_ref()
                                .clone(),
                        );
                    }
                    for _ in 0..(d - i - 1) {
                        col.push(ring_mod_p.zero());
                    }
                    debug_assert_eq!(col.len(), h_deg);
                    cols.push(col);
                }
            }
            debug_assert_eq!(cols.len(), h_deg);
            cols
        });
        debug_assert_eq!(mat.rows(), h_deg);
        debug_assert_eq!(mat.cols(), h_deg);
        println!("inv mat");
        mats_mod_p.inv(mat.clone()).unwrap()
    }

    fn compute_t_poly(
        ring: &Ring,
        ring_mod_p: &RingModP,
        polys_mod_p: &PolynomialStructure<RingModP, RingModP>,
        p: &Ring::Elem,
        k: &Natural,
        h: &Polynomial<Ring::Elem>,
        fs: &Vec<Polynomial<Ring::Elem>>,
        lc: &Ring::Elem,
        lc_inv_mod_p: &RingModP::Elem,
    ) -> Polynomial<RingModP::Elem> {
        let pk0 = ring.nat_pow(p, k);
        let pk1 = ring.mul(&pk0, p);

        // polynomials mod p^(k+1)
        let polys_mod_pk1 = ring.quotient_ring(pk1.clone()).unwrap().into_polynomials();

        // compute t = (h - lc * (f1 * f2 * f3)) / p^k mod p
        polys_mod_pk1
            .sub(
                h,
                &polys_mod_pk1
                    .scalar_mul(&polys_mod_pk1.product(&fs.iter().collect::<Vec<_>>()), lc),
            )
            .apply_map_into(|c| ring_mod_p.project(ring.try_divide(&c, &pk0).unwrap()))
    }

    pub fn check(&self) -> Result<(), String> {
        assert_eq!(self.h_deg, self.polys.degree(&self.h).unwrap());
        assert_eq!(self.fs_count, self.fs.len());
        assert_eq!(self.fs_count, self.fs_deg.len());
        for i in 0..self.fs_count {
            assert_eq!(self.fs_deg[i], self.polys.degree(&self.fs[i]).unwrap());
        }

        assert!(self.mats_mod_p.equal(
            &self.lifting_mat_inv_mod_p,
            &Self::compute_lifting_mat_inv_mod_p(
                self.fs_count,
                &self.fs,
                &self.fs_deg,
                self.h_deg,
                &self.ring_mod_p,
                &self.polys_mod_p,
                &self.mats_mod_p,
            )
        ));

        if !self.ring().is_irreducible(&self.p) {
            return Err(format!("p={:?} is not irreducible", self.p));
        }
        if self.k == Natural::ZERO {
            return Err("k=0 is not valid".to_string());
        }

        let ring_mod_p = self.ring().quotient_field(self.p.clone()).unwrap();
        if !ring_mod_p.is_unit(self.lc()) {
            return Err("leading coeff is not invertible modulo p".to_string());
        }

        let ring_poly = self.ring().polynomials();
        for f in &self.fs {
            if !ring_poly.is_monic(f) {
                return Err("modular factor is not monic".to_string());
            }
        }
        for i in 0..self.fs.len() {
            for j in 0..i {
                let fi = &self.fs[i];
                let fj = &self.fs[j];
                if ring_poly
                    .degree(&ring_poly.gcd_by_primitive_subresultant(fi.clone(), fj.clone()))
                    != Some(0)
                {
                    return Err("modular factors are not pairwise coprime".to_string());
                }
            }
        }

        let ring_poly_mod_pk = self
            .ring()
            .quotient_ring(self.modulus())
            .unwrap()
            .into_polynomials();
        if !ring_poly_mod_pk.equal(
            &self.h,
            &ring_poly_mod_pk.mul(
                &Polynomial::constant(self.lc().clone()),
                &ring_poly_mod_pk.product(&self.fs.iter().collect::<Vec<_>>()),
            ),
        ) {
            return Err("factorization mod p^k is not correct".to_string());
        }

        let mut deg_sum = 0;
        for f in &self.fs {
            deg_sum += ring_poly.degree(f).unwrap();
        }
        debug_assert_eq!(ring_poly.degree(&self.h).unwrap(), deg_sum);

        if !self.polys_mod_p.equal(
            &self.t,
            &Self::compute_t_poly(
                self.ring(),
                &self.ring_mod_p,
                &self.polys_mod_p,
                &self.p,
                &self.k,
                &self.h,
                &self.fs,
                self.lc(),
                &self.lc_inv_mod_p,
            ),
        ) {
            return Err("Incorrect self.t".to_string());
        }

        if !self
            .ring()
            .equal(&self.ring().nat_pow(&self.p, &self.k), &self.pk)
        {
            return Err("Incorrect self.pk".to_string());
        }

        Ok(())
    }

    pub fn new_unchecked(
        ring_mod_p: RingModP,
        k: Natural,
        h: Polynomial<Ring::Elem>,
        fs: Vec<Polynomial<Ring::Elem>>,
    ) -> Self {
        let p = ring_mod_p.modulus().into_owned();
        let ring = ring_mod_p.pre_quotient_set();
        let polys = ring.clone().into_polynomials();
        let polys_mod_p = ring_mod_p.clone().into_polynomials();
        let mats_mod_p = MatrixStructure::new(ring_mod_p.clone());

        let pk = ring.nat_pow(&p, &k);

        let h_deg = polys.degree(&h).unwrap();
        let fs_count = fs.len();
        let fs_deg = fs
            .iter()
            .map(|f| polys.degree(f).unwrap())
            .collect::<Vec<_>>();

        let lifting_mat_inv_mod_p = Self::compute_lifting_mat_inv_mod_p(
            fs_count,
            &fs,
            &fs_deg,
            h_deg,
            &ring_mod_p,
            &polys_mod_p,
            &mats_mod_p,
        );

        let lc = polys.leading_coeff(&h).unwrap();
        let lc_inv_mod_p = ring_mod_p
            .try_reciprocal(&ring_mod_p.project_ref(lc))
            .unwrap();

        let t = Self::compute_t_poly(
            ring,
            &ring_mod_p,
            &polys_mod_p,
            &p,
            &k,
            &h,
            &fs,
            lc,
            &lc_inv_mod_p,
        );

        let s = Self {
            ring_mod_p,
            polys,
            polys_mod_p,
            mats_mod_p,
            p,
            k,
            pk,
            h,
            h_deg,
            fs_count,
            fs,
            fs_deg,
            lifting_mat_inv_mod_p,
            lc_inv_mod_p,
            t,
        };
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }

    // lift the factorization from mod p^k to mod p^(k+1)
    pub fn linear_lift(&mut self) {
        /*
        We start with
           h = lc * f1 * f2 * f3 mod p^k
        and search for g1, g2, g3 mod p such that
           h = lc * (f1 + p^k * g1) * (f2 + p^k * g2) * (f3 + p^k * g3) mod p^(k+1)
        Do some rearranging...
             lc * (f1 + p^k * g1) * (f2 + p^k * g2) * (f3 + p^k * g3) = h mod p^(k+1)
         iff lc * (f1 * f2 * f3) + lc * p^k * (g1 * f2 * f3 + f1 * g2 * g3 + f1 * f2 * g3) + p^(2k) * (stuff) = h mod p^(k+1)
         iff lc * p^k * (g1 * f2 * f3 + f1 * g2 * g3 + f1 * f2 * g3) = h - lc * (f1 * f2 * f3) mod p^(k+1)
        at this point both sides are 0 mod p^k, so we can divide through by p^k to get
         iff lc * (g1 * f2 * f3 + f1 * g2 * g3 + f1 * f2 * g3) = (h - lc * (f1 * f2 * f3)) / p^k mod p
        one of our assumptions is that lc is invertible mod p, so
         iff g1 * f2 * f3 + f1 * g2 * g3 + f1 * f2 * g3 = (h - lc * (f1 * f2 * f3)) / (lc * p^k) mod p
        It's a theorem (Hensel lifting) that there exist unique g1, g2, g3 satisfying this
        with deg(g1) < deg(f1), deg(g2) < deg(f2), deg(g3) < deg(f3).
        We can solve for g1, g2, g3 using a system of linear equations mod p:
        If e.g.
            deg(f1)=2, deg(f2)=3, deg(f3)=2
            f12 = f1 * f2 = f12_0 + f12_1 * λ + f12_2 * λ^2 + f12_3 * λ^3 + f12_4 * λ^4 + f12_5 * λ^5
            f13 = f1 * f3 = f13_0 + f13_1 * λ + f13_2 * λ^2 + f13_3 * λ^3 + f13_4 * λ^4
            f23 = f2 * f3 = f23_0 + f23_1 * λ + f23_2 * λ^2 + f23_3 * λ^3 + f23_4 * λ^4 + f23_5 * λ^5
            (h - lc * (f1 * f2 * f3)) / (lc * p^k) = t0 + t1 * λ + t2 * λ^2 + t3 * λ^3 + t4 * λ^4 + t5 * λ^5 + t6 * λ^6
        and we are solving for
            g1 = x0 + x1 * λ
            g2 = y0 + y1 * λ + y2 * λ^2
            g3 = z0 + z0 * λ
        then the linear system equivalent to
            g1 * f2 * f3 + f1 * g2 * g3 + f1 * f2 * g3 = (h - lc * (f1 * f2 * f3)) / (lc * p^k) mod p
        looks like

        / f23_0   0   f13_0   0     0   f12_0   0   \ / x0 \   / t0 \
        | f23_1 f23_0 f13_1 f13_0   0   f12_1 f12_0 | | x1 |   | t1 |
        | f23_2 f23_1 f13_2 f13_1 f13_0 f12_2 f12_1 | | y0 |   | t2 |
        | f23_3 f23_2 f13_3 f13_2 f13_1 f12_3 f12_2 | | y1 | = | t3 |
        | f23_4 f23_3 f13_4 f13_3 f13_2 f12_4 f12_3 | | y2 |   | t4 |
        | f23_5 f23_4   0   f13_4 f13_3 f12_5 f12_4 | | z0 |   | t5 |
        \   0   f23_5   0     0   f13_4   0   f12_5 / \ z1 /   \ t6 /

        We've pre-computed the inverse of the matrix on the left. It does not change between lifts since f1, f2, f3 do not change mod p between lifts.
        */

        #[cfg(debug_assertions)]
        self.check().unwrap();

        let t_poly = &self.t;
        let t_coeffs = (0..self.h_deg)
            .map(|i| self.polys_mod_p.coeff(t_poly, i).as_ref().clone())
            .collect::<Vec<_>>();

        // solve the linear system for g1, g2, g3
        let delta_sol = self
            .mats_mod_p
            .apply_col(&self.lifting_mat_inv_mod_p, &t_coeffs);

        // extract the coefficients of g1, g2, g3 from the single vector solution to the linear system and add them to self.fs to lift the factorization from mod p^k to mod p^(k+1)
        let mut offset = 0;

        let mut sum_t_thing = self.polys_mod_p.zero();
        let fs_mod_p = self
            .fs
            .iter()
            .map(|f| f.apply_map(|c| self.ring_mod_p.project_ref(c)))
            .collect::<Vec<_>>();
        let fs_prod = self.polys_mod_p.product(&fs_mod_p);
        let fs_prod_excluding_each_mod_p = (0..self.fs_count)
            .map(|i| self.polys_mod_p.try_divide(&fs_prod, &fs_mod_p[i]).unwrap())
            .collect::<Vec<_>>();

        for i in 0..self.fs.len() {
            let f_deg = self.fs_deg[i];
            let mut delta_mod_p = vec![];
            for j in 0..f_deg {
                delta_mod_p.push(
                    self.ring_mod_p
                        .mul(&delta_sol[offset + j], &self.lc_inv_mod_p),
                );
            }
            let delta = delta_mod_p
                .iter()
                .map(|d| {
                    self.ring_mod_p
                        .unproject_ref(d)
                })
                .collect::<Vec<_>>();

            self.polys_mod_p.add_mut(
                &mut sum_t_thing,
                &self.polys_mod_p.mul(
                    &Polynomial::from_coeffs(delta_mod_p),
                    &fs_prod_excluding_each_mod_p[i],
                ),
            );

            self.polys.add_mut(
                self.fs.get_mut(i).unwrap(),
                &self
                    .polys
                    .scalar_mul(&Polynomial::from_coeffs(delta), &self.pk),
            );
            offset += f_deg;
        }

        self.k += Natural::ONE;
        self.pk = self.ring().mul(&self.pk, &self.p);

        let new_t = self.polys_mod_p.sub(
            &self.t,
            &self
                .polys_mod_p
                .scalar_mul(&sum_t_thing, &self.ring_mod_p.project_ref(self.lc())),
        );

        self.t = Self::compute_t_poly(
            self.ring(),
            &self.ring_mod_p,
            &self.polys_mod_p,
            &self.p,
            &self.k,
            &self.h,
            &self.fs,
            self.lc(),
            &self.lc_inv_mod_p,
        );

        println!(
            "{} {:?} {:?} {:?}",
            self.k,
            self.t,
            new_t,
            self.polys_mod_p.equal(&self.t, &new_t)
        );

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    pub fn factors(&self) -> &Vec<Polynomial<Ring::Elem>> {
        &self.fs
    }

    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn from_mod_p_factorization<
        RingModPB: BorrowedStructure<RingModP>,
        RingModPPolyB: BorrowedStructure<PolynomialStructure<RingModP, RingModPB>>,
        NatB: BorrowedStructure<NaturalCanonicalStructure>,
    >(
        fs_structure: &FactoringStructure<
            PolynomialStructure<RingModP, RingModPB>,
            RingModPPolyB,
            NaturalCanonicalStructure,
            NatB,
        >,
        fs: NonZeroFactored<Polynomial<RingModP::Elem>, Natural>,
        h: Polynomial<Ring::Elem>,
    ) -> Option<Self> {
        let poly_ring_mod = fs_structure.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring().clone();
        let mut fs_vec = vec![];
        let (_unit, factors) = fs.into_unit_and_powers();
        for (factor, power) in factors {
            if power == Natural::ONE {
                fs_vec.push(factor.apply_map_into(|c| ring_mod.unproject(c)));
            } else {
                return None;
            }
        }
        let hensel_factorization = Self::new_unchecked(ring_mod, Natural::ONE, h, fs_vec);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        Some(hensel_factorization)
    }

    pub fn from_mod_p_factorization_unchecked<
        RingModPB: BorrowedStructure<RingModP>,
        RingModPPolyB: BorrowedStructure<PolynomialStructure<RingModP, RingModPB>>,
        NatB: BorrowedStructure<NaturalCanonicalStructure>,
    >(
        fs_structure: &FactoringStructure<
            PolynomialStructure<RingModP, RingModPB>,
            RingModPPolyB,
            NaturalCanonicalStructure,
            NatB,
        >,
        fs: NonZeroFactored<Polynomial<RingModP::Elem>, Natural>,
        h: Polynomial<Ring::Elem>,
    ) -> Self {
        let poly_ring_mod = fs_structure.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring().clone();
        let mut fs_vec = vec![];
        let (_unit, factors) = fs.into_unit_and_powers();
        for (factor, power) in factors {
            debug_assert_eq!(power, Natural::ONE);
            fs_vec.push(factor.apply_map_into(|c| ring_mod.unproject(c)));
        }
        let hensel_factorization = Self::new_unchecked(ring_mod, Natural::ONE, h, fs_vec);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        hensel_factorization
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{
        MetaEuclideanDivisionSignature, MetaMultiplicationSignature,
        MetaMultiplicativeMonoidSignature,
    };
    use algebraeon_nzq::*;
    use algebraeon_sets::structure::MetaType;

    #[test]
    fn test_hensel_factorization_example1_mod5() {
        let f1 =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-1), Integer::from(1)]);
        let f2 = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(3),
            Integer::from(2),
            Integer::from(1),
            Integer::from(1),
        ]);
        let f3 = Polynomial::from_coeffs(vec![
            Integer::from(4),
            Integer::from(3),
            Integer::from(2),
            Integer::from(1),
        ]);
        let h = Polynomial::product(&[&Polynomial::constant(Integer::from(8)), &f1, &f2, &f3]);
        let h = h.apply_map(|c| Integer::rem(c, &Integer::from(5)));
        //h = prod fs mod 5
        //fs are coprime
        //h and fs are monic

        //set up bezout coefficients for hensel lifting the factorization modulo 25, 125, ...
        let mut hensel_fact = HenselFactorization::new_unchecked(
            Integer::structure()
                .into_quotient_field(Integer::from(5))
                .unwrap(),
            Natural::from(1u8),
            h.clone(),
            vec![f1, f2, f3],
        );

        hensel_fact.check().unwrap();
        println!("5^1: {:?}", hensel_fact.factors());
        for i in 2..20 {
            hensel_fact.linear_lift();
            hensel_fact.check().unwrap();
            println!("5^{}: {:?}", i, hensel_fact.factors());
            let lifted_product = Polynomial::mul(
                &Polynomial::constant(Polynomial::leading_coeff(&h).unwrap().clone()),
                &Polynomial::product(hensel_fact.factors()),
            )
            .apply_map(|c| Integer::rem(c, &hensel_fact.modulus()));
            println!("{:?} {:?}", lifted_product, h);
            assert_eq!(lifted_product, h);
        }
    }
}
