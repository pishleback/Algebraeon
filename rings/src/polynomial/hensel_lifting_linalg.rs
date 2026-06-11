use crate::{
    matrix::{Matrix, RingMatricesSignature},
    polynomial::{Polynomial, PolynomialStructure, ToPolynomialSignature},
    structure::{
        AdditionSignature, AdditiveGroupSignature, CancellativeMultiplicationSignature,
        EuclideanDomainSignature, FactoringMonoidSignature, FactoringStructure, FieldSignature,
        GreatestCommonDivisorSignature, MultiplicationSignature, MultiplicativeMonoidSignature,
        NonZeroFactored, QuotientRingGetPrincipalIdealSignature, RingToQuotientRingSignature,
        SemiModuleSignature, ZeroEqSignature,
    },
};
use algebraeon_structures::*;
use std::sync::Arc;

fn compute_lifting_mat_mod<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingMod: QuotientRingGetPrincipalIdealSignature<Ring> + EqSignature,
>(
    ring_mod_p: &RingMod,
    fs: &[Polynomial<Ring::Elem>],
) -> Matrix<RingMod::Elem> {
    let ring = ring_mod_p.pre_quotient_set();

    let fs_count = fs.len();
    let fs_deg = fs
        .iter()
        .map(|f| {
            ring_mod_p
                .pre_quotient_set()
                .polynomials()
                .degree(f)
                .unwrap()
        })
        .collect::<Vec<_>>();
    let h_deg = fs_deg.iter().sum::<usize>();

    let fs_prod = ring.polynomials().product(fs);
    let fs_prod_excluding_each = (0..fs_count)
        .map(|i| ring.polynomials().try_divide(&fs_prod, &fs[i]).unwrap())
        .collect::<Vec<_>>();
    let fs_prod_excluding_each_mod_p = fs_prod_excluding_each
        .iter()
        .map(|f| f.apply_map(|c| ring_mod_p.project_ref(c)))
        .collect::<Vec<_>>();

    let mat = Matrix::<RingMod::Elem>::from_cols({
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
                        ring_mod_p
                            .polynomials()
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
    mat
}

// A factorization h = h.lc() * f1 * f2 * ... * fn mod p^k
#[derive(Clone)]
pub struct HenselFactorization<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingMod: QuotientRingGetPrincipalIdealSignature<Ring> + EqSignature,
> {
    make_ring_mod: Arc<dyn Fn(&Ring::Elem) -> RingMod + Send + Sync>,
    polys: PolynomialStructure<Ring, Ring>,
    p: Ring::Elem, // An irreducible element of Ring
    k: Natural,
    h: Polynomial<Ring::Elem>, // A polynomial over Ring
    h_deg: usize,
    fs_count: usize,
    fs: Vec<Polynomial<Ring::Elem>>, // monic Polynomials over Ring mod p^k
    fs_deg: Vec<usize>,
    ring_mod_p: RingMod,
    // p^k
    pk: Ring::Elem,
    // The inverse of the leading coefficient of h modulo p
    lc_inv_mod_p: RingMod::Elem,
    // A matrix used to compute lifts. We cache it here since it is unchanged during lifts.
    lifting_mat_mod_p: Matrix<RingMod::Elem>,
    lifting_mat_inv_mod_p: Matrix<RingMod::Elem>,
}

impl<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingMod: QuotientRingGetPrincipalIdealSignature<Ring> + EqSignature,
> std::fmt::Debug for HenselFactorization<Ring, RingMod>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HenselFactorization")
            .field("polys", &self.polys)
            .field("p", &self.p)
            .field("k", &self.k)
            .field("h", &self.h)
            .field("h_deg", &self.h_deg)
            .field("fs_count", &self.fs_count)
            .field("fs", &self.fs)
            .field("fs_deg", &self.fs_deg)
            .field("ring_mod_p", &self.ring_mod_p)
            .field("pk", &self.pk)
            .field("lc_inv_mod_p", &self.lc_inv_mod_p)
            .field("lifting_mat_mod_p", &self.lifting_mat_mod_p)
            .field("lifting_mat_inv_mod_p", &self.lifting_mat_inv_mod_p)
            .finish()
    }
}

impl<
    Ring: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
    RingMod: QuotientRingGetPrincipalIdealSignature<Ring> + EqSignature,
> HenselFactorization<Ring, RingMod>
{
    pub fn ring(&self) -> &Ring {
        self.ring_mod_p.pre_quotient_set()
    }

    pub fn base_modulus(&self) -> &Ring::Elem {
        &self.p
    }

    pub fn modulus(&self) -> &Ring::Elem {
        &self.pk
    }

    fn lc(&self) -> &Ring::Elem {
        self.ring().polynomials().leading_coeff(&self.h).unwrap()
    }

    fn compute_err_poly_mod_p_poly(
        ring: &Ring,
        ring_mod_p: &RingMod,
        p: &Ring::Elem,
        k: &Natural,
        h: &Polynomial<Ring::Elem>,
        fs: &[Polynomial<Ring::Elem>],
        lc: &Ring::Elem,
    ) -> Polynomial<RingMod::Elem> {
        let pk0 = ring.nat_pow(p, k);
        let pk1 = ring.mul(&pk0, p);

        // polynomials mod p^(k+1)
        let polys_mod_pk1 = ring
            .euclidean_quotient_ring(pk1.clone())
            .unwrap()
            .into_polynomials();

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

        assert_eq!(self.lifting_mat_inv_mod_p.rows(), self.h_deg);
        assert_eq!(self.lifting_mat_inv_mod_p.cols(), self.h_deg);
        assert!(
            self.ring_mod_p.matrix_structure().equal(
                &self
                    .ring_mod_p
                    .matrix_structure()
                    .mul(&self.lifting_mat_mod_p, &self.lifting_mat_inv_mod_p)
                    .unwrap(),
                &self.ring_mod_p.matrix_structure().ident(self.h_deg)
            )
        );
        let mat = compute_lifting_mat_mod(&self.ring_mod_p, &self.fs);
        assert!(
            self.ring_mod_p
                .matrix_structure()
                .equal(&self.lifting_mat_mod_p, &mat)
        );

        if !self.ring().equal(&self.p, &self.ring_mod_p.modulus()) {
            return Err("Incorrect self.ring_mod_p".to_string());
        }
        if self.k == Natural::ZERO {
            return Err("k=0 is not valid".to_string());
        }

        if !self.ring_mod_p.equal(
            &self
                .ring_mod_p
                .mul(&self.ring_mod_p.project_ref(self.lc()), &self.lc_inv_mod_p),
            &self.ring_mod_p.one(),
        ) {
            return Err("leading coeff bad inverse modulo p".to_string());
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
            .euclidean_quotient_ring(self.pk.clone())
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

        if !self
            .ring()
            .equal(&self.ring().nat_pow(&self.p, &self.k), &self.pk)
        {
            return Err("Incorrect self.pk".to_string());
        }

        Ok(())
    }

    pub fn new_unchecked(
        make_ring_mod: impl Fn(&Ring::Elem) -> RingMod + Send + Sync + 'static,
        ring_mod_p: RingMod,
        k: Natural,
        h: Polynomial<Ring::Elem>,
        fs: Vec<Polynomial<Ring::Elem>>,
        lifting_mat_inv_mod_p: Matrix<RingMod::Elem>,
        lc_inv_mod_p: RingMod::Elem,
    ) -> Self {
        let p = ring_mod_p.modulus().into_owned();
        let ring = ring_mod_p.pre_quotient_set();
        let polys = ring.clone().into_polynomials();
        let pk = ring.nat_pow(&p, &k);
        let h_deg = polys.degree(&h).unwrap();
        let fs_count = fs.len();
        let fs_deg = fs
            .iter()
            .map(|f| polys.degree(f).unwrap())
            .collect::<Vec<_>>();
        let lifting_mat_mod_p = compute_lifting_mat_mod(&ring_mod_p, &fs);

        let s = Self {
            make_ring_mod: Arc::new(make_ring_mod),
            ring_mod_p,
            polys,
            p,
            k,
            pk,
            h,
            h_deg,
            fs_count,
            fs,
            fs_deg,
            lifting_mat_mod_p,
            lifting_mat_inv_mod_p,
            lc_inv_mod_p,
        };
        #[cfg(debug_assertions)]
        s.check().unwrap();
        s
    }

    pub fn new_unchecked_mod_field<
        FieldModP: QuotientRingGetPrincipalIdealSignature<Ring> + FieldSignature,
    >(
        make_ring_mod: impl Fn(&Ring::Elem) -> RingMod + Send + Sync + 'static,
        field_mod_p: FieldModP,
        k: Natural,
        h: Polynomial<Ring::Elem>,
        fs: Vec<Polynomial<Ring::Elem>>,
    ) -> Self {
        let ring = field_mod_p.pre_quotient_set();
        let ring_mod_p = make_ring_mod(field_mod_p.modulus().as_ref());
        let lifting_mat_inv_mod_p = field_mod_p
            .matrix_structure()
            .inv(compute_lifting_mat_mod(&field_mod_p, &fs))
            .unwrap()
            .apply_map_into(|x| ring_mod_p.project(field_mod_p.unproject(x)));
        let lc = ring.polynomials().leading_coeff(&h).unwrap();
        let lc_inv_mod_p = ring_mod_p.project(
            field_mod_p.unproject(
                field_mod_p
                    .try_reciprocal(&field_mod_p.project_ref(lc))
                    .unwrap(),
            ),
        );
        Self::new_unchecked(
            make_ring_mod,
            ring_mod_p,
            k,
            h,
            fs,
            lifting_mat_inv_mod_p,
            lc_inv_mod_p,
        )
    }

    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn from_mod_field_factorization<
        FieldModP: QuotientRingGetPrincipalIdealSignature<Ring> + FieldSignature,
        RingModPB: BorrowedStructure<FieldModP>,
        RingModPPolyB: BorrowedStructure<PolynomialStructure<FieldModP, RingModPB>>,
        NatB: BorrowedStructure<NaturalCanonicalStructure>,
    >(
        make_ring_mod: impl Fn(&Ring::Elem) -> RingMod + Send + Sync + 'static,
        fs_structure: &FactoringStructure<
            PolynomialStructure<FieldModP, RingModPB>,
            RingModPPolyB,
            NaturalCanonicalStructure,
            NatB,
        >,
        fs: NonZeroFactored<Polynomial<FieldModP::Elem>, Natural>,
        h: Polynomial<Ring::Elem>,
    ) -> Option<Self> {
        let polys_mod_p = fs_structure.objects().clone();
        let field_mod_p = polys_mod_p.coeff_ring().clone();
        let mut fs_vec = vec![];
        let (_unit, factors) = fs.into_unit_and_powers();
        for (factor, power) in factors {
            if power == Natural::ONE {
                fs_vec.push(factor.apply_map_into(|c| field_mod_p.unproject(c)));
            } else {
                return None;
            }
        }
        let hensel_factorization =
            Self::new_unchecked_mod_field(make_ring_mod, field_mod_p, Natural::ONE, h, fs_vec);
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        Some(hensel_factorization)
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

        // It might be possible to update some data each lift so that this is much cheaper to compute and update per lift
        let err_poly_mod_p = Self::compute_err_poly_mod_p_poly(
            self.ring(),
            &self.ring_mod_p,
            &self.p,
            &self.k,
            &self.h,
            &self.fs,
            self.lc(),
        );
        let t_poly = &err_poly_mod_p;
        let t_coeffs = (0..self.h_deg)
            .map(|i| {
                self.ring_mod_p
                    .polynomials()
                    .coeff(t_poly, i)
                    .as_ref()
                    .clone()
            })
            .collect::<Vec<_>>();

        // solve the linear system for g1, g2, g3
        let delta_sol = self
            .ring_mod_p
            .matrix_structure()
            .apply_col(&self.lifting_mat_inv_mod_p, &t_coeffs);

        // extract the coefficients of g1, g2, g3 from the single vector solution to the linear system and add them to self.fs to lift the factorization from mod p^k to mod p^(k+1)
        let mut offset = 0;

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
                .map(|d| self.ring_mod_p.unproject_ref(d))
                .collect::<Vec<_>>();

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

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    pub fn quadratic_lift(&mut self) {
        assert_eq!(self.k, Natural::ONE);
        /*
        Start with
            h = lc * f1 * f2 * f3 mod p
        and aim for
            h = lc * f1 * f2 * f3 mod q
        where q = p^2  and lift
            lc_inv_mod_p -> lc_inv_mod_q
            lifting_mat_inv_mod_p -> lifting_mat_inv_mod_q
         */

        self.linear_lift();
        debug_assert_eq!(self.k, Natural::TWO);

        let q = self.pk.clone();
        let ring_mod_q = (self.make_ring_mod)(&q);

        /*
        For lifting the lifting matrix inverse:
            `mat_mod_q * mat_inv_mod_p = I` mod p
            `mat_mod_q * mat_inv_mod_p = I + p * error` mod q for some `error` mod p
        want `correction` mod p such that
            `mat_mod_q * correction = -error` mod p so that
        i.e.
            `correction := -mat_mod_p_inv * error` mod p
        so that
            `mat_mod_q * (mat_inv_mod_p + p * correction) = I + p * error - p * error = I` mod q
         */
        let lifting_mat_mod_q = compute_lifting_mat_mod(&ring_mod_q, &self.fs);
        let lifting_mat_inv_mod_p = self
            .lifting_mat_inv_mod_p
            .apply_map(|x| ring_mod_q.project(self.ring_mod_p.unproject_ref(x)));

        let error_times_p = ring_mod_q
            .matrix_structure()
            .add(
                &ring_mod_q
                    .matrix_structure()
                    .mul(&lifting_mat_mod_q, &lifting_mat_inv_mod_p)
                    .unwrap(),
                &ring_mod_q
                    .matrix_structure()
                    .neg(ring_mod_q.matrix_structure().ident(self.h_deg)),
            )
            .unwrap();

        debug_assert!(
            self.ring_mod_p.matrix_structure().equal(
                &error_times_p.apply_map(|x| self.ring_mod_p.project(ring_mod_q.unproject_ref(x))),
                &self
                    .ring_mod_p
                    .matrix_structure()
                    .zero(self.h_deg, self.h_deg)
            )
        ); // error = 0 mod p
        let error_mod_p = error_times_p.apply_map(|x| {
            self.ring_mod_p.project(
                self.ring()
                    .try_divide(&ring_mod_q.unproject_ref(x), &self.p)
                    .unwrap(),
            )
        });
        debug_assert!(
            ring_mod_q.matrix_structure().equal(
                &ring_mod_q
                    .matrix_structure()
                    .mul(&lifting_mat_mod_q, &lifting_mat_inv_mod_p)
                    .unwrap(),
                &ring_mod_q
                    .matrix_structure()
                    .add(
                        &ring_mod_q.matrix_structure().ident(self.h_deg),
                        &error_times_p
                    )
                    .unwrap()
            )
        );
        let correction_mod_p = self.ring_mod_p.matrix_structure().neg(
            self.ring_mod_p
                .matrix_structure()
                .mul(&self.lifting_mat_inv_mod_p, &error_mod_p)
                .unwrap(),
        );
        let lifting_mat_inv_mod_q = ring_mod_q
            .matrix_structure()
            .add(
                &lifting_mat_inv_mod_p,
                &correction_mod_p.apply_map_into(|x| {
                    ring_mod_q.project(self.ring().mul(&self.p, &self.ring_mod_p.unproject(x)))
                }),
            )
            .unwrap();

        /*
        For lifting the lc inverse:
            `lc * lc_inv_mod_p = 1` mod p
            `lc * lc_inv_mod_p = 1 + p * error` mod q for some `error` mod p
        want `correction` mod p such that
            `lc * correction = -error` mod p so that
        i.e.
            `correction := -lc_inv_mod_p * error` mod p
        so that
            `lc * (lc_inv_mod_p + p * correction) = 1 + p * error - p * error = 1` mod q
         */
        let lc = self.lc();
        let lc_mod_q = ring_mod_q.project_ref(lc);
        let lc_inv_mod_p = ring_mod_q.project(self.ring_mod_p.unproject_ref(&self.lc_inv_mod_p));
        let error_times_p =
            ring_mod_q.sub(&ring_mod_q.mul(&lc_mod_q, &lc_inv_mod_p), &ring_mod_q.one());
        debug_assert!(
            self.ring_mod_p.is_zero(
                &self
                    .ring_mod_p
                    .project(ring_mod_q.unproject_ref(&error_times_p))
            )
        ); // error = 0 mod p
        let error_mod_p = self.ring_mod_p.project(
            self.ring()
                .try_divide(&ring_mod_q.unproject(error_times_p), &self.p)
                .unwrap(),
        );
        let correction_mod_p = self
            .ring_mod_p
            .neg(&self.ring_mod_p.mul(&self.lc_inv_mod_p, &error_mod_p));
        let lc_inv_mod_q = ring_mod_q.add(
            &lc_inv_mod_p,
            &ring_mod_q.project(
                self.ring()
                    .mul(&self.p, &self.ring_mod_p.unproject(correction_mod_p)),
            ),
        );

        self.pk = q.clone();
        self.p = q;
        self.k = Natural::ONE;
        self.ring_mod_p = ring_mod_q;
        self.lifting_mat_mod_p = lifting_mat_mod_q;
        self.lifting_mat_inv_mod_p = lifting_mat_inv_mod_q;
        self.lc_inv_mod_p = lc_inv_mod_q;

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    pub fn factors(&self) -> &Vec<Polynomial<Ring::Elem>> {
        &self.fs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::{
        MetaEuclideanDivisionSignature, MetaMultiplicationSignature,
        MetaMultiplicativeMonoidSignature, RingToQuotientFieldSignature,
    };

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
        let mod5 = Integer::structure()
            .into_quotient_field(Integer::from(5))
            .unwrap();

        let fs = vec![f1, f2, f3];
        let mut hensel_fact = HenselFactorization::new_unchecked_mod_field(
            |x| {
                Integer::structure()
                    .into_euclidean_quotient_ring(x.clone())
                    .unwrap()
            },
            mod5,
            Natural::from(1u8),
            h.clone(),
            fs,
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
            .apply_map(|c| Integer::rem(c, hensel_fact.modulus()));
            println!("{:?} {:?}", lifted_product, h);
            assert_eq!(lifted_product, h);
        }
    }
}
