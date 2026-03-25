use crate::{
    matrix::{Matrix, MatrixStructure},
    polynomial::{Polynomial, PolynomialStructure, ToPolynomialSignature},
    structure::{
        AdditionSignature, AdditiveGroupSignature, CancellativeMultiplicationSignature,
        EuclideanDomainSignature, EuclideanRemainderQuotientStructure, FactoringMonoidSignature,
        FactoringStructure, GreatestCommonDivisorSignature, MultiplicationSignature,
        MultiplicativeMonoidSignature, NonZeroFactored, QuotientRingGetPrincipalIdealSignature,
        RingToQuotientFieldSignature, RingToQuotientRingSignature, SemiModuleSignature,
        TryReciprocalSignature,
    },
};
use algebraeon_nzq::{Natural, NaturalCanonicalStructure};
use algebraeon_sets::structure::{BorrowedStructure, EqSignature};

#[derive(Debug, Clone)]
pub struct HenselFactorization<
    RS: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
> {
    // A factorization h = h.lc() * f1 * f2 * ... * fn mod p^k
    ring: RS,
    polys: PolynomialStructure<RS, RS>,
    polys_mod_p: PolynomialStructure<
        EuclideanRemainderQuotientStructure<RS, RS, true>,
        EuclideanRemainderQuotientStructure<RS, RS, true>,
    >,
    mats_mod_p: MatrixStructure<
        EuclideanRemainderQuotientStructure<RS, RS, true>,
        EuclideanRemainderQuotientStructure<RS, RS, true>,
    >,
    p: RS::Elem, // An irreducible element of Ring
    k: Natural,
    h: Polynomial<RS::Elem>, // A polynomial over Ring
    h_deg: usize,
    fs_count: usize,
    fs: Vec<Polynomial<RS::Elem>>, // monic Polynomials over Ring mod p^k
    fs_deg: Vec<usize>,
    // A matrix used to compute lifts. We cache it here since it is unchanged during lifts.
    lifting_mat_inv_mod_p: Matrix<RS::Elem>,
    // The inverse of the leading coefficient of h modulo p
    lc_inv_mod_p: RS::Elem,
}

impl<
    RS: EuclideanDomainSignature
        + GreatestCommonDivisorSignature
        + FactoringMonoidSignature<FactoredExponent = NaturalCanonicalStructure>,
> HenselFactorization<RS>
{
    pub fn ring(&self) -> &RS {
        &self.ring
    }

    pub fn modulus(&self) -> RS::Elem {
        self.ring().nat_pow(&self.p, &self.k)
    }

    fn lc(&self) -> &RS::Elem {
        self.ring.polynomials().leading_coeff(&self.h).unwrap()
    }

    fn compute_lifting_mat_inv_mod_p(
        fs_count: usize,
        fs: &[Polynomial<RS::Elem>],
        fs_deg: &[usize],
        h_deg: usize,
        ring: &RS,
        polys_mod_p: &PolynomialStructure<
            EuclideanRemainderQuotientStructure<RS, RS, true>,
            EuclideanRemainderQuotientStructure<RS, RS, true>,
        >,
        mats_mod_p: &MatrixStructure<
            EuclideanRemainderQuotientStructure<RS, RS, true>,
            EuclideanRemainderQuotientStructure<RS, RS, true>,
        >,
    ) -> Matrix<RS::Elem> {
        debug_assert_eq!(fs_count, fs.len());
        debug_assert_eq!(fs_count, fs_deg.len());
        debug_assert_eq!(fs_deg.iter().cloned().sum::<usize>(), h_deg);

        println!("make prods");

        let fs_prod = polys_mod_p.product(fs);
        let fs_prod_excluding_each_mod_p = (0..fs_count)
            .map(|i| polys_mod_p.try_divide(&fs_prod, &fs[i]).unwrap())
            .collect::<Vec<_>>();

        println!("make {}x{} mat", h_deg, h_deg);

        let mat = Matrix::<RS::Elem>::from_cols({
            let mut cols = vec![];
            for idx in 0..fs_count {
                let d = fs_deg[idx];
                for i in 0..d {
                    let mut col = vec![];
                    for _ in 0..i {
                        col.push(ring.zero());
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
                        col.push(ring.zero());
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
                &self.ring,
                &self.polys_mod_p,
                &self.mats_mod_p,
            )
        ));

        if !self.ring.is_irreducible(&self.p) {
            return Err(format!("p={:?} is not irreducible", self.p));
        }
        if self.k == Natural::ZERO {
            return Err("k=0 is not valid".to_string());
        }

        let ring_mod_p = self.ring.quotient_field(self.p.clone()).unwrap();
        if !ring_mod_p.is_unit(self.lc()) {
            return Err("leading coeff is not invertible modulo p".to_string());
        }

        let ring_poly = self.ring.polynomials();
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
            .ring
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

        Ok(())
    }

    pub fn new_unchecked(
        ring: RS,
        p: RS::Elem,
        k: Natural,
        h: Polynomial<RS::Elem>,
        fs: Vec<Polynomial<RS::Elem>>,
    ) -> Self {
        let ring_mod_p = ring.clone().into_quotient_field(p.clone()).unwrap();
        let polys = ring.clone().into_polynomials();
        let polys_mod_p = ring_mod_p.clone().into_polynomials();
        let mats_mod_p = MatrixStructure::new(ring_mod_p.clone());

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
            &ring,
            &polys_mod_p,
            &mats_mod_p,
        );

        let lc_inv_mod_p = ring_mod_p
            .try_reciprocal(polys.leading_coeff(&h).unwrap())
            .unwrap();

        let s = Self {
            ring,
            polys,
            polys_mod_p,
            mats_mod_p,
            p,
            k,
            h,
            h_deg,
            fs_count,
            fs,
            fs_deg,
            lifting_mat_inv_mod_p,
            lc_inv_mod_p,
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

        // mod p^k
        let pk0 = self.modulus();
        // mod p^(k+1)
        let pk1 = self.ring().mul(&pk0, &self.p);
        // polynomials mod p^(k+1)
        let polys_mod_pk1 = self
            .ring()
            .quotient_ring(pk1.clone())
            .unwrap()
            .into_polynomials();

        // compute t = (h - lc * (f1 * f2 * f3)) / (lc * p^k)
        let t_poly = self.polys_mod_p.scalar_mul(
            &polys_mod_pk1
                .sub(
                    &self.h,
                    &polys_mod_pk1.scalar_mul(
                        &polys_mod_pk1.product(&self.fs.iter().collect::<Vec<_>>()),
                        self.lc(),
                    ),
                )
                .apply_map_into(|c| self.ring().try_divide(&c, &pk0).unwrap()),
            &self.lc_inv_mod_p,
        );
        let t_coeffs = (0..self.h_deg)
            .map(|i| self.polys.coeff(&t_poly, i).as_ref().clone())
            .collect::<Vec<_>>();

        // solve the linear system for g1, g2, g3
        let delta_sol = self
            .mats_mod_p
            .apply_col(&self.lifting_mat_inv_mod_p, &t_coeffs);

        // extract the coefficients of g1, g2, g3 from the single vector solution to the linear system and add them to self.fs to lift the factorization from mod p^k to mod p^(k+1)
        let mut offset = 0;
        for i in 0..self.fs.len() {
            let f_deg = self.fs_deg[i];
            let mut delta = vec![];
            for j in 0..f_deg {
                delta.push(delta_sol[offset + j].clone());
            }
            self.polys.add_mut(
                self.fs.get_mut(i).unwrap(),
                &self.polys.scalar_mul(&Polynomial::from_coeffs(delta), &pk0),
            );
            offset += f_deg;
        }

        self.k += Natural::ONE;

        #[cfg(debug_assertions)]
        self.check().unwrap();
    }

    pub fn factors(&self) -> &Vec<Polynomial<RS::Elem>> {
        &self.fs
    }

    /// If the polynomial is squarefree return a hensel factorization, otherwise return None
    pub fn from_mod_p_factorization<
        RSB: BorrowedStructure<RS>,
        RSQB: BorrowedStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>>,
        RSQPB: BorrowedStructure<
            PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
        >,
        NB: BorrowedStructure<NaturalCanonicalStructure>,
    >(
        fs_structure: &FactoringStructure<
            PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
            RSQPB,
            NaturalCanonicalStructure,
            NB,
        >,
        fs: NonZeroFactored<Polynomial<RS::Elem>, Natural>,
        h: Polynomial<RS::Elem>,
    ) -> Option<Self> {
        let poly_ring_mod = fs_structure.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring();
        let ring = ring_mod.ring().clone();
        let mut fs_vec = vec![];
        let (_unit, factors) = fs.into_unit_and_powers();
        for (factor, power) in factors {
            if power == Natural::ONE {
                fs_vec.push(factor);
            } else {
                return None;
            }
        }
        let hensel_factorization = Self::new_unchecked(
            ring,
            ring_mod.modulus().into_owned(),
            Natural::ONE,
            h,
            fs_vec,
        );
        #[cfg(debug_assertions)]
        hensel_factorization.check().unwrap();
        Some(hensel_factorization)
    }

    pub fn from_mod_p_factorization_unchecked<
        RSB: BorrowedStructure<RS>,
        RSQB: BorrowedStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>>,
        RSQPB: BorrowedStructure<
            PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
        >,
        NB: BorrowedStructure<NaturalCanonicalStructure>,
    >(
        fs_structure: &FactoringStructure<
            PolynomialStructure<EuclideanRemainderQuotientStructure<RS, RSB, true>, RSQB>,
            RSQPB,
            NaturalCanonicalStructure,
            NB,
        >,
        fs: NonZeroFactored<Polynomial<RS::Elem>, Natural>,
        h: Polynomial<RS::Elem>,
    ) -> Self {
        let poly_ring_mod = fs_structure.objects().clone();
        let ring_mod = poly_ring_mod.coeff_ring();
        let ring = ring_mod.ring().clone();
        let mut fs_vec = vec![];
        let (_unit, factors) = fs.into_unit_and_powers();
        for (factor, power) in factors {
            debug_assert_eq!(power, Natural::ONE);
            fs_vec.push(factor);
        }
        let hensel_factorization = Self::new_unchecked(
            ring,
            ring_mod.modulus().into_owned(),
            Natural::ONE,
            h,
            fs_vec,
        );
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
            Integer::structure(),
            Integer::from(5),
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
