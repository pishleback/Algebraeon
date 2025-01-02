use algebraeon_structure::*;
use itertools::Itertools;
use malachite_nz::natural::Natural;

use crate::{
    linear::{matrix::*, subspace::*},
    polynomial::polynomial::*,
    ring_structure::{factorization::*, structure::*},
};

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// Reduce a factorization problem for polynomials over a finite field to factorizations of squarefree polynomials over the finite field
    // https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
    pub fn factorize_by_sqfree_factorize_over_finite_field(
        &self,
        f: Polynomial<FS::Set>,
        sqfree_factorize: &impl Fn(Polynomial<FS::Set>) -> Factored<PolynomialStructure<FS>>,
    ) -> Factored<PolynomialStructure<FS>>
    where
        PolynomialStructure<FS>: EuclideanDivisionStructure
            + GreatestCommonDivisorStructure
            + UniqueFactorizationStructure,
    {
        let (unit, f) = self.factor_fav_assoc(&f);
        // println!("f = {:?}", f);
        let mut factors = Factored::new_unchecked(self.clone().into(), unit, vec![]);
        //let w be the squarefree product of all factors of f of multiplicity not divisible by p
        let mut c = self.gcd(&f, &self.derivative(f.clone()));
        let mut w = self.div(&f, &c).unwrap();
        // println!("c = {:?}", c);
        // println!("w = {:?}", w);
        // println!("{:?}", self.divisible(&f, &c));
        // println!("{:?}", self.divisible(&f, &w));
        // println!("f' = {:?}", self.derivative(f.clone()));

        //find all the factors in w
        let mut i = Natural::from(1u8);
        while !self.equal(&w, &self.one()) {
            let y = self.gcd(&w, &c);
            factors.mul_mut(sqfree_factorize(self.div(&w, &y).unwrap()).pow(&i));
            let c_over_y = self.div(&c, &y).unwrap();
            (w, c) = (y, c_over_y);
            i += Natural::from(1u8);
        }

        // c is now the product, with multiplicity, of the remaining factors of f
        if !self.equal(&c, &self.one()) {
            //c = c^{1/p}
            let p = self.coeff_ring().characteristic_and_power().0;
            // println!("p = {:?} c = {:?}", p, c);
            let mut reduced_c_coeffs = vec![];
            for (k, coeff) in c.coeffs().into_iter().enumerate() {
                if Natural::from(k) % &p == 0 {
                    reduced_c_coeffs.push(coeff.clone());
                } else {
                    debug_assert!(self.coeff_ring().is_zero(coeff));
                }
            }
            let reduced_c = Polynomial::from_coeffs(reduced_c_coeffs);
            // println!("reduced_c = {}", reduced_c);
            factors.mul_mut(
                self.factorize_by_sqfree_factorize_over_finite_field(reduced_c, sqfree_factorize)
                    .pow(&p),
            );
        }
        factors
    }
}

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    fn find_factor_by_berlekamps_algorithm(
        &self,
        f: Polynomial<FS::Set>,
    ) -> FindFactorResult<PolynomialStructure<FS>> {
        //f is squarefree
        let f_deg = self.degree(&f).unwrap();
        let all_elems = self.coeff_ring().all_elements();
        let q = all_elems.len();

        let row_polys = (0..f_deg)
            .map(|i| {
                self.rem(
                    &self.add(&self.var_pow(i * q), &self.neg(&self.var_pow(i))),
                    &f,
                )
            })
            .collect::<Vec<_>>();
        let mat = Matrix::construct(f_deg, f_deg, |i, j| self.coeff(&row_polys[i], j).clone());
        // mat.pprint();
        //the column kernel gives a basis of berlekamp subalgebra - all polynomials g such that g^q=g
        let mat_struct = MatrixStructure::<FS>::new(self.coeff_ring().into());
        let linlat_struct = LinearLatticeStructure::<FS>::new(self.coeff_ring().into());
        let ker = mat_struct.row_kernel(mat);
        // ker.pprint();
        let ker_rank = linlat_struct.rank(&ker);
        let ker_basis = linlat_struct
            .basis_matrices(&ker)
            .into_iter()
            .map(|col| {
                Polynomial::from_coeffs((0..f_deg).map(|c| col.at(0, c).unwrap().clone()).collect())
            })
            .collect::<Vec<_>>();

        for berlekamp_subspace_coeffs in (0..ker_rank)
            .map(|_i| all_elems.clone().into_iter())
            .multi_cartesian_product()
        {
            let h = self.sum(
                (0..ker_rank)
                    .map(|i| {
                        self.mul(
                            &Polynomial::constant(berlekamp_subspace_coeffs[i].clone()),
                            &ker_basis[i],
                        )
                    })
                    .collect(),
            );
            //g is a possible non-trivial factor
            let g = self.gcd(&h, &f);
            // println!("g = {}", g);
            let g_deg = self.degree(&g).unwrap();
            if g_deg != 0 && g_deg != f_deg {
                match self.div(&f, &g) {
                    Ok(g_prime) => {
                        return FindFactorResult::Composite(g, g_prime);
                    }
                    Err(_) => panic!(),
                }
            }
        }
        FindFactorResult::Irreducible
    }

    pub fn factorize_by_berlekamps_algorithm(
        &self,
        f: Polynomial<FS::Set>,
    ) -> Option<Factored<PolynomialStructure<FS>>>
    where
        PolynomialStructure<FS>: UniqueFactorizationStructure,
    {
        if self.is_zero(&f) {
            None
        } else {
            Some(
                self.factorize_by_sqfree_factorize_over_finite_field(f, &|f| {
                    factorize_by_find_factor(self, f, &|f| {
                        self.find_factor_by_berlekamps_algorithm(f)
                    })
                }),
            )
        }
    }
}

impl<F: MetaType> Polynomial<F>
where
    F::Structure: FiniteFieldStructure,
    PolynomialStructure<F::Structure>:
        Structure<Set = Polynomial<F>> + UniqueFactorizationStructure,
{
    pub fn factorize_by_berlekamps_algorithm(
        &self,
    ) -> Option<Factored<PolynomialStructure<F::Structure>>> {
        Self::structure().factorize_by_berlekamps_algorithm(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use crate::{
        number::finite_fields::{modulo::*, quaternary_field::QuaternaryField}, ring_structure::elements::IntoRingElem,
    };

    use super::*;

    #[test]
    fn test_factorize_over_f2_example1() {
        let x = &Polynomial::<Modulo<2>>::var().into_ring();
        let p = (1 + x.pow(4) + x.pow(5)).pow(12).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<2>>::structure().into(),
            Polynomial::one(),
            vec![
                ((1 + x + x.pow(2)).into_set(), Natural::from(12u32)),
                ((1 + x + x.pow(3)).into_set(), Natural::from(12u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f2_example2() {
        let x = &Polynomial::<Modulo<2>>::var().into_ring();
        let p = ((1 + x.pow(4) + x.pow(7)).pow(6) * (1 + x.pow(6) + x.pow(7)).pow(4)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<2>>::structure().into(),
            Polynomial::one(),
            vec![
                ((1 + x.pow(4) + x.pow(7)).into_set(), Natural::from(6u32)),
                ((1 + x.pow(6) + x.pow(7)).into_set(), Natural::from(4u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example1() {
        let x = &Polynomial::<Modulo<5>>::var().into_ring();
        let p = (1 + x.pow(4)).pow(5).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<5>>::structure().into(),
            Polynomial::one(),
            vec![
                ((2 + x.pow(2)).into_set(), Natural::from(5u8)),
                ((3 + x.pow(2)).into_set(), Natural::from(5u8)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f5_example2() {
        let x = &Polynomial::<Modulo<5>>::var().into_ring();
        let p = (3 + 2 * x.pow(2) + x.pow(4) + x.pow(6)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<5>>::structure().into(),
            Polynomial::one(),
            vec![((2 + x.pow(2)).into_set(), Natural::from(3u8))],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f31_example1() {
        let x = &Polynomial::<Modulo<31>>::var().into_ring();
        let p = (1 + x.pow(27) + 8 * x.pow(30)).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<Modulo<31>>::structure().into(),
            Polynomial::constant(Modulo::from_int(&Integer::from(8))),
            vec![
                ((12 + x.pow(3)).into_set(), Natural::from(1u32)),
                (
                    (25 + 27 * x + 3 * x.pow(2) + 3 * x.pow(3) + 29 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (1 + 24 * x + 3 * x.pow(2) + 15 * x.pow(3) + 12 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (21 + 12 * x.pow(3) + 22 * x.pow(6) + 4 * x.pow(9) + x.pow(12))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
                (
                    (5 + 11 * x + 3 * x.pow(2) + 13 * x.pow(3) + 21 * x.pow(4) + x.pow(5))
                        .into_set()
                        .clone(),
                    Natural::from(1u32),
                ),
            ],
        );

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }

    #[test]
    fn test_factorize_over_f4_example1() {
        let x = &Polynomial::<QuaternaryField>::var().into_ring();

        let a = x - Polynomial::constant(QuaternaryField::One).into_ring();
        let b = x - Polynomial::constant(QuaternaryField::Alpha).into_ring();
        let c = x - Polynomial::constant(QuaternaryField::Beta).into_ring();
        let p = (1 - x.pow(3)).pow(48).into_set();
        let ans = Factored::new_unchecked(
            Polynomial::<QuaternaryField>::structure().into(),
            Polynomial::one(),
            vec![
                (a.into_set(), Natural::from(48u32)),
                (b.into_set(), Natural::from(48u32)),
                (c.into_set(), Natural::from(48u32)),
            ],
        );

        let f = p.factorize_by_trying_all_factors().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));

        let f = p.factorize_by_berlekamps_algorithm().unwrap();
        println!("{} = {}", p, f);
        assert!(Factored::equal(&f, &ans));
    }
}
