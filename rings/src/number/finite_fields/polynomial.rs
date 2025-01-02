use std::rc::Rc;

use algebraeon_structure::*;
use itertools::Itertools;
use malachite_nz::natural::Natural;

use crate::{
    linear::{matrix::*, subspace::*},
    polynomial::polynomial::*,
    ring_structure::{factorization::*, structure::*},
};

// Useful: https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
/*
Factorization over finite fields typically happens in the following steps:
1. Monic factorization: Remove a unit so that the polynomial to be factored is monic.
2. Squarefree factorization: Factor into squarefree factors.
3. Distinct degree factorization: Factor a squarefree polynomial into polynomials such that all irreducible factors of each polynomial have equal degree.
4. Equal degree factorization: Factor a squarefree polynomial each of whose irreducible factors all have the same known degree d.

1 is trivial.
There is a standard algorithm for 2.
Berlekamps algorithm does 3 and 4 at the same time.
There is a standard algorithm for 3.
Cantor–Zassenhaus algorithm does 4.
*/

/// Store a monic factorization
#[derive(Debug, Clone)]
pub struct MonicFactored<FS: FiniteFieldStructure> {
    poly_ring: PolynomialStructure<FS>,
    unit: FS::Set,              // a unit
    monic: Polynomial<FS::Set>, // a monic polynomial
}

/// Store a squarefree factorization
#[derive(Debug, Clone)]
pub struct SquarefreeFactored<FS: FiniteFieldStructure> {
    poly_ring: PolynomialStructure<FS>,
    unit: FS::Set,                                // a unit
    squarefree_factors: Vec<Polynomial<FS::Set>>, // squarefree monic polynomials
}

#[derive(Debug, Clone)]
struct DistinctDegreeFactor<FS: FiniteFieldStructure> {
    irreducible_factor_degree: usize, // the degree of the irreducible factors
    polynomial: Polynomial<FS::Set>, // a monic squarefree polynomial equal to a product of irreducibles all of the same degree
}
/// Store a distinct degree factorization
#[derive(Debug, Clone)]
pub struct DistinctDegreeFactored<FS: FiniteFieldStructure> {
    poly_ring: PolynomialStructure<FS>,
    unit: FS::Set, // a unit
    distinct_degree_factors: Vec<DistinctDegreeFactor<FS>>,
}

#[derive(Debug, Clone)]
pub struct EqualDegreeFactored<FS: FiniteFieldStructure> {
    poly_ring: PolynomialStructure<FS>,
    unit: FS::Set,                                 // a unit
    irreducible_factors: Vec<Polynomial<FS::Set>>, // a list of irreducible factors
}

impl<FS: FiniteFieldStructure> MonicFactored<FS> {
    pub fn into_polynomial(self) -> Polynomial<FS::Set> {
        self.poly_ring
            .mul(&Polynomial::constant(self.unit), &self.monic)
    }
}
impl<FS: FiniteFieldStructure> SquarefreeFactored<FS> {
    pub fn into_monic_factored(self) -> MonicFactored<FS> {
        let monic = self.poly_ring.product(self.squarefree_factors);
        MonicFactored {
            poly_ring: self.poly_ring,
            unit: self.unit,
            monic,
        }
    }
}
impl<FS: FiniteFieldStructure> DistinctDegreeFactored<FS> {
    pub fn into_squarefree_factored(self) -> SquarefreeFactored<FS> {
        let squarefree_factors = self
            .distinct_degree_factors
            .into_iter()
            .map(|ddf| ddf.polynomial)
            .collect();
        SquarefreeFactored {
            poly_ring: self.poly_ring,
            unit: self.unit,
            squarefree_factors,
        }
    }
}
impl<FS: FiniteFieldStructure> EqualDegreeFactored<FS> {
    pub fn into_distinct_degree_factored(self) -> DistinctDegreeFactored<FS> {
        let distinct_degree_factors = self
            .irreducible_factors
            .into_iter()
            .map(|f| DistinctDegreeFactor {
                irreducible_factor_degree: self.poly_ring.degree(&f).unwrap(),
                polynomial: f,
            })
            .collect();
        DistinctDegreeFactored {
            poly_ring: self.poly_ring,
            unit: self.unit,
            distinct_degree_factors,
        }
    }
    pub fn into_factored(self) -> Factored<PolynomialStructure<FS>>
    where
        PolynomialStructure<FS>:
            Structure<Set = Polynomial<FS::Set>> + UniqueFactorizationStructure,
    {
        let poly_ring: Rc<_> = self.poly_ring.into();
        let mut factored =
            Factored::factored_unit_unchecked(poly_ring.clone(), Polynomial::constant(self.unit));
        for factor in self.irreducible_factors {
            factored.mul_mut(Factored::factored_irreducible_unchecked(
                poly_ring.clone(),
                factor,
            ));
        }
        factored
    }
}
impl<FS: FiniteFieldStructure> Factored<PolynomialStructure<FS>>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>> + UniqueFactorizationStructure,
{
    fn into_equal_degree_factored(self) -> EqualDegreeFactored<FS> {
        let poly_ring = self.ring().as_ref().clone();
        let unit = poly_ring.as_constant(self.unit()).unwrap();
        let irreducible_factors = self.factors_list();
        EqualDegreeFactored {
            poly_ring: poly_ring,
            unit,
            irreducible_factors,
        }
    }
}

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// monic factorization
    pub fn factorize_monic(&self, poly: &Polynomial<FS::Set>) -> MonicFactored<FS> {
        let (unit, monic) = self.factor_fav_assoc(poly);
        let unit = self.as_constant(&unit).unwrap();
        MonicFactored {
            poly_ring: self.clone(),
            unit,
            monic,
        }
    }
}

impl<FS: FiniteFieldStructure> MonicFactored<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// squarefree factorization
    pub fn factorize_squarefree(&self) -> SquarefreeFactored<FS> {
        todo!()
    }
}

impl<FS: FiniteFieldStructure> SquarefreeFactored<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// Berlekamps algorithm for finding a factor of a squarefree polynomial or proving it is irreducible
    pub fn find_factor_berlekamps(&self) -> FindFactorResult<FS> {
        todo!()
    }

    /// use Berlekamps algorithm for a full factorization from a squarefree
    pub fn factorize_berlekamps(&self) -> EqualDegreeFactored<FS> {
        todo!()
    }
}

impl<FS: FiniteFieldStructure> SquarefreeFactored<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// distinct degree factorization
    pub fn factorize_distinct_degree(&self) -> DistinctDegreeFactored<FS> {
        todo!()
    }
}

impl<FS: FiniteFieldStructure> DistinctDegreeFactored<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// Cantor–Zassenhaus algorithm for equal degree factorization
    pub fn factorize_cantor_zassenhaus(&self) -> EqualDegreeFactored<FS> {
        todo!()
    }
}

// old

impl<FS: FiniteFieldStructure> PolynomialStructure<FS>
where
    PolynomialStructure<FS>: Structure<Set = Polynomial<FS::Set>>,
{
    /// Reduce a factorization problem for polynomials over a finite field to factorizations of squarefree polynomials over the finite field
    pub fn factorize_using_sqfree_factorize_over_finite_field(
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
        let mut factors = Factored::new_unchecked(self.clone().into(), unit, vec![]);
        //let w be the squarefree product of all factors of f of multiplicity not divisible by p
        let mut c = self.gcd(&f, &self.derivative(f.clone()));
        let mut w = self.div(&f, &c).unwrap();

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
            let mut reduced_c_coeffs = vec![];
            for (k, coeff) in c.coeffs().into_iter().enumerate() {
                if Natural::from(k) % &p == 0 {
                    reduced_c_coeffs.push(coeff.clone());
                } else {
                    debug_assert!(self.coeff_ring().is_zero(coeff));
                }
            }
            let reduced_c: Polynomial<<FS as Structure>::Set> =
                Polynomial::from_coeffs(reduced_c_coeffs);
            factors.mul_mut(
                self.factorize_using_sqfree_factorize_over_finite_field(
                    reduced_c,
                    sqfree_factorize,
                )
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
        println!("berk ff");

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
                self.factorize_using_sqfree_factorize_over_finite_field(f, &|f| {
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
        number::finite_fields::{modulo::*, quaternary_field::QuaternaryField},
        ring_structure::elements::IntoRingElem,
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
