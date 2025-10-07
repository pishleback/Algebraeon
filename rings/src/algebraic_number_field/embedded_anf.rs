use super::polynomial_quotient_number_field::AlgebraicNumberFieldPolynomialQuotientStructure;
use crate::isolated_algebraic::ComplexAlgebraicCanonicalStructure;
use crate::structure::*;
use crate::{
    isolated_algebraic::{ComplexAlgebraic, RealAlgebraic},
    polynomial::*,
};
use algebraeon_nzq::traits::Fraction;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::rc::Rc;

#[derive(Debug, Clone)]
pub struct EmbeddedAnf {
    //anf.modulus() == gen.min_poly()
    anf: Rc<AlgebraicNumberFieldPolynomialQuotientStructure>,
    generator: ComplexAlgebraic,
}

impl AlgebraicNumberFieldPolynomialQuotientStructure {
    pub fn all_complex_embeddings(&self) -> Vec<EmbeddedAnf> {
        self.modulus()
            .primitive_part_fof()
            .all_complex_roots()
            .into_iter()
            .map(|generator| EmbeddedAnf {
                anf: self.clone().into(),
                generator,
            })
            .collect()
    }

    pub fn all_real_embeddings(&self) -> Vec<EmbeddedAnf> {
        self.modulus()
            .primitive_part_fof()
            .all_real_roots()
            .into_iter()
            .map(|generator| EmbeddedAnf {
                anf: self.clone().into(),
                generator: ComplexAlgebraic::from(generator),
            })
            .collect()
    }

    /// Return (r, s) where r is the number of real embeddings and s is the number of pairs of complex embeddings
    pub fn signature(&self) -> (usize, usize) {
        let poly = self.modulus();
        let d = poly.degree().unwrap();
        let r = poly.count_real_roots();
        let two_s = d - r;
        debug_assert_eq!(two_s % 2, 0);
        (r, two_s / 2)
    }

    pub fn is_totally_real(&self) -> bool {
        self.signature().1 == 0
    }

    pub fn is_totally_positive(&self, a: Polynomial<Rational>) -> bool {
        self.all_real_embeddings()
            .into_iter()
            .all(|mut v| match v.generator.apply_poly(&a) {
                ComplexAlgebraic::Real(v_image_real) => v_image_real > RealAlgebraic::zero(),
                _ => false,
            })
    }
}

impl ComplexAlgebraic {
    pub fn generated_algebraic_number_field(
        &self,
    ) -> AlgebraicNumberFieldPolynomialQuotientStructure {
        self.min_poly().algebraic_number_field_unchecked()
    }

    pub fn embedded_generated_algebraic_number_field(self) -> EmbeddedAnf {
        EmbeddedAnf {
            anf: self.generated_algebraic_number_field().into(),
            generator: self,
        }
    }
}

//TODO
// use complex::*;
// use  real::*;

impl AlgebraicClosureSignature for ComplexAlgebraicCanonicalStructure {
    type BFS = <Rational as algebraeon_sets::structure::MetaType>::Signature;

    fn base_field(&self) -> Self::BFS {
        Rational::structure()
    }

    fn base_field_inclusion(&self, x: &Rational) -> Self::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(x.clone()))
    }

    fn all_roots_list(&self, poly: &Polynomial<Rational>) -> Option<Vec<Self::Set>> {
        if poly.is_zero() {
            None
        } else {
            Some(poly.primitive_part_fof().all_complex_roots())
        }
    }
}

impl EmbeddedAnf {
    pub fn is_real(&self) -> bool {
        matches!(self.generator, ComplexAlgebraic::Real(_))
    }
    pub fn is_rational(&self) -> bool {
        self.is_real()
            && match &self.generator {
                ComplexAlgebraic::Real(real_algebraic) => match real_algebraic {
                    RealAlgebraic::Rational(_) => true,
                    RealAlgebraic::Real(_) => false,
                },
                ComplexAlgebraic::Complex(_) => false,
            }
    }

    pub fn embed(
        &self,
        element: &<AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set,
    ) -> ComplexAlgebraic {
        let complex_alg_canonical = ComplexAlgebraicCanonicalStructure {};
        let mut answer: ComplexAlgebraic =
            ComplexAlgebraic::Real(RealAlgebraic::Rational(element.coeff(0).into_owned()));
        let mut x_to_idx = self.generator.clone();
        for coeff in element.coeffs().into_iter().skip(1) {
            let cur_contribution: ComplexAlgebraic = complex_alg_canonical.mul(
                &x_to_idx,
                &ComplexAlgebraic::Real(RealAlgebraic::Rational(coeff.clone())),
            );
            complex_alg_canonical.add_mut(&mut answer, &cur_contribution);
            complex_alg_canonical.mul_mut(&mut x_to_idx, &self.generator);
        }
        answer
    }
}

#[cfg(any())]
impl EmbeddedAnf {
    pub fn intersect_pair(field1: &Self, field2: &Self) -> Self {
        todo!()
    }

    pub fn intersect_list(fields: Vec<impl Borrow<Self>>) -> Self {
        todo!()
    }

    pub fn generated_pair(field1: &Self, field2: &Self) -> Self {
        todo!()
    }

    pub fn generated_list(fields: Vec<impl Borrow<Self>>) -> Self {
        todo!()
    }
}

//write target as a polynomial expression of generator if possible, otherwise return None
pub fn as_poly_expr(
    target: &ComplexAlgebraic,
    generator: &ComplexAlgebraic,
) -> Option<Polynomial<Rational>> {
    //idea: factor the minimal polynomial of target in the algebraic number field generated by the generator
    //loop through the linear factors (x - a) and check if a = target

    //let K = Q[generator]
    let gen_anf = generator.generated_algebraic_number_field();
    let gen_anf_poly = gen_anf.polynomial_ring();

    //the minimal polynomial of target in K[x]
    let target_min_poly = target
        .min_poly()
        .apply_map(|c| Polynomial::constant(c.clone()));

    let target_min_poly_factored = gen_anf_poly.factor(&target_min_poly).unwrap();
    let mut generator = generator.clone();
    for (factor, _factor_mult) in gen_anf_poly
        .factorizations()
        .to_powers(&target_min_poly_factored)
    {
        //the factor should be monic
        debug_assert!(gen_anf.equal(gen_anf_poly.leading_coeff(factor).unwrap(), &gen_anf.one()));
        if gen_anf_poly.degree(factor).unwrap() == 1 {
            let possible_embedded_target = gen_anf.neg(gen_anf_poly.coeff(factor, 0).as_ref());
            if generator.apply_poly(&possible_embedded_target) == *target {
                return Some(possible_embedded_target);
            }
        }
    }
    None
}

pub fn anf_pair_primitive_element_theorem(
    a: &ComplexAlgebraic,
    b: &ComplexAlgebraic,
) -> (
    ComplexAlgebraic,
    Integer,
    Integer,
    Polynomial<Rational>,
    Polynomial<Rational>,
) {
    //try g = a
    if let Some(q) = as_poly_expr(b, a) {
        return (a.clone(), Integer::ONE, Integer::ZERO, Polynomial::var(), q);
    }

    //try g = b
    if let Some(p) = as_poly_expr(a, b) {
        return (b.clone(), Integer::ZERO, Integer::ONE, p, Polynomial::var());
    }

    let mut nontrivial_linear_combinations = Rational::exhaustive_rationals().map(|r| {
        let (n, d) = r.numerator_and_denominator();
        (n, Integer::from(d))
    });
    nontrivial_linear_combinations.next().unwrap();
    for (x, y) in nontrivial_linear_combinations {
        let generator = ComplexAlgebraic::add(
            &ComplexAlgebraic::mul(
                &ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(x.clone()))),
                a,
            ),
            &ComplexAlgebraic::mul(
                &ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(y.clone()))),
                b,
            ),
        );

        if let Some(a_rel_gen) = as_poly_expr(a, &generator) {
            let anf = generator.generated_algebraic_number_field();
            //gen = xa + yb
            //so b = (gen - xa) / y
            let b_rel_gen = anf.mul(
                &anf.add(
                    &Polynomial::var(),
                    &anf.mul(&a_rel_gen, &Polynomial::constant(Rational::from(-&x))),
                ),
                &Polynomial::constant(Rational::from_integers(Integer::from(1), y.clone())),
            );
            #[cfg(debug_assertions)]
            {
                let mut gen_mut = generator.clone();
                assert_eq!(a, &gen_mut.apply_poly(&a_rel_gen));
                assert_eq!(b, &gen_mut.apply_poly(&b_rel_gen));
            }
            return (generator, x, y, a_rel_gen, b_rel_gen);
        }
    }
    unreachable!()
}

/*
input: non-empty list of complex algebraic numbers (a_1, a_2, ..., a_n)
output: (g, p_1, p_2, ..., p_n) such that Q[a_1, a_2, ..., a_n] = Q[g]
        moreover a_i=p_i(g)
*/
pub fn anf_multi_primitive_element_theorem(
    nums: Vec<&ComplexAlgebraic>,
) -> (ComplexAlgebraic, Vec<Polynomial<Rational>>) {
    #[cfg(debug_assertions)]
    let orig_nums = nums.clone();

    assert!(!nums.is_empty());
    let mut nums = nums.into_iter();
    let mut g = nums.next().unwrap().clone();
    let mut p = vec![Polynomial::var()];
    for num in nums {
        let (new_g, _x, _y, old_g_poly, num_poly) = anf_pair_primitive_element_theorem(&g, num);
        let new_g_anf = new_g.generated_algebraic_number_field();
        p = p
            .into_iter()
            .map(|old_p| new_g_anf.reduce(Polynomial::compose(&old_p, &old_g_poly)))
            .collect();
        p.push(num_poly);
        g = new_g;
    }
    #[cfg(debug_assertions)]
    {
        let n = orig_nums.len();
        assert_eq!(n, p.len());
        for i in 0..n {
            assert_eq!(orig_nums[i], &g.apply_poly(&p[i]));
        }
    }
    (g, p)
}

#[cfg(test)]
mod tests {
    use crate::{
        algebraic_number_field::structure::AlgebraicNumberFieldSignature, structure::IntoErgonomic,
    };

    use super::*;

    #[test]
    fn test_embedded_anf_integral_basis() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (x.pow(2) + 7).into_verbose();
        for root in f.all_complex_roots() {
            println!(
                "{:?}",
                root.generated_algebraic_number_field()
                    .into_ring_of_integers_extension()
            );
        }
    }

    #[test]
    fn test_as_poly_expr() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();

        let two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)));
        let sqrt_two = two.nth_root(2).unwrap();
        let fourrt_two = sqrt_two.nth_root(2).unwrap();
        assert_eq!(
            as_poly_expr(&sqrt_two, &fourrt_two),
            Some(x.pow(2).into_verbose())
        );
        assert_eq!(as_poly_expr(&fourrt_two, &sqrt_two), None);

        let f = (x.pow(3) - 2 * x.pow(2) - 2 * x - 1).into_verbose();
        for root in f.primitive_part_fof().all_complex_roots() {
            assert_eq!(as_poly_expr(&root, &root), Some(x.pow(1).into_verbose()));
        }
    }

    #[test]
    fn test_pair_generated_anf() {
        // let x = &Polynomial::<Rational>::var().into_ergonomic();

        let sqrt_two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)))
            .nth_root(2)
            .unwrap();
        let sqrt_three = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(3)))
            .nth_root(2)
            .unwrap();
        let sqrt_six = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(6)))
            .nth_root(2)
            .unwrap();

        println!("{}", sqrt_two);
        println!("{}", sqrt_three);
        println!("{}", sqrt_six);

        let (generator, _, _, x, y) = anf_pair_primitive_element_theorem(&sqrt_two, &sqrt_three);
        println!(
            "gen = {} min_poly = {} deg = {}",
            generator,
            generator.min_poly(),
            generator.min_poly().degree().unwrap()
        );
        assert_eq!(
            generator,
            (sqrt_two.clone().into_ergonomic() + sqrt_three.clone().into_ergonomic())
                .into_verbose()
        );
        let oof = (sqrt_two.clone().into_ergonomic() + sqrt_three.clone().into_ergonomic())
            .into_verbose();
        println!("{} {}", oof, oof.min_poly());
        println!("x = {}", x);
        println!("y = {}", y);
    }

    #[test]
    fn test_multi_generated_anf() {
        // let x = &Polynomial::<Rational>::var().into_ergonomic();

        let sqrt_two = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(2)))
            .nth_root(2)
            .unwrap();
        let sqrt_three = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(3)))
            .nth_root(2)
            .unwrap();
        let sqrt_six = ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from(6)))
            .nth_root(2)
            .unwrap();

        println!("{}", sqrt_two);
        println!("{}", sqrt_three);
        println!("{}", sqrt_six);

        let (mut g, p) =
            anf_multi_primitive_element_theorem(vec![&sqrt_two, &sqrt_three, &sqrt_six]);

        println!("g = {} of degree {}", g, g.degree());
        println!("sqrt(2) = {} evaluated at {}", p[0], g);
        println!("sqrt(3) = {} evaluated at {}", p[1], g);
        println!("sqrt(6) = {} evaluated at {}", p[2], g);

        assert_eq!(sqrt_two, g.apply_poly(&p[0]));
        assert_eq!(sqrt_three, g.apply_poly(&p[1]));
        assert_eq!(sqrt_six, g.apply_poly(&p[2]));
    }
}
