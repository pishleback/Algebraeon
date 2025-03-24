use std::{borrow::Borrow, collections::HashSet};

use itertools::Itertools;

use crate::structure::structure::*;
use algebraeon_sets::structure::*;

use super::multipoly::*;

pub fn ss_num(n: usize) -> String {
    let mut ss = String::new();
    for x in n.to_string().chars() {
        match x {
            '0' => ss += "₀",
            '1' => ss += "₁",
            '2' => ss += "₂",
            '3' => ss += "₃",
            '4' => ss += "₄",
            '5' => ss += "₅",
            '6' => ss += "₆",
            '7' => ss += "₇",
            '8' => ss += "₈",
            '9' => ss += "₉",
            c => ss += &mut c.to_string(),
        }
    }
    ss
}

//express poly as a polynomial in the elementary symmetric polynomials in the given variables
//return none if poly is not symmetric in the given variables
impl<RS: IntegralDomainStructure> MultiPolynomialStructure<RS>
//TODO: replace integral domain with division ring structure
where
    MultiPolynomialStructure<RS>: Structure<Set = MultiPolynomial<RS::Set>> + ToStringStructure,
{
    pub fn is_symmetric(
        &self,
        vars: Vec<impl Borrow<Variable>>,
        poly: &MultiPolynomial<RS::Set>,
    ) -> bool {
        let vars_hash: HashSet<_> = vars.iter().map(|v| v.borrow().clone()).collect();
        for v in poly.free_vars() {
            if !vars_hash.contains(&v) {
                return false;
            }
        }

        let n = vars.len();
        //find a generating set of Sn
        let perms = match n {
            0 => {
                //trivial
                vec![]
            }
            1 => {
                //trivial
                vec![vec![0]]
            }
            2 => {
                // (1 2)
                vec![
                    (0..n)
                        .map(|i| match i {
                            0 => 1,
                            1 => 0,
                            x => x,
                        })
                        .collect(),
                ]
            }
            n => {
                // (1 2) and (1 2 ... n)
                vec![
                    (0..n)
                        .map(|i| match i {
                            0 => 1,
                            1 => 0,
                            x => x,
                        })
                        .collect(),
                    (0..n).map(|i| (i + 1) % n).collect(),
                ]
            }
        };

        perms.into_iter().all(|perm| {
            let perm_poly = poly.clone().apply_map_vars(
                perm.into_iter()
                    .enumerate()
                    .map(|(i, j)| (vars[i].borrow().clone(), vars[j].borrow().clone()))
                    .collect(),
            );
            self.equal(poly, &perm_poly)
        })
    }

    pub fn elementary_symmetric(
        &self,
        n: usize,
        vars: &Vec<impl Borrow<Variable>>,
    ) -> MultiPolynomial<RS::Set> {
        let mp_ring = MultiPolynomialStructure::new(self.coeff_ring().clone());
        let mut e = mp_ring.zero();
        for vs in algebraeon_sets::subset::subsets_vec(vars.iter().collect(), n) {
            let mut p = mp_ring.one();
            for v in vs {
                mp_ring.mul_mut(&mut p, &mp_ring.var(v.borrow().clone()));
            }
            mp_ring.add_mut(&mut e, &p);
        }
        e
    }

    fn as_elementary_symmetric_polynomials_homogeneous_impl(
        &self,
        vars: &Vec<Variable>,
        p: &MultiPolynomial<RS::Set>,
        e: &Vec<Variable>,
    ) -> MultiPolynomial<RS::Set> {
        // println!("e = {:?}", e);
        // println!("p = {}", self.elem_to_string(p));

        let (last_var, first_vars) = vars.split_last().unwrap();

        let p_tilde = p.clone().evaluate_var_zero(last_var);
        let r_sym = self.as_elementary_symmetric_polynomials_impl(
            &first_vars.iter().map(|v| v.clone()).collect(),
            &p_tilde,
            e,
        );
        // println!("first_vars={:?} last_var={:?}", first_vars, last_var);
        // println!("r_sym = {}", self.elem_to_string(&r_sym));
        let r = MultiPolynomialStructure::new(self.clone().into()).evaluate(
            &r_sym.apply_map(|x| MultiPolynomial::constant(x.clone())),
            (0..first_vars.len())
                .map(|i| (e[i].clone(), self.elementary_symmetric(i + 1, vars)))
                .collect(),
        );
        // println!("p_tilde = {}", self.elem_to_string(&p_tilde));
        // println!("p_tilde_sym = {}", self.elem_to_string(&r_sym));
        // println!("r = {}", self.elem_to_string(&r));
        //q = p - r  is divisible by x1*x2*...*xn = en
        let q = self.add(p, &self.neg(&r));
        // println!("q = {}", self.elem_to_string(&q));
        let en = self.product(vars.iter().map(|v| self.var(v.clone())).collect());
        let en_sym = self.var(e.get(vars.len() - 1).unwrap().clone());
        // println!("en = {}", self.elem_to_string(&en));
        // println!("en_sym = {}", self.elem_to_string(&en_sym));
        //p = r + en * q

        let p_sym = self.add(
            &r_sym,
            &self.mul(
                &en_sym,
                &self.as_elementary_symmetric_polynomials_impl(
                    vars,
                    &self.div(&q, &en).unwrap(),
                    e,
                ),
            ),
        );

        // println!("p_sym = {}", self.elem_to_string(&p_sym));

        p_sym
    }

    fn as_elementary_symmetric_polynomials_impl(
        &self,
        vars: &Vec<Variable>,
        poly: &MultiPolynomial<RS::Set>,
        e: &Vec<Variable>,
    ) -> MultiPolynomial<RS::Set> {
        let mut total = self.zero();
        for (d, hom_poly) in self.split_by_degree(poly.clone()) {
            if d == 0 {
                self.add_mut(&mut total, &hom_poly);
            } else {
                self.add_mut(
                    &mut total,
                    &self.as_elementary_symmetric_polynomials_homogeneous_impl(vars, &hom_poly, &e),
                );
            }
        }
        total
    }

    //assume input is symmetrical
    pub fn as_elementary_symmetric_polynomials_unchecked(
        &self,
        vars: Vec<impl Borrow<Variable>>,
        poly: &MultiPolynomial<RS::Set>,
    ) -> (Vec<Variable>, MultiPolynomial<RS::Set>) {
        let e = (0..vars.len())
            .map(|i| Variable::new(format!("e{}", ss_num(i + 1))))
            .collect_vec();

        #[cfg(debug_assertions)]
        let vars_clone = vars.iter().map(|v| v.borrow().clone()).collect();

        let poly_sym = self.as_elementary_symmetric_polynomials_impl(
            &vars.iter().map(|v| v.borrow().clone()).collect(),
            poly,
            &e,
        );

        #[cfg(debug_assertions)]
        {
            let poly_check = MultiPolynomialStructure::new(self.clone().into()).evaluate(
                &poly_sym.apply_map(|x| MultiPolynomial::constant(x.clone())),
                (0..e.len())
                    .map(|i| (e[i].clone(), self.elementary_symmetric(i + 1, &vars_clone)))
                    .collect(),
            );
            assert!(self.equal(poly, &poly_check));
        }

        (e, poly_sym)
    }

    //return None if not symmetrical
    pub fn as_elementary_symmetric_polynomials(
        &self,
        vars: Vec<impl Borrow<Variable>>,
        poly: &MultiPolynomial<RS::Set>,
    ) -> Option<(Vec<Variable>, MultiPolynomial<RS::Set>)> {
        if !self.is_symmetric(vars.iter().map(|v| v.borrow().clone()).collect(), poly) {
            None
        } else {
            Some(self.as_elementary_symmetric_polynomials_unchecked(vars, poly))
        }
    }
}

impl<R: MetaType> MultiPolynomial<R>
where
    R::Structure: IntegralDomainStructure,
    MultiPolynomialStructure<R::Structure>: Structure<Set = Self> + ToStringStructure,
{
    pub fn is_symmetric(&self, vars: Vec<impl Borrow<Variable>>) -> bool {
        Self::structure().is_symmetric(vars, self)
    }

    pub fn elementary_symmetric(n: usize, vars: Vec<impl Borrow<Variable>>) -> MultiPolynomial<R> {
        Self::structure().elementary_symmetric(n, &vars)
    }

    pub fn as_elementary_symmetric_polynomials(
        &self,
        vars: Vec<impl Borrow<Variable>>,
    ) -> Option<(Vec<Variable>, MultiPolynomial<R>)> {
        Self::structure().as_elementary_symmetric_polynomials(vars, self)
    }

    pub fn as_elementary_symmetric_polynomials_unchecked(
        &self,
        vars: Vec<impl Borrow<Variable>>,
    ) -> (Vec<Variable>, MultiPolynomial<R>) {
        Self::structure().as_elementary_symmetric_polynomials_unchecked(vars, self)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::structure::elements::*;

    use algebraeon_nzq::integer::*;

    #[test]
    fn test_ffosp() {
        let x_var = Variable::new("x");
        let y_var = Variable::new("y");
        let z_var = Variable::new("z");

        // let mut espm =
        //     MultiPolynomial::<Integer>::symmetric_polynomial_manager(vec![&x_var, &y_var, &z_var]);

        // println!(
        //     "e0 = {}",
        //     MultiPolynomial::<Integer>::elementary_symmetric(0, vec![&x_var, &y_var, &z_var])
        // );
        // println!(
        //     "e1 = {}",
        //     MultiPolynomial::<Integer>::elementary_symmetric(1, vec![&x_var, &y_var, &z_var])
        // );
        // println!(
        //     "e2 = {}",
        //     MultiPolynomial::<Integer>::elementary_symmetric(2, vec![&x_var, &y_var, &z_var])
        // );
        // println!(
        //     "e3 = {}",
        //     MultiPolynomial::<Integer>::elementary_symmetric(3, vec![&x_var, &y_var, &z_var])
        // );

        let x = &MultiPolynomial::<Integer>::var(x_var.clone()).into_ergonomic();
        let y = &MultiPolynomial::<Integer>::var(y_var.clone()).into_ergonomic();
        let z = &MultiPolynomial::<Integer>::var(z_var.clone()).into_ergonomic();

        let f =
            (x * x * y + x * x * z + y * y * x + y * y * z + z * z * x + z * z * y).into_verbose();

        println!(
            "f={} {}",
            f,
            f.is_symmetric(
                vec![x_var.clone(), y_var.clone(), z_var.clone()]
                    .into_iter()
                    .collect()
            )
        );

        let (e, g) = f
            .as_elementary_symmetric_polynomials(vec![&x_var, &y_var, &z_var])
            .unwrap();

        println!("{:?}", e);
        println!("{}", g);
    }

    #[test]
    fn run() {
        let x_var = Variable::new("x");
        let y_var = Variable::new("y");
        let z_var = Variable::new("z");

        let x = &MultiPolynomial::<Integer>::var(x_var.clone()).into_ergonomic();
        let y = &MultiPolynomial::<Integer>::var(y_var.clone()).into_ergonomic();
        let z = &MultiPolynomial::<Integer>::var(z_var.clone()).into_ergonomic();

        let f = (x.pow(12) + y.pow(12) + z.pow(12)).into_verbose();

        let (e, g) = f
            .as_elementary_symmetric_polynomials(vec![&x_var, &y_var, &z_var])
            .unwrap();

        println!("{:?}", e);
        println!("{}", g);
    }
}
