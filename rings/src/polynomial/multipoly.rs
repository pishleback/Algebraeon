use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use std::rc::Rc;
use std::sync::atomic::AtomicUsize;

use malachite_nz::natural::Natural;

use super::super::ring_structure::cannonical::*;
use super::super::ring_structure::factorization::*;
use super::super::ring_structure::structure::*;
use super::super::structure::*;
use super::polynomial::*;

#[derive(Debug, Hash, Clone)]
pub struct Variable {
    ident: usize,
    name: String,
}

impl PartialEq for Variable {
    fn eq(&self, other: &Self) -> bool {
        self.ident == other.ident
    }
}

impl Eq for Variable {}

impl Variable {
    pub fn new<S: Into<String>>(name: S) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let name = name.into();
        assert!(name.len() >= 1);
        Self {
            ident: COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            name,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct VariablePower {
    var: Variable,
    pow: usize,
}

#[derive(Debug, Clone)]
pub struct Monomial {
    prod: Vec<VariablePower>,            //should be sorted by variable ident
    ident_lookup: HashMap<usize, usize>, //point from variable ident to index of that variables power in self.prod
}

impl PartialEq for Monomial {
    fn eq(&self, other: &Self) -> bool {
        self.prod == other.prod
    }
}

impl Eq for Monomial {}

impl Hash for Monomial {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.prod.hash(state);
    }
}

impl Display for Monomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.prod.len() == 0 {
            write!(f, "1")
        } else {
            for VariablePower { var, pow } in &self.prod {
                write!(f, "{}", &var.name);
                write!(f, "^");
                write!(f, "{}", pow);
            }
            Ok(())
        }
    }
}

impl Monomial {
    fn check_invariants(&self) -> Result<(), &'static str> {
        let mut vars = HashSet::new();
        for VariablePower { var, pow } in &self.prod {
            if pow == &0 {
                return Err("shouldn't have a variable to the power of zero");
            }
            if vars.contains(var) {
                return Err("each var should appear at most once");
            }
            vars.insert(var);
        }
        for (ident, idx) in &self.ident_lookup {
            if &self.prod[*idx].var.ident != ident {
                return Err("bad ident_lookup");
            }
        }
        let mut ordered_prod = self.prod.clone();
        ordered_prod.sort_by_key(|VariablePower { var, pow: _pow }| var.ident);
        if self.prod != ordered_prod {
            return Err("var powers are not sorted");
        }
        Ok(())
    }

    fn new(mut prod: Vec<VariablePower>) -> Self {
        prod = prod
            .into_iter()
            .filter(|VariablePower { var: _var, pow }| pow != &0)
            .collect();
        prod.sort_by_key(|vpow| vpow.var.ident);
        let mut ident_lookup = HashMap::new();
        for (idx, VariablePower { var, pow: _pow }) in prod.iter().enumerate() {
            ident_lookup.insert(var.ident, idx);
        }
        Self { prod, ident_lookup }
    }

    fn one() -> Self {
        Monomial {
            prod: vec![],
            ident_lookup: HashMap::new(),
        }
    }

    fn degree(&self) -> usize {
        let mut d = 0;
        for VariablePower { var: _var, pow } in &self.prod {
            d += pow;
        }
        d
    }

    fn homogenize(&self, target_degree: usize, v: &Variable) -> Self {
        let self_degree = self.degree();
        if target_degree >= self_degree {
            let mut prod = self.prod.clone();
            prod.push(VariablePower {
                var: v.clone(),
                pow: target_degree - self_degree,
            });
            Self::new(prod)
        } else {
            panic!();
        }
    }

    fn get_var_pow(&self, v: &Variable) -> usize {
        if self.ident_lookup.contains_key(&v.ident) {
            self.prod[*self.ident_lookup.get(&v.ident).unwrap()].pow
        } else {
            0
        }
    }

    fn free_vars(&self) -> HashSet<Variable> {
        self.prod
            .iter()
            .map(|VariablePower { var, pow: _pow }| var.clone())
            .collect()
    }

    pub fn evaluate<RS: RingStructure>(
        &self,
        ring: Rc<RS>,
        values: &HashMap<Variable, impl Borrow<RS::Set>>,
    ) -> RS::Set {
        ring.product(
            self.prod
                .iter()
                .map(|VariablePower { var, pow }| {
                    ring.nat_pow(values.get(var).unwrap().borrow(), &Natural::from(*pow))
                })
                .collect(),
        )
    }

    fn mul(a: &Self, b: &Self) -> Self {
        Self::new({
            let mut prod = HashMap::new();
            for VariablePower { var: v, pow: k } in &a.prod {
                *prod.entry(v.clone()).or_insert(0) += k;
            }
            for VariablePower { var: v, pow: k } in &b.prod {
                *prod.entry(v.clone()).or_insert(0) += k;
            }
            let prod: Vec<VariablePower> = prod
                .into_iter()
                .map(|(v, k)| VariablePower { var: v, pow: k })
                .collect();
            prod
        })
    }

    fn lexicographic_order(a: &Self, b: &Self) -> std::cmp::Ordering {
        let mut i = 0;
        while i < std::cmp::min(a.prod.len(), b.prod.len()) {
            if a.prod[i].var.ident < b.prod[i].var.ident {
                return std::cmp::Ordering::Less;
            } else if a.prod[i].var.ident > b.prod[i].var.ident {
                return std::cmp::Ordering::Greater;
            } else {
                if a.prod[i].pow > b.prod[i].pow {
                    return std::cmp::Ordering::Less;
                } else if a.prod[i].pow < b.prod[i].pow {
                    return std::cmp::Ordering::Greater;
                } else {
                    i += 1;
                }
            }
        }
        if a.prod.len() > b.prod.len() {
            return std::cmp::Ordering::Less;
        } else if a.prod.len() < b.prod.len() {
            return std::cmp::Ordering::Greater;
        } else {
            return std::cmp::Ordering::Equal;
        }
    }

    fn graded_lexicographic_order(a: &Self, b: &Self) -> std::cmp::Ordering {
        if a.degree() < b.degree() {
            std::cmp::Ordering::Greater
        } else if a.degree() > b.degree() {
            std::cmp::Ordering::Less
        } else {
            Self::lexicographic_order(a, b)
        }
    }
}

#[derive(Debug, Clone)]
pub struct Term<ElemT: Clone> {
    coeff: ElemT,
    monomial: Monomial,
}

impl<ElemT: Clone> Term<ElemT> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        self.monomial.check_invariants()
    }

    pub fn degree(&self) -> usize {
        self.monomial.degree()
    }
}

#[derive(Debug, Clone)]
pub struct MultiPolynomial<R: Clone> {
    terms: Vec<Term<R>>, //sorted by monomial ordering
}

impl<R: Clone> MultiPolynomial<R> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        for term in &self.terms {
            match term.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            }
            // if self.coeff_ring.is_zero(&term.coeff) {
            //     return Err("coeff should not be zero");
            // }
        }

        if !(0..self.terms.len() - 1).all(|i| {
            Monomial::lexicographic_order(&self.terms[i].monomial, &self.terms[i + 1].monomial)
                .is_le()
        }) {
            return Err("terms are not sorted");
        }

        Ok(())
    }

    fn new(mut terms: Vec<Term<R>>) -> Self {
        terms.sort_by(|t1, t2| Monomial::lexicographic_order(&t1.monomial, &t2.monomial));
        Self { terms }
    }

    pub fn constant(c: R) -> MultiPolynomial<R> {
        MultiPolynomial {
            terms: vec![Term {
                coeff: c,
                monomial: Monomial::one(),
            }],
        }
    }

    pub fn term(t: Term<R>) -> MultiPolynomial<R> {
        MultiPolynomial { terms: vec![t] }
    }

    pub fn free_vars(&self) -> HashSet<Variable> {
        let mut vars = HashSet::new();
        for term in &self.terms {
            vars.extend(term.monomial.free_vars());
        }
        vars
    }

    pub fn apply_map<ImgSet: Clone>(&self, f: impl Fn(&R) -> ImgSet) -> MultiPolynomial<ImgSet> {
        MultiPolynomial {
            terms: self
                .terms
                .iter()
                .map(|Term { coeff, monomial }| Term {
                    coeff: f(coeff),
                    monomial: monomial.clone(),
                })
                .collect(),
        }
    }

    pub fn apply_map_vars(self, f: HashMap<Variable, Variable>) -> MultiPolynomial<R> {
        MultiPolynomial {
            terms: self
                .terms
                .into_iter()
                .map(
                    |Term {
                         coeff,
                         monomial: Monomial { prod, ident_lookup },
                     }| Term {
                        coeff: coeff,
                        monomial: Monomial::new(
                            prod.into_iter()
                                .map(|VariablePower { var, pow }| VariablePower {
                                    var: f.get(&var).unwrap().clone(),
                                    pow,
                                })
                                .collect(),
                        ),
                    },
                )
                .collect(),
        }
    }

    pub fn evaluate_var_zero(self, v: &Variable) -> MultiPolynomial<R> {
        MultiPolynomial {
            terms: self
                .terms
                .into_iter()
                .filter(|Term { coeff, monomial }| monomial.get_var_pow(v) == 0)
                .collect(),
        }
    }

    //return the terms where every var in vars is present
    pub fn has_all_vars_parts(self, vars: Vec<&Variable>) -> MultiPolynomial<R> {
        MultiPolynomial {
            terms: self
                .terms
                .into_iter()
                .filter(|Term { coeff, monomial }| vars.iter().all(|v| monomial.get_var_pow(v) > 0))
                .collect(),
        }
    }
}

impl<RS: RingStructure + DisplayableStructure> DisplayableStructure
    for MultiPolynomialStructure<RS>
{
    fn elem_to_string(&self, p: &Self::Set) -> String {
        if p.terms.len() == 0 {
            "0".into()
        } else {
            let mut s = String::new();
            for (idx, term) in p.terms.iter().enumerate() {
                if idx != 0 {
                    s += "+";
                }
                s += "(";
                s += &self.coeff_ring.elem_to_string(&term.coeff);
                s += ")";
                s += &term.monomial.to_string();
            }
            s
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MultiPolynomialStructure<RS: RingStructure> {
    coeff_ring: Rc<RS>,
}

impl<RS: RingStructure> MultiPolynomialStructure<RS> {
    pub fn new(coeff_ring: Rc<RS>) -> Self {
        Self { coeff_ring }
    }

    pub fn coeff_ring(&self) -> Rc<RS> {
        self.coeff_ring.clone()
    }
}

impl<RS: RingStructure> Structure for MultiPolynomialStructure<RS> {
    type Set = MultiPolynomial<RS::Set>;
}

impl<RS: RingStructure> EqualityStructure for MultiPolynomialStructure<RS> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        let a = self.reduce(a.clone());
        let b = self.reduce(b.clone());

        let n = a.terms.len();
        if n != b.terms.len() {
            false
        } else {
            (0..n).all(|i| {
                self.coeff_ring.equal(&a.terms[i].coeff, &b.terms[i].coeff)
                    && &a.terms[i].monomial == &b.terms[i].monomial
            })
        }
    }
}

impl<RS: RingStructure> RingStructure for MultiPolynomialStructure<RS> {
    fn zero(&self) -> Self::Set {
        MultiPolynomial { terms: vec![] }
    }

    fn one(&self) -> Self::Set {
        MultiPolynomial {
            terms: vec![Term {
                coeff: self.coeff_ring.one(),
                monomial: Monomial::one(),
            }],
        }
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        MultiPolynomial {
            terms: a
                .terms
                .iter()
                .map(
                    |Term {
                         coeff: c,
                         monomial: m,
                     }| Term {
                        coeff: self.coeff_ring.neg(&c),
                        monomial: m.clone(),
                    },
                )
                .collect(),
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let mut a = a.clone();
        let mut existing_monomials: HashMap<Monomial, usize> = HashMap::new(); //the index of each monomial
        for (
            idx,
            Term {
                coeff: _coeff,
                monomial,
            },
        ) in a.terms.clone().into_iter().enumerate()
        {
            existing_monomials.insert(monomial, idx);
        }
        for Term { coeff, monomial } in &b.terms {
            if existing_monomials.contains_key(&monomial) {
                self.coeff_ring.add_mut(
                    &mut a.terms[*existing_monomials.get(&monomial).unwrap()].coeff,
                    &coeff,
                );
            } else {
                a.terms.push(Term {
                    coeff: coeff.clone(),
                    monomial: monomial.clone(),
                });
            }
        }
        self.reduce(a) //sort the coeffs
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let mut terms: HashMap<Monomial, RS::Set> = HashMap::new();
        for Term {
            coeff: a_coeff,
            monomial: a_monomial,
        } in a.terms.iter()
        {
            for Term {
                coeff: b_coeff,
                monomial: b_monomial,
            } in b.terms.iter()
            {
                let mon = Monomial::mul(a_monomial, b_monomial);
                let coeff = self.coeff_ring.mul(a_coeff, b_coeff);
                self.coeff_ring
                    .add_mut(terms.entry(mon).or_insert(self.coeff_ring.zero()), &coeff);
            }
        }
        self.reduce(MultiPolynomial::new(
            terms
                .into_iter()
                .map(|(monomial, coeff)| Term { coeff, monomial })
                .collect(),
        ))
    }
}

impl<RS: IntegralDomainStructure> IntegralDomainStructure for MultiPolynomialStructure<RS> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut vars = HashSet::new();
        vars.extend(a.free_vars());
        vars.extend(b.free_vars());
        if vars.len() == 0 {
            //a and b are constants
            debug_assert!(a.terms.len() <= 1);
            debug_assert!(b.terms.len() <= 1);
            if b.terms.len() == 0 {
                Err(RingDivisionError::DivideByZero)
            } else if a.terms.len() == 0 {
                Ok(self.zero())
            } else {
                debug_assert!(a.terms.len() == 1);
                debug_assert!(b.terms.len() == 1);
                match self.coeff_ring.div(&a.terms[0].coeff, &b.terms[0].coeff) {
                    Ok(c) => Ok(MultiPolynomial::constant(c)),
                    Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
        } else {
            let var = vars.iter().next().unwrap();
            let a_poly = self.expand(a, var);
            let b_poly = self.expand(b, var);
            let poly_ring = PolynomialStructure::<Self>::new(self.clone().into());
            match poly_ring.div(&a_poly, &b_poly) {
                Ok(c_poly) => Ok(poly_ring.evaluate(&c_poly, &self.var(var.clone()))),
                Err(e) => Err(e),
            }
        }
    }
}

impl<RS: RingStructure> MultiPolynomialStructure<RS> {
    pub fn reduce(&self, p: MultiPolynomial<RS::Set>) -> MultiPolynomial<RS::Set> {
        MultiPolynomial::new(
            p.terms
                .into_iter()
                .filter(|Term { monomial, coeff }| !self.coeff_ring.is_zero(coeff))
                .collect(),
        )
    }

    pub fn var_pow(&self, v: Variable, k: usize) -> MultiPolynomial<RS::Set> {
        MultiPolynomial {
            terms: vec![Term {
                coeff: self.coeff_ring.one(),
                monomial: Monomial::new(vec![VariablePower { var: v, pow: k }]),
            }],
        }
    }

    pub fn var(&self, v: Variable) -> MultiPolynomial<RS::Set> {
        self.var_pow(v, 1)
    }

    pub fn as_constant(&self, p: &MultiPolynomial<RS::Set>) -> Option<RS::Set> {
        if p.terms.len() == 0 {
            Some(self.coeff_ring.zero())
        } else if p.terms.len() == 1 {
            let Term { coeff, monomial } = &p.terms[0];
            if monomial == &Monomial::one() {
                Some(coeff.clone())
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn degree(&self, p: &MultiPolynomial<RS::Set>) -> Option<usize> {
        p.terms.iter().map(|t| t.degree()).max()
    }

    pub fn split_by_degree(
        &self,
        p: MultiPolynomial<RS::Set>,
    ) -> HashMap<usize, MultiPolynomial<RS::Set>> {
        let mut p_by_deg = HashMap::new();
        for term in p.terms {
            // let term = MultiPolynomial::new(vec![term]);
            let deg = term.degree();
            if !p_by_deg.contains_key(&deg) {
                p_by_deg.insert(deg, vec![]);
            }
            p_by_deg.get_mut(&deg).unwrap().push(term);
        }
        p_by_deg
            .into_iter()
            .map(|(d, t)| (d, MultiPolynomial { terms: t }))
            .collect()
    }

    pub fn homogenize(
        &self,
        p: &MultiPolynomial<RS::Set>,
        v: &Variable,
    ) -> MultiPolynomial<RS::Set> {
        if self.is_zero(p) {
            self.zero()
        } else {
            let d = self.degree(p).unwrap();
            let h = MultiPolynomial::new(
                p.terms
                    .iter()
                    .map(|Term { coeff, monomial }| Term {
                        coeff: coeff.clone(),
                        monomial: monomial.homogenize(d, v),
                    })
                    .collect(),
            );
            debug_assert!(h.check_invariants().is_ok());
            h
        }
    }

    pub fn expand(
        &self,
        p: &MultiPolynomial<RS::Set>,
        v: &Variable,
    ) -> Polynomial<MultiPolynomial<RS::Set>> {
        let mut coeffs = vec![];
        for Term { coeff, monomial } in &p.terms {
            let k = monomial.get_var_pow(v);
            while coeffs.len() <= k {
                coeffs.push(self.zero())
            }
            self.add_mut(
                &mut coeffs[k],
                &MultiPolynomial {
                    terms: vec![Term {
                        coeff: coeff.clone(),
                        monomial: Monomial::new(
                            monomial
                                .clone()
                                .prod
                                .into_iter()
                                .filter(|VariablePower { var, pow: _pow }| var != v)
                                .collect(),
                        ),
                    }],
                },
            );
        }
        Polynomial::from_coeffs(coeffs)
    }

    pub fn evaluate(
        &self,
        poly: &MultiPolynomial<RS::Set>,
        values: &HashMap<Variable, impl Borrow<RS::Set>>,
    ) -> RS::Set {
        self.coeff_ring().sum(
            poly.terms
                .iter()
                .map(|Term { coeff, monomial }| {
                    self.coeff_ring()
                        .mul(coeff, &monomial.evaluate(self.coeff_ring(), &values))
                })
                .collect(),
        )
    }
}

impl<R: StructuredType> StructuredType for MultiPolynomial<R>
where
    R::Structure: RingStructure,
{
    type Structure = MultiPolynomialStructure<R::Structure>;

    fn structure() -> Rc<Self::Structure> {
        MultiPolynomialStructure::new(R::structure()).into()
    }
}

impl<R: StructuredType> Display for MultiPolynomial<R>
where
    R::Structure: RingStructure + DisplayableStructure,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Self::structure().elem_to_string(self))
    }
}

impl<R: StructuredType> PartialEq for MultiPolynomial<R>
where
    R::Structure: RingStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: StructuredType> Eq for MultiPolynomial<R> where R::Structure: RingStructure {}

impl<R: StructuredType> MultiPolynomial<R>
where
    R::Structure: RingStructure,
{
    pub fn reduce(self) -> Self {
        Self::structure().reduce(self)
    }

    pub fn var(v: Variable) -> Self {
        Self::structure().var(v)
    }

    pub fn as_constant(&self) -> Option<R> {
        Self::structure().as_constant(self)
    }

    pub fn degree(&self) -> Option<usize> {
        Self::structure().degree(self)
    }

    pub fn homogenize(&self, v: &Variable) -> Self {
        Self::structure().homogenize(self, v)
    }

    pub fn expand(&self, v: &Variable) -> Polynomial<Self> {
        Self::structure().expand(self, v)
    }

    pub fn evaluate(&self, values: &HashMap<Variable, &R>) -> R {
        Self::structure().evaluate(self, values)
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::*;

    #[test]
    fn test_monomial_ordering() {
        let xv = Variable::new(String::from("x"));
        let yv = Variable::new(String::from("y"));

        let xx = Monomial::new(vec![VariablePower {
            var: xv.clone(),
            pow: 2,
        }]);
        let yy = Monomial::new(vec![VariablePower {
            var: yv.clone(),
            pow: 2,
        }]);
        let xy = Monomial::new(vec![
            VariablePower {
                var: xv.clone(),
                pow: 1,
            },
            VariablePower {
                var: yv.clone(),
                pow: 1,
            },
        ]);
        let x = Monomial::new(vec![VariablePower {
            var: xv.clone(),
            pow: 1,
        }]);
        let y = Monomial::new(vec![VariablePower {
            var: yv.clone(),
            pow: 1,
        }]);
        let one = Monomial::one();

        {
            let terms = vec![
                xx.clone(),
                xy.clone(),
                x.clone(),
                yy.clone(),
                y.clone(),
                one.clone(),
            ];
            let mut sorted_terms = terms.clone();
            sorted_terms.sort_by(|a, b| Monomial::lexicographic_order(a, b));

            assert_eq!(terms, sorted_terms);
        }

        {
            let terms = vec![
                xx.clone(),
                xy.clone(),
                yy.clone(),
                x.clone(),
                y.clone(),
                one.clone(),
            ];
            let mut sorted_terms = terms.clone();
            sorted_terms.sort_by(|a, b| Monomial::graded_lexicographic_order(a, b));

            assert_eq!(terms, sorted_terms);
        }
    }

    #[test]
    fn test_reduction() {
        let x = MultiPolynomial::<Integer>::var(Variable::new(String::from("x")));
        let f = MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg(&x)]);
        assert_eq!(f.terms.len(), 0);

        let x = MultiPolynomial::<Integer>::var(Variable::new(String::from("x")));
        let y = MultiPolynomial::<Integer>::var(Variable::new(String::from("y")));
        let g = MultiPolynomial::sum(vec![&x, &y]);
        let h = MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg(&y)]);
        let f = MultiPolynomial::product(vec![&g, &h]);
        println!("g = {}", g);
        println!("h = {}", h);
        println!("f = {}", f);
        assert_eq!(f.terms.len(), 2);
    }

    #[test]
    fn test_division() {
        let x = MultiPolynomial::<Integer>::var(Variable::new(String::from("x")));
        let y = MultiPolynomial::<Integer>::var(Variable::new(String::from("y")));

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(&MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg(&y)]);
        match MultiPolynomial::div(&f, &g) {
            Ok(h) => {
                assert_eq!(f, MultiPolynomial::mul(&g, &h));
            }
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => panic!(),
        }

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(&MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::zero();
        match MultiPolynomial::div(&f, &g) {
            Ok(_) => panic!(),
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
        }

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(&MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::sum(vec![&x]);
        match MultiPolynomial::div(&f, &g) {
            Ok(_) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(RingDivisionError::DivideByZero) => panic!(),
        }
    }

    #[test]
    fn test_elems() {
        let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ring();
        let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ring();
        let z = &MultiPolynomial::<Integer>::var(Variable::new("z")).into_ring();

        let f = x + y + z;
        let g = x - y + z;

        let h = (&f * &g) / &f;
        h.into_set().check_invariants().unwrap();

        println!("f = {}", f);
        println!("g = {}", g);
        println!("fg = {}", &f * &g);
        println!("fg/f = {}", (&f * &g) / &f);

        assert_eq!((&f * &g) / &f, g);
    }
}
