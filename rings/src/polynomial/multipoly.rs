use crate::structure::*;
use algebraeon_nzq::*;
use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use std::sync::atomic::AtomicUsize;

#[derive(Debug, Clone)]
pub struct Variable {
    pub(crate) ident: usize,
    pub(crate) name: String,
}

impl PartialEq for Variable {
    fn eq(&self, other: &Self) -> bool {
        self.ident == other.ident
    }
}

impl Eq for Variable {}

impl Hash for Variable {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.ident.hash(state);
    }
}

impl Variable {
    pub fn new<S: Into<String>>(name: S) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let name = name.into();
        assert!(!name.is_empty());
        Self {
            ident: COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed),
            name,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct VariablePower {
    pub(crate) var: Variable,
    pub(crate) pow: usize,
}

#[derive(Debug, Clone)]
pub struct Monomial {
    pub(crate) prod: Vec<VariablePower>, //should be sorted by variable ident
    pub(crate) ident_lookup: HashMap<usize, usize>, //point from variable ident to index of that variables power in self.prod
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
        if self.prod.is_empty() {
            write!(f, "1")
        } else {
            for VariablePower { var, pow } in &self.prod {
                write!(f, "{}", &var.name)?;
                if *pow != 1 {
                    write!(f, "^")?;
                    write!(f, "{}", pow)?;
                }
            }
            Ok(())
        }
    }
}

impl Monomial {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
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

    pub fn new(mut prod: Vec<VariablePower>) -> Self {
        prod.retain(|VariablePower { var: _var, pow }| pow != &0);
        prod.sort_by_key(|vpow| vpow.var.ident);
        let mut ident_lookup = HashMap::new();
        for (idx, VariablePower { var, pow: _pow }) in prod.iter().enumerate() {
            ident_lookup.insert(var.ident, idx);
        }
        Self { prod, ident_lookup }
    }

    pub fn one() -> Self {
        Monomial {
            prod: vec![],
            ident_lookup: HashMap::new(),
        }
    }

    pub fn degree(&self) -> usize {
        let mut d = 0;
        for VariablePower { var: _var, pow } in &self.prod {
            d += pow;
        }
        d
    }

    pub fn homogenize(&self, target_degree: usize, v: &Variable) -> Self {
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

    pub fn get_var_pow(&self, v: &Variable) -> usize {
        if self.ident_lookup.contains_key(&v.ident) {
            self.prod[*self.ident_lookup.get(&v.ident).unwrap()].pow
        } else {
            0
        }
    }

    pub fn free_vars(&self) -> HashSet<Variable> {
        self.prod
            .iter()
            .map(|VariablePower { var, pow: _pow }| var.clone())
            .collect()
    }

    pub fn evaluate<RS: RingSignature>(
        &self,
        ring: &RS,
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

    pub fn mul(a: &Self, b: &Self) -> Self {
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

    pub fn lexicographic_order(a: &Self, b: &Self) -> std::cmp::Ordering {
        let mut i = 0;
        while i < std::cmp::min(a.prod.len(), b.prod.len()) {
            #[allow(clippy::comparison_chain)]
            if a.prod[i].var.ident < b.prod[i].var.ident {
                return std::cmp::Ordering::Less;
            } else if a.prod[i].var.ident > b.prod[i].var.ident {
                return std::cmp::Ordering::Greater;
            }
            if a.prod[i].pow > b.prod[i].pow {
                return std::cmp::Ordering::Less;
            }
            if a.prod[i].pow < b.prod[i].pow {
                return std::cmp::Ordering::Greater;
            }
            i += 1;
        }
        #[allow(clippy::comparison_chain)]
        return if a.prod.len() > b.prod.len() {
            std::cmp::Ordering::Less
        } else if a.prod.len() < b.prod.len() {
            std::cmp::Ordering::Greater
        } else {
            std::cmp::Ordering::Equal
        };
    }

    pub fn graded_lexicographic_order(a: &Self, b: &Self) -> std::cmp::Ordering {
        #[allow(clippy::comparison_chain)]
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
    pub(crate) coeff: ElemT,
    pub(crate) monomial: Monomial,
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
    pub(crate) terms: Vec<Term<R>>, //sorted by monomial ordering
}

impl<R: Clone> MultiPolynomial<R> {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
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

    pub fn new(mut terms: Vec<Term<R>>) -> Self {
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
                         monomial: Monomial { prod, .. },
                     }| Term {
                        coeff,
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
                .filter(|Term { monomial, .. }| monomial.get_var_pow(v) == 0)
                .collect(),
        }
    }

    //return the terms where every var in vars is present
    pub fn has_all_vars_parts(self, vars: Vec<&Variable>) -> MultiPolynomial<R> {
        MultiPolynomial {
            terms: self
                .terms
                .into_iter()
                .filter(|Term { monomial, .. }| vars.iter().all(|v| monomial.get_var_pow(v) > 0))
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HomogeneousOfDegreeResult {
    No,
    Zero,
    Homogeneous(usize),
}

impl<R: Clone> MultiPolynomial<R> {
    pub fn homogeneous_of_degree(&self) -> HomogeneousOfDegreeResult {
        let mut terms = self.terms.iter().collect::<Vec<_>>();
        match terms.pop() {
            Some(first_term) => {
                let d = first_term.degree();
                for term in terms {
                    if d != term.degree() {
                        return HomogeneousOfDegreeResult::No;
                    }
                }
                HomogeneousOfDegreeResult::Homogeneous(d)
            }
            None => {
                // Zero polynomial is homogeneous of degree infinity
                HomogeneousOfDegreeResult::Zero
            }
        }
    }
}
