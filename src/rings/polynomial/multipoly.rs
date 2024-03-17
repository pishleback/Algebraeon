use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use std::sync::atomic::AtomicUsize;

use super::super::numbers::nzq::*;
use super::super::polynomial::poly::*;
use super::super::ring::*;

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
    pub fn new(name: String) -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
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
}

#[derive(Debug, Clone)]
pub struct MultiPolynomial<Ring: ComRing> {
    terms: Vec<Term<Ring>>, //sorted by monomial ordering
}

impl<Ring: ComRing> PartialEq for MultiPolynomial<Ring> {
    fn eq(&self, other: &Self) -> bool {
        let n = self.terms.len();

        if n != other.terms.len() {
            false
        } else {
            (0..n).all(|i| {
                self.terms[i].coeff == other.terms[i].coeff
                    && self.terms[i].monomial == other.terms[i].monomial
            })
        }
    }
}

impl<Ring: ComRing> Eq for MultiPolynomial<Ring> {}

impl<Ring: ComRing + Display> Display for MultiPolynomial<Ring> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.terms.len() == 0 {
            write!(f, "0")
        } else {
            for (idx, term) in self.terms.iter().enumerate() {
                if idx != 0 {
                    write!(f, "+");
                }
                write!(f, "(");
                write!(f, "{}", term.coeff);
                write!(f, ")");
                write!(f, "{}", term.monomial);
            }
            Ok(())
        }
    }
}

impl<Ring: ComRing> ComRing for MultiPolynomial<Ring> {
    fn zero() -> Self {
        MultiPolynomial { terms: vec![] }
    }

    fn one() -> Self {
        MultiPolynomial {
            terms: vec![Term {
                coeff: Ring::one(),
                monomial: Monomial::one(),
            }],
        }
    }

    fn neg_mut(&mut self) {
        for Term {
            coeff,
            monomial: _monomial,
        } in &mut self.terms
        {
            Ring::neg_mut(coeff);
        }
    }

    fn add_mut(&mut self, offset: &Self) {
        self.clone_from(&Self::add(self.clone(), offset.clone()))
    }

    fn add(mut elem: Self, offset: Self) -> Self {
        let mut existing_monomials: HashMap<Monomial, usize> = HashMap::new(); //the index of each monomial
        for (
            idx,
            Term {
                coeff: _coeff,
                monomial,
            },
        ) in elem.terms.clone().into_iter().enumerate()
        {
            existing_monomials.insert(monomial, idx);
        }
        for Term { coeff, monomial } in offset.terms {
            if existing_monomials.contains_key(&monomial) {
                Ring::add_mut(
                    &mut elem.terms[*existing_monomials.get(&monomial).unwrap()].coeff,
                    &coeff,
                );
            } else {
                elem.terms.push(Term { coeff, monomial });
            }
        }
        MultiPolynomial::new(
            elem.terms
                .into_iter()
                .filter(|term| term.coeff != Ring::zero())
                .collect(),
        ) //sort the coeffs
    }

    fn mul_mut(&mut self, offset: &Self) {
        self.clone_from(&Self::mul_refs(&self, &offset))
    }

    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut terms: HashMap<Monomial, Ring> = HashMap::new();
        for Term {
            coeff: a_coeff,
            monomial: a_monomial,
        } in &a.terms
        {
            for Term {
                coeff: b_coeff,
                monomial: b_monomial,
            } in &b.terms
            {
                let mon = Monomial::mul(a_monomial, b_monomial);
                let coeff = Ring::mul_refs(a_coeff, b_coeff);
                Ring::add_mut(terms.entry(mon).or_insert(Ring::zero()), &coeff);
            }
        }
        MultiPolynomial::new(
            terms
                .into_iter()
                .filter(|(_monomial, coeff)| coeff != &Ring::zero())
                .map(|(monomial, coeff)| Term { coeff, monomial })
                .collect(),
        )
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
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
                Ok(Self::zero())
            } else {
                debug_assert!(a.terms.len() == 1);
                debug_assert!(b.terms.len() == 1);
                match Ring::div_refs(&a.terms[0].coeff, &b.terms[0].coeff) {
                    Ok(c) => Ok(Self::constant(c)),
                    Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
        } else {
            let var = vars.iter().next().unwrap();
            let a_poly = a.expand(var);
            let b_poly = b.expand(var);
            match Polynomial::div(a_poly, b_poly) {
                Ok(c_poly) => Ok(c_poly.evaluate(&Self::var(var.clone()))),

                Err(e) => Err(e),
            }
        }
    }
}

impl<Ring: IntegralDomain> IntegralDomain for MultiPolynomial<Ring> {}

impl<Ring: ComRing> MultiPolynomial<Ring> {
    fn new(mut terms: Vec<Term<Ring>>) -> Self {
        terms.sort_by(|t1, t2| Monomial::lexicographic_order(&t1.monomial, &t2.monomial));
        Self { terms }
    }

    fn check_invariants(&self, poly: MultiPolynomial<Ring>) -> Result<(), &'static str> {
        for term in &poly.terms {
            match term.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            }
            if term.coeff == Ring::zero() {
                return Err("coeff should not be zero");
            }
        }

        if !(0..poly.terms.len() - 1).all(|i| {
            Monomial::lexicographic_order(&poly.terms[i].monomial, &poly.terms[i + 1].monomial)
                .is_le()
        }) {
            return Err("terms are not sorted");
        }

        Ok(())
    }

    pub fn var_pow(v: Variable, k: usize) -> MultiPolynomial<Ring> {
        MultiPolynomial {
            terms: vec![Term {
                coeff: Ring::one(),
                monomial: Monomial::new(vec![VariablePower { var: v, pow: k }]),
            }],
        }
    }

    pub fn var(v: Variable) -> MultiPolynomial<Ring> {
        Self::var_pow(v, 1)
    }

    pub fn constant(c: Ring) -> MultiPolynomial<Ring> {
        if c == Ring::zero() {
            Self::zero()
        } else {
            MultiPolynomial {
                terms: vec![Term {
                    coeff: c,
                    monomial: Monomial::one(),
                }],
            }
        }
    }

    pub fn as_constant(&self) -> Option<Ring> {
        if self.terms.len() == 0 {
            Some(Ring::zero())
        } else if self.terms.len() == 1 {
            let Term { coeff, monomial } = &self.terms[0];
            if monomial == &Monomial::one() {
                Some(coeff.clone())
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn degree(&self) -> Option<usize> {
        if self.terms.len() == 0 {
            None
        } else {
            let mut d = 0;
            for Term {
                coeff: _coeff,
                monomial,
            } in &self.terms
            {
                d = std::cmp::max(d, monomial.degree())
            }
            Some(d)
        }
    }

    pub fn term(t: Term<Ring>) -> MultiPolynomial<Ring> {
        MultiPolynomial { terms: vec![t] }
    }

    pub fn free_vars(&self) -> HashSet<Variable> {
        let mut vars = HashSet::new();
        for term in &self.terms {
            vars.extend(term.monomial.free_vars());
        }
        vars
    }

    pub fn homogenize(&self, v: &Variable) -> MultiPolynomial<Ring> {
        if self == &Self::zero() {
            Self::zero()
        } else {
            let d = self.degree().unwrap();
            let h = MultiPolynomial::new(
                self.terms
                    .iter()
                    .map(|Term { coeff, monomial }| Term {
                        coeff: coeff.clone(),
                        monomial: monomial.homogenize(d, v),
                    })
                    .collect(),
            );
            debug_assert!(self.check_invariants(h.clone()).is_ok());
            h
        }
    }

    pub fn expand(&self, v: &Variable) -> Polynomial<MultiPolynomial<Ring>> {
        let mut coeffs = vec![];
        for Term { coeff, monomial } in &self.terms {
            let k = monomial.get_var_pow(v);
            while coeffs.len() <= k {
                coeffs.push(Self::zero())
            }
            Self::add_mut(
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
        let f = MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg(x.clone())]);
        assert_eq!(f.terms.len(), 0);

        let x = MultiPolynomial::<Integer>::var(Variable::new(String::from("x")));
        let y = MultiPolynomial::<Integer>::var(Variable::new(String::from("y")));
        let f = MultiPolynomial::product(vec![
            &MultiPolynomial::sum(vec![&x, &y]),
            &MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg_ref(&y)]),
        ]);
        println!("{}", f.to_string());
        assert_eq!(f.terms.len(), 2);
    }

    #[test]
    fn test_division() {
        let x = MultiPolynomial::<Integer>::var(Variable::new(String::from("x")));
        let y = MultiPolynomial::<Integer>::var(Variable::new(String::from("y")));

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::sum(vec![&x, &MultiPolynomial::neg_ref(&y)]);
        match MultiPolynomial::div_refs(&f, &g) {
            Ok(h) => {
                assert_eq!(f, MultiPolynomial::mul_refs(&g, &h));
            }
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => panic!(),
        }

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::zero();
        match MultiPolynomial::div_refs(&f, &g) {
            Ok(_) => panic!(),
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
        }

        let f = MultiPolynomial::sum(vec![
            &MultiPolynomial::product(vec![&x, &x]),
            &MultiPolynomial::neg(MultiPolynomial::product(vec![&y, &y])),
        ]);
        let g = MultiPolynomial::sum(vec![&x]);
        match MultiPolynomial::div_refs(&f, &g) {
            Ok(_) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(RingDivisionError::DivideByZero) => panic!(),
        }
    }
}
