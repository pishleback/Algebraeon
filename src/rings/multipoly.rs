#![allow(dead_code)]

use std::collections::HashMap;
use std::sync::atomic::AtomicUsize;

use itertools::Itertools;
use malachite_nz::natural::Natural;

use super::nzq::*;
use super::poly::*;
use super::ring::*;

pub const ZZ_MULTIPOLY: MultiPolynomialRing<IntegerRing> = MultiPolynomialRing { ring: &ZZ };
pub const QQ_MULTIPOLY: MultiPolynomialRing<RationalField> = MultiPolynomialRing { ring: &QQ };

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

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
struct VariablePower {
    var: Variable,
    pow: Natural,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
struct Monomial {
    prod: Vec<VariablePower>, //should be sorted by variable ident
}

impl ToString for Monomial {
    fn to_string(&self) -> String {
        if self.prod.len() == 0 {
            String::from("1")
        } else {
            let mut ans = String::from("");
            for VariablePower { var, pow } in &self.prod {
                ans += &var.name;
                ans += "^";
                ans += pow.to_string().as_str();
            }
            ans
        }
    }
}

impl Monomial {
    fn new(mut prod: Vec<VariablePower>) -> Self {
        prod.sort_by_key(|vpow| vpow.var.ident);
        Self { prod }
    }

    fn mul(a: &Self, b: &Self) -> Self {
        Self {
            prod: {
                let mut prod = HashMap::new();
                for VariablePower { var: v, pow: k } in &a.prod {
                    *prod.entry(v.clone()).or_insert(Natural::from(0u8)) += k;
                }
                for VariablePower { var: v, pow: k } in &b.prod {
                    *prod.entry(v.clone()).or_insert(Natural::from(0u8)) += k;
                }
                prod.into_iter()
                    .map(|(v, k)| VariablePower { var: v, pow: k })
                    .collect()
            },
        }
    }

    fn lexicographic_order(a: &Self, b: &Self) -> std::cmp::Ordering {
        let mut i = 0;
        while i < std::cmp::min(a.prod.len(), b.prod.len()) {
            if a.prod[i].var.ident < b.prod[i].var.ident {
                return std::cmp::Ordering::Greater;
            } else if a.prod[i].var.ident > b.prod[i].var.ident {
                return std::cmp::Ordering::Greater;
            } else {
                if a.prod[i].pow > b.prod[i].pow {
                    return std::cmp::Ordering::Greater;
                } else if a.prod[i].pow < b.prod[i].pow {
                    return std::cmp::Ordering::Greater;
                } else {
                    i += 1;
                }
            }
        }
        if a.prod.len() > b.prod.len() {
            return std::cmp::Ordering::Greater;
        } else if a.prod.len() < b.prod.len() {
            return std::cmp::Ordering::Less;
        } else {
            return std::cmp::Ordering::Equal;
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Term<ElemT: Clone + PartialEq + Eq> {
    coeff: ElemT,
    monomial: Monomial,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct MultiPolynomial<ElemT: Clone + PartialEq + Eq> {
    terms: Vec<Term<ElemT>>, //sorted by monomial ordering
}

impl<ElemT: Clone + PartialEq + Eq> MultiPolynomial<ElemT> {
    fn new(mut terms: Vec<Term<ElemT>>) -> Self {
        terms.sort_by(|t1, t2| Monomial::lexicographic_order(&t1.monomial, &t2.monomial));
        Self { terms }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct MultiPolynomialRing<'a, R: ComRing> {
    ring: &'a R,
}

impl<'a, R: ComRing> MultiPolynomialRing<'a, R> {}

impl<'a, R: ComRing> ComRing for MultiPolynomialRing<'a, R> {
    type ElemT = MultiPolynomial<R::ElemT>;

    fn to_string(&self, elem: &Self::ElemT) -> String {
        if elem.terms.len() == 0 {
            String::from("0")
        } else {
            let mut ans = String::from("");
            for (idx, term) in elem.terms.iter().enumerate() {
                if idx != 0 {
                    ans += "+";
                }
                ans += "(";
                ans += self.ring.to_string(&term.coeff).as_str();
                ans += ")";
                ans += term.monomial.to_string().as_str();
            }
            ans
        }
    }

    fn zero(&self) -> Self::ElemT {
        MultiPolynomial { terms: vec![] }
    }

    fn one(&self) -> Self::ElemT {
        MultiPolynomial {
            terms: vec![Term {
                coeff: self.ring.one(),
                monomial: Monomial { prod: vec![] },
            }],
        }
    }

    fn neg_mut(&self, elem: &mut Self::ElemT) {
        for Term { coeff, monomial } in &mut elem.terms {
            self.ring.neg_mut(coeff);
        }
    }

    fn add_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT) {
        elem.clone_from(&self.add(elem.clone(), offset.clone()))
    }

    fn add(&self, mut elem: Self::ElemT, offset: Self::ElemT) -> Self::ElemT {
        let mut existing_monomials: HashMap<Monomial, usize> = HashMap::new(); //the index of each monomial
        for (idx, Term { coeff, monomial }) in elem.terms.clone().into_iter().enumerate() {
            existing_monomials.insert(monomial, idx);
        }
        for Term { coeff, monomial } in offset.terms {
            if existing_monomials.contains_key(&monomial) {
                self.ring.add_mut(
                    &mut elem.terms[*existing_monomials.get(&monomial).unwrap()].coeff,
                    &coeff,
                );
            } else {
                elem.terms.push(Term { coeff, monomial });
            }
        }
        MultiPolynomial::new(elem.terms) //sort the coeffs
    }

    fn mul_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT) {
        elem.clone_from(&self.mul_refs(&elem, &offset))
    }

    fn mul_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        let mut ans = self.zero();
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
                self.add_mut(
                    &mut ans,
                    &self.term(Term {
                        coeff: self.ring.mul_refs(a_coeff, b_coeff),
                        monomial: Monomial::mul(a_monomial, b_monomial),
                    }),
                );
            }
        }
        ans
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        todo!()
    }
}

impl<'a, R: ComRing> MultiPolynomialRing<'a, R> {
    pub fn var_pow(&self, v: Variable, k: Natural) -> MultiPolynomial<R::ElemT> {
        MultiPolynomial {
            terms: vec![Term {
                coeff: self.ring.one(),
                monomial: Monomial {
                    prod: vec![VariablePower { var: v, pow: k }],
                },
            }],
        }
    }

    pub fn var(&self, v: Variable) -> MultiPolynomial<R::ElemT> {
        self.var_pow(v, Natural::from(1u8))
    }

    pub fn term(&self, t: Term<R::ElemT>) -> MultiPolynomial<R::ElemT> {
        MultiPolynomial { terms: vec![t] }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn idk() {
        let x = Variable::new(String::from("x"));
        let y = Variable::new(String::from("y"));
        let z = Variable::new(String::from("z"));

        let a = ZZ_MULTIPOLY.var_pow(x, Natural::from(2u8));
        let b = ZZ_MULTIPOLY.var_pow(y, Natural::from(3u8));
        let c = ZZ_MULTIPOLY.var_pow(z, Natural::from(4u8));

        let g = ZZ_MULTIPOLY.sum(vec![ZZ_MULTIPOLY.product(vec![a, b]), c]);

        println!("{}", ZZ_MULTIPOLY.to_string(&g));
    }
}
