use super::polynomial::*;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt::Display;
use std::hash::Hash;
use std::rc::Rc;
use std::sync::atomic::AtomicUsize;

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
                         monomial: Monomial { prod, .. },
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

impl<RS: RingStructure> Structure for MultiPolynomialStructure<RS> {}

impl<RS: RingStructure + ToStringStructure> ToStringStructure for MultiPolynomialStructure<RS> {
    fn to_string(&self, p: &Self::Set) -> String {
        if p.terms.len() == 0 {
            "0".into()
        } else {
            let mut s = String::new();
            for (idx, term) in p.terms.iter().enumerate() {
                if idx != 0 {
                    s += "+";
                }
                if !self.coeff_ring.equal(&term.coeff, &self.coeff_ring.one()) {
                    s += "(";
                    s += &self.coeff_ring.to_string(&term.coeff);
                    s += ")";
                }
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

impl<RS: RingStructure> SetStructure for MultiPolynomialStructure<RS> {
    type Set = MultiPolynomial<RS::Set>;
}

impl<RS: RingStructure> PartialEqStructure for MultiPolynomialStructure<RS> {
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
impl<RS: RingStructure> EqStructure for MultiPolynomialStructure<RS> {}

impl<RS: RingStructure> SemiRingStructure for MultiPolynomialStructure<RS> {
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

impl<RS: RingStructure> RingStructure for MultiPolynomialStructure<RS> {
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
}

impl<RS: IntegralDomainStructure> UnitsStructure for MultiPolynomialStructure<RS> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div(&self.one(), a)
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

impl<RS: FavoriteAssociateStructure> FavoriteAssociateStructure for MultiPolynomialStructure<RS> {
    fn factor_fav_assoc(&self, mpoly: &Self::Set) -> (Self::Set, Self::Set) {
        match mpoly.terms.first() {
            None => {
                debug_assert!(self.is_zero(mpoly));
                (self.one(), self.zero())
            }
            Some(first_term) => {
                debug_assert!(!self.is_zero(mpoly));
                let (unit, _) = self.coeff_ring().factor_fav_assoc(&first_term.coeff);
                let unit_inv = MultiPolynomial::constant(self.coeff_ring().inv(&unit).unwrap());
                let unit = MultiPolynomial::constant(unit);
                (unit, self.mul(&unit_inv, mpoly))
            }
        }
    }
}

impl<RS: CharZeroStructure> CharZeroStructure for MultiPolynomialStructure<RS> {}

impl<RS: FiniteUnitsStructure> FiniteUnitsStructure for MultiPolynomialStructure<RS> {
    fn all_units(&self) -> Vec<Self::Set> {
        self.coeff_ring()
            .all_units()
            .into_iter()
            .map(|u| MultiPolynomial::constant(u))
            .collect()
    }
}

impl<RS: GreatestCommonDivisorStructure> GreatestCommonDivisorStructure
    for PolynomialStructure<MultiPolynomialStructure<RS>>
{
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl<RS: GreatestCommonDivisorStructure> GreatestCommonDivisorStructure
    for MultiPolynomialStructure<RS>
where
    PolynomialStructure<MultiPolynomialStructure<RS>>:
        SetStructure<Set = Polynomial<MultiPolynomial<RS::Set>>>,
{
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        match x.free_vars().into_iter().chain(y.free_vars()).next() {
            Some(free_var) => {
                let poly_over_self = PolynomialStructure::new(self.clone().into());
                let x_poly = self.expand(x, &free_var);
                let y_poly = self.expand(y, &free_var);
                let g_poly = poly_over_self.gcd(&x_poly, &y_poly);
                poly_over_self.evaluate(&g_poly, &self.var(free_var))
            }
            None => {
                let x = self.as_constant(x).unwrap();
                let y = self.as_constant(y).unwrap();
                let g = self.coeff_ring().gcd(&x, &y);
                MultiPolynomial::constant(g)
            }
        }
    }
}

impl<RS: UniqueFactorizationStructure> UniqueFactorizationStructure
    for MultiPolynomialStructure<RS>
{
}

impl<
    RS: UniqueFactorizationStructure
        + GreatestCommonDivisorStructure
        + CharZeroStructure
        + FiniteUnitsStructure
        + 'static,
> PolynomialStructure<MultiPolynomialStructure<RS>>
where
    PolynomialStructure<MultiPolynomialStructure<RS>>:
        SetStructure<Set = Polynomial<MultiPolynomial<RS::Set>>> + UniqueFactorizationStructure,
    PolynomialStructure<RS>: SetStructure<Set = Polynomial<RS::Set>> + UniqueFactorizationStructure,
    MultiPolynomialStructure<RS>: SetStructure<Set = MultiPolynomial<RS::Set>>
        + UniqueFactorizationStructure
        + GreatestCommonDivisorStructure,
{
    pub fn factor_by_yuns_and_kroneckers_inductively(
        &self,
        factor_poly: impl Fn(&Polynomial<RS::Set>) -> Option<Factored<PolynomialStructure<RS>>>,
        factor_multipoly_coeff: impl Fn(
            &MultiPolynomial<RS::Set>,
        ) -> Option<Factored<MultiPolynomialStructure<RS>>>,
        mpoly: &<Self as SetStructure>::Set,
    ) -> Option<Factored<PolynomialStructure<MultiPolynomialStructure<RS>>>> {
        match |mpoly: &<Self as SetStructure>::Set| -> Option<Polynomial<RS::Set>> {
            let mut const_coeffs = vec![];
            for coeff in mpoly.coeffs() {
                const_coeffs.push(self.coeff_ring().as_constant(&coeff)?);
            }
            Some(Polynomial::from_coeffs(const_coeffs))
        }(mpoly)
        {
            // It is a polynomial with multipolynomial coefficients where all coefficients are constant
            // So we can defer to a univariate factoring algorithm
            Some(poly) => {
                let (unit, factors) = factor_poly(&poly)?.unit_and_factors();
                Some(Factored::new_unchecked(
                    self.clone().into(),
                    unit.apply_map_into(|c| MultiPolynomial::constant(c)),
                    factors
                        .into_iter()
                        .map(|(factor, power)| {
                            (
                                factor.apply_map_into(|c| MultiPolynomial::constant(c)),
                                power,
                            )
                        })
                        .collect(),
                ))
            }
            None => {
                self.factorize_by_yuns_and_kroneckers_method(mpoly, |c| factor_multipoly_coeff(c))
            }
        }
    }
}

impl<
    RS: UniqueFactorizationStructure
        + GreatestCommonDivisorStructure
        + CharZeroStructure
        + FiniteUnitsStructure
        + 'static,
> MultiPolynomialStructure<RS>
where
    MultiPolynomialStructure<RS>:
        SetStructure<Set = MultiPolynomial<RS::Set>> + UniqueFactorizationStructure,
    PolynomialStructure<RS>: SetStructure<Set = Polynomial<RS::Set>> + UniqueFactorizationStructure,
    PolynomialStructure<MultiPolynomialStructure<RS>>:
        SetStructure<Set = Polynomial<MultiPolynomial<RS::Set>>> + UniqueFactorizationStructure,
{
    pub fn factor_by_yuns_and_kroneckers_inductively(
        &self,
        factor_coeff: Rc<dyn Fn(&RS::Set) -> Option<Factored<RS>>>,
        factor_poly: Rc<dyn Fn(&Polynomial<RS::Set>) -> Option<Factored<PolynomialStructure<RS>>>>,
        mpoly: &<Self as SetStructure>::Set,
    ) -> Option<Factored<MultiPolynomialStructure<RS>>> {
        if self.is_zero(mpoly) {
            None
        } else {
            match mpoly.free_vars().into_iter().next() {
                Some(free_var) => {
                    match mpoly.homogeneous_of_degree() {
                        HomogeneousOfDegreeResult::Homogeneous(_) => {
                            // If mpoly is homogeneous we can eliminate a variable right away
                            let dehom_mpoly = self.partial_evaluate(
                                mpoly,
                                HashMap::from([(free_var.clone(), self.coeff_ring().one())]),
                            );
                            let (unit, factors) = self
                                .factor_by_yuns_and_kroneckers_inductively(
                                    factor_coeff.clone(),
                                    factor_poly.clone(),
                                    &dehom_mpoly,
                                )
                                .unwrap()
                                .unit_and_factors();
                            Some(Factored::new_unchecked(
                                self.clone().into(),
                                self.homogenize(&unit, &free_var),
                                factors
                                    .into_iter()
                                    .map(|(factor, power)| {
                                        (self.homogenize(&factor, &free_var), power)
                                    })
                                    .collect(),
                            ))
                        }
                        HomogeneousOfDegreeResult::No => {
                            // Not homogeneous but
                            // There exists a free variable
                            // So turn ourself into a polynomial with respect to that free variable
                            let poly_over_self = PolynomialStructure::new(self.clone().into());
                            let expanded_poly = self.expand(mpoly, &free_var);
                            let free_var = self.var(free_var);
                            let (unit, factors) = poly_over_self
                                .factor_by_yuns_and_kroneckers_inductively(
                                    factor_poly.as_ref(),
                                    |c| {
                                        self.factor_by_yuns_and_kroneckers_inductively(
                                            factor_coeff.clone(),
                                            factor_poly.clone(),
                                            c,
                                        )
                                    },
                                    &expanded_poly,
                                )
                                .unwrap()
                                .unit_and_factors();
                            Some(Factored::new_unchecked(
                                self.clone().into(),
                                poly_over_self.evaluate(&unit, &free_var),
                                factors
                                    .into_iter()
                                    .map(|(factor, power)| {
                                        (poly_over_self.evaluate(&factor, &free_var), power)
                                    })
                                    .collect(),
                            ))
                        }
                        HomogeneousOfDegreeResult::Zero => unreachable!(),
                    }
                }
                None => {
                    // Just an element of the coefficient ring
                    let value = self.as_constant(mpoly).unwrap();
                    let factored = factor_coeff(&value)?;
                    let (unit, factors) = factored.unit_and_factors();
                    Some(Factored::new_unchecked(
                        self.clone().into(),
                        MultiPolynomial::constant(unit),
                        factors
                            .into_iter()
                            .map(|(factor, power)| (MultiPolynomial::constant(factor), power))
                            .collect(),
                    ))
                }
            }
        }
    }
}

impl<RS: RingStructure> MultiPolynomialStructure<RS> {
    pub fn reduce(&self, p: MultiPolynomial<RS::Set>) -> MultiPolynomial<RS::Set> {
        MultiPolynomial::new(
            p.terms
                .into_iter()
                .filter(|Term { coeff, .. }| !self.coeff_ring.is_zero(coeff))
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

    pub fn partial_evaluate(
        &self,
        poly: &MultiPolynomial<RS::Set>,
        values: HashMap<Variable, impl Borrow<RS::Set>>,
    ) -> MultiPolynomial<RS::Set> {
        let mut values = values
            .into_iter()
            .map(|(v, a)| (v, MultiPolynomial::constant(a.borrow().clone())))
            .collect::<HashMap<_, _>>();

        let free_vars = poly.free_vars();
        for free_var in free_vars {
            if !values.contains_key(&free_var) {
                values.insert(free_var.clone(), self.var(free_var));
            }
        }
        MultiPolynomialStructure::new(self.clone().into()).evaluate(
            &poly.apply_map(|x| MultiPolynomial::constant(x.clone())),
            values,
        )
    }

    pub fn evaluate(
        &self,
        poly: &MultiPolynomial<RS::Set>,
        values: HashMap<Variable, impl Borrow<RS::Set>>,
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

impl<R: MetaType> MetaType for MultiPolynomial<R>
where
    R::Structure: RingStructure,
{
    type Structure = MultiPolynomialStructure<R::Structure>;

    fn structure() -> std::rc::Rc<Self::Structure> {
        MultiPolynomialStructure::new(R::structure()).into()
    }
}

impl<R: MetaType> Display for MultiPolynomial<R>
where
    R::Structure: RingStructure + ToStringStructure,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Self::structure().to_string(self))
    }
}

impl<R: MetaType> PartialEq for MultiPolynomial<R>
where
    R::Structure: RingStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: MetaType> Eq for MultiPolynomial<R> where R::Structure: RingStructure {}

impl<R: MetaType> MultiPolynomial<R>
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

    pub fn evaluate(&self, values: HashMap<Variable, impl Borrow<R>>) -> R {
        Self::structure().evaluate(self, values)
    }
}

#[cfg(test)]
mod tests {
    use crate::structure::IntoErgonomic;

    use algebraeon_nzq::*;

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
        let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
        let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();
        let z = &MultiPolynomial::<Integer>::var(Variable::new("z")).into_ergonomic();

        let f = x + y + z;
        let g = x - y + z;

        let h = (&f * &g) / &f;
        h.into_verbose().check_invariants().unwrap();

        println!("f = {}", f);
        println!("g = {}", g);
        println!("fg = {}", &f * &g);
        println!("fg/f = {}", (&f * &g) / &f);

        assert_eq!((&f * &g) / &f, g);
    }

    // #[test]
    // fn test_gcd_and_factor() {
    //     let x = &MultiPolynomial::<Integer>::var(Variable::new("x")).into_ergonomic();
    //     let y = &MultiPolynomial::<Integer>::var(Variable::new("y")).into_ergonomic();

    //     let a = y - x;
    //     let b = y.pow(2) - x.pow(3);
    //     let c = y.pow(3) * x + 5 + y;
    //     let d = x.pow(5) * y.pow(4) - 1;
    //     let e = y.pow(4) - x.pow(4);
    //     let f = 24 * (y - x);

    //     assert!(a.clone().into_verbose().is_irreducible());
    //     assert!(b.clone().into_verbose().is_irreducible());
    //     assert!(c.clone().into_verbose().is_irreducible());
    //     assert!(d.clone().into_verbose().is_irreducible());
    //     assert!(!e.clone().into_verbose().is_irreducible());
    //     assert!(!f.clone().into_verbose().is_irreducible());

    //     assert!(MultiPolynomial::are_associate(
    //         &b.clone().into_verbose(),
    //         &MultiPolynomial::gcd(&(&a * &b).into_verbose(), &(&b * &c).into_verbose())
    //     ));

    //     assert!(MultiPolynomial::are_associate(
    //         &(&b * &c).into_verbose(),
    //         &MultiPolynomial::gcd(
    //             &(&a * &b * &c).into_verbose(),
    //             &(&b * &c * &d).into_verbose()
    //         )
    //     ));
    // }
}
