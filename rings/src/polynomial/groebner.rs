//! Gröbner bases for ideals of multivariate polynomial rings over a field.
//!
//! A *Gröbner basis* of an ideal `I` of a polynomial ring `k[x_1, ..., x_n]`
//! (with respect to a fixed monomial ordering `>`) is a finite generating set
//! `G = {g_1, ..., g_t}` of `I` such that the leading monomials of `G` generate
//! the leading-monomial ideal of the whole ideal:
//! `(lm(g_1), ..., lm(g_t)) = lm(I)`.
//!
//! Gröbner bases turn many otherwise intractable questions about polynomial
//! ideals into algorithmic ones: ideal membership, ideal equality, elimination
//! of variables, intersection of ideals, and computation in the quotient ring
//! `k[x_1, ..., x_n] / I`. This module provides Buchberger's algorithm to
//! compute them, a verification routine, multivariate division, and the FGLM
//! algorithm to convert a Gröbner basis between monomial orderings for
//! zero-dimensional ideals.

use crate::polynomial::Monomial;
use crate::polynomial::MultiPolynomial;
use crate::polynomial::MultiPolynomialStructure;
use crate::polynomial::Term;
use crate::polynomial::Variable;
use crate::polynomial::VariablePower;
use crate::structure::*;
use algebraeon_structures::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;

/// The three classical monomial orderings supported for Gröbner basis
/// computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MonomialOrderingKind {
    /// Lexicographic order. Monomials are compared like words in a dictionary:
    /// the monomial with the larger exponent in the most significant variable
    /// is larger. Lexicographic Gröbner bases are well suited to elimination
    /// because they "triangularise" a system, but are typically the most
    /// expensive to compute.
    Lexicographic,
    /// Graded (degree) lexicographic order. Monomials are compared first by
    /// total degree (higher degree is larger) and ties are broken
    /// lexicographically.
    GradedLexicographic,
    /// Graded reverse lexicographic order. Monomials are compared first by total
    /// degree; ties are broken by looking at the *least* significant variable
    /// and declaring the monomial with the *smaller* exponent there to be the
    /// larger one. This is usually the cheapest ordering for Buchberger's
    /// algorithm and the default choice when no elimination is required.
    GradedReverseLexicographic,
}

/// A monomial ordering: a total order on monomials that is compatible with
/// multiplication (`a > b` implies `ac > bc`) and for which `1` is the smallest
/// monomial. Such an ordering is what makes the notion of "leading term"
/// well defined and is the data on which every Gröbner basis depends.
///
/// The ordering carries an explicit list of variables in *descending*
/// significance (`variables[0]` is the largest variable). Any variable that
/// appears in a compared monomial but is absent from this list is treated as
/// less significant than every listed variable and is ordered by its internal
/// identifier, which keeps the comparison a total order on all monomials.
#[derive(Debug, Clone)]
pub struct MonomialOrdering {
    kind: MonomialOrderingKind,
    variables: Vec<Variable>,
}

impl MonomialOrdering {
    /// Build a monomial ordering of the given kind. `variables` lists the
    /// variables in descending significance.
    pub fn new(kind: MonomialOrderingKind, variables: Vec<Variable>) -> Self {
        Self { kind, variables }
    }

    /// Lexicographic ordering with the given variables in descending
    /// significance.
    pub fn lex(variables: Vec<Variable>) -> Self {
        Self::new(MonomialOrderingKind::Lexicographic, variables)
    }

    /// Graded lexicographic ordering with the given variables in descending
    /// significance.
    pub fn grlex(variables: Vec<Variable>) -> Self {
        Self::new(MonomialOrderingKind::GradedLexicographic, variables)
    }

    /// Graded reverse lexicographic ordering with the given variables in
    /// descending significance.
    pub fn grevlex(variables: Vec<Variable>) -> Self {
        Self::new(MonomialOrderingKind::GradedReverseLexicographic, variables)
    }

    pub fn kind(&self) -> MonomialOrderingKind {
        self.kind
    }

    pub fn variables(&self) -> &[Variable] {
        &self.variables
    }

    /// The full list of variables used to compare `a` and `b`: the declared
    /// variables in order, followed by any further variables occurring in `a`
    /// or `b`, sorted by their internal identifier. This guarantees the
    /// comparison is total even when a monomial uses a variable the ordering was
    /// not explicitly told about.
    fn resolved_variables(&self, a: &Monomial, b: &Monomial) -> Vec<Variable> {
        let mut order = self.variables.clone();
        let known: HashSet<Variable> = order.iter().cloned().collect();
        let mut extra: Vec<Variable> = a
            .free_vars()
            .into_iter()
            .chain(b.free_vars())
            .filter(|v| !known.contains(v))
            .collect();
        extra.sort_by_key(Variable::ident);
        extra.dedup();
        order.extend(extra);
        order
    }

    /// Compare two monomials, returning [`Ordering::Greater`] when `a` is the
    /// larger monomial in this ordering. The leading term of a polynomial is the
    /// term whose monomial is the maximum under this comparison.
    pub fn cmp(&self, a: &Monomial, b: &Monomial) -> Ordering {
        match self.kind {
            MonomialOrderingKind::Lexicographic => {
                for v in self.resolved_variables(a, b) {
                    let c = a.get_var_pow(&v).cmp(&b.get_var_pow(&v));
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                Ordering::Equal
            }
            MonomialOrderingKind::GradedLexicographic => {
                let c = a.degree().cmp(&b.degree());
                if c != Ordering::Equal {
                    return c;
                }
                for v in self.resolved_variables(a, b) {
                    let c = a.get_var_pow(&v).cmp(&b.get_var_pow(&v));
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                Ordering::Equal
            }
            MonomialOrderingKind::GradedReverseLexicographic => {
                let c = a.degree().cmp(&b.degree());
                if c != Ordering::Equal {
                    return c;
                }
                // Tie-break: the monomial with the smaller exponent in the last
                // variable where they differ (reading from least significant) is
                // the larger one.
                for v in self.resolved_variables(a, b).into_iter().rev() {
                    let c = b.get_var_pow(&v).cmp(&a.get_var_pow(&v));
                    if c != Ordering::Equal {
                        return c;
                    }
                }
                Ordering::Equal
            }
        }
    }
}

impl Variable {
    pub(crate) fn ident(&self) -> usize {
        self.ident
    }
}

impl<FS: FieldSignature, FSB: BorrowedStructure<FS>> MultiPolynomialStructure<FS, FSB> {
    /// The leading term of `p` with respect to `ord`: the term whose monomial is
    /// the largest. Returns `None` for the zero polynomial, which has no leading
    /// term. The polynomial is reduced first so that terms with zero
    /// coefficients are ignored.
    pub fn leading_term(
        &self,
        ord: &MonomialOrdering,
        p: &MultiPolynomial<FS::Elem>,
    ) -> Option<Term<FS::Elem>> {
        self.reduce(p.clone())
            .terms
            .into_iter()
            .max_by(|s, t| ord.cmp(&s.monomial, &t.monomial))
    }

    /// The leading monomial of `p`: the monomial of its leading term.
    pub fn leading_monomial(
        &self,
        ord: &MonomialOrdering,
        p: &MultiPolynomial<FS::Elem>,
    ) -> Option<Monomial> {
        self.leading_term(ord, p).map(|t| t.monomial)
    }

    /// The leading coefficient of `p`: the coefficient of its leading term.
    pub fn leading_coefficient(
        &self,
        ord: &MonomialOrdering,
        p: &MultiPolynomial<FS::Elem>,
    ) -> Option<FS::Elem> {
        self.leading_term(ord, p).map(|t| t.coeff)
    }

    /// The single-term polynomial `coeff * monomial`.
    fn term_poly(&self, coeff: FS::Elem, monomial: Monomial) -> MultiPolynomial<FS::Elem> {
        if self.coeff_ring().is_zero(&coeff) {
            self.zero()
        } else {
            MultiPolynomial {
                terms: vec![Term { coeff, monomial }],
            }
        }
    }

    /// Scale `p` so that its leading coefficient is `1`. A monic polynomial is
    /// the canonical representative of its associate class; reduced Gröbner
    /// bases are required to be monic so that they are unique.
    pub fn make_monic(
        &self,
        ord: &MonomialOrdering,
        p: &MultiPolynomial<FS::Elem>,
    ) -> MultiPolynomial<FS::Elem> {
        match self.leading_coefficient(ord, p) {
            None => self.zero(),
            Some(lc) => {
                let inv = self
                    .coeff_ring()
                    .try_divide(&self.coeff_ring().one(), &lc)
                    .unwrap();
                self.mul(&MultiPolynomial::constant(inv), p)
            }
        }
    }

    /// Multivariate division: reduce `f` modulo the ordered list `divisors`,
    /// returning the remainder `r`. The remainder satisfies
    /// `f = q_1 g_1 + ... + q_s g_s + r` where no monomial of `r` is divisible
    /// by the leading monomial of any divisor. At each step the leading term of
    /// the running polynomial is cancelled by whichever divisor has a leading
    /// monomial dividing it; if none does, that leading term is final and is
    /// moved into the remainder.
    ///
    /// When `divisors` is a Gröbner basis the remainder is the unique *normal
    /// form* of `f` modulo the ideal and does not depend on the order of the
    /// divisors; in particular `f` lies in the ideal exactly when the remainder
    /// is zero.
    pub fn reduce_modulo(
        &self,
        ord: &MonomialOrdering,
        f: &MultiPolynomial<FS::Elem>,
        divisors: &[MultiPolynomial<FS::Elem>],
    ) -> MultiPolynomial<FS::Elem> {
        let divisors: Vec<MultiPolynomial<FS::Elem>> = divisors
            .iter()
            .map(|g| self.reduce(g.clone()))
            .filter(|g| !self.is_zero(g))
            .collect();
        let leads: Vec<Term<FS::Elem>> = divisors
            .iter()
            .map(|g| self.leading_term(ord, g).unwrap())
            .collect();

        let mut p = self.reduce(f.clone());
        let mut remainder = self.zero();
        while !self.is_zero(&p) {
            let lt_p = self.leading_term(ord, &p).unwrap();
            let mut divided = false;
            for (g, lead) in divisors.iter().zip(leads.iter()) {
                if let Some(monomial_quotient) = lt_p.monomial.try_div(&lead.monomial) {
                    let coeff_quotient = self
                        .coeff_ring()
                        .try_divide(&lt_p.coeff, &lead.coeff)
                        .unwrap();
                    let factor = self.term_poly(coeff_quotient, monomial_quotient);
                    p = self.sub(&p, &self.mul(&factor, g));
                    divided = true;
                    break;
                }
            }
            if !divided {
                let lt_poly = self.term_poly(lt_p.coeff.clone(), lt_p.monomial.clone());
                remainder = self.add(&remainder, &lt_poly);
                p = self.sub(&p, &lt_poly);
            }
        }
        remainder
    }

    /// The S-polynomial ("syzygy polynomial") of `f` and `g`. Writing
    /// `L = lcm(lm(f), lm(g))`, it is
    /// `S(f, g) = (L / lt(f)) * f - (L / lt(g)) * g`,
    /// the combination engineered so that the leading terms of the two summands
    /// cancel. S-polynomials are the only obstructions to a generating set being
    /// a Gröbner basis: a set is a Gröbner basis precisely when every
    /// S-polynomial reduces to zero modulo it (Buchberger's criterion).
    pub fn s_polynomial(
        &self,
        ord: &MonomialOrdering,
        f: &MultiPolynomial<FS::Elem>,
        g: &MultiPolynomial<FS::Elem>,
    ) -> MultiPolynomial<FS::Elem> {
        let (Some(lt_f), Some(lt_g)) = (self.leading_term(ord, f), self.leading_term(ord, g))
        else {
            return self.zero();
        };
        let lcm = Monomial::lcm(&lt_f.monomial, &lt_g.monomial);
        let one = self.coeff_ring().one();
        let f_factor = self.term_poly(
            self.coeff_ring().try_divide(&one, &lt_f.coeff).unwrap(),
            lcm.try_div(&lt_f.monomial).unwrap(),
        );
        let g_factor = self.term_poly(
            self.coeff_ring().try_divide(&one, &lt_g.coeff).unwrap(),
            lcm.try_div(&lt_g.monomial).unwrap(),
        );
        self.sub(&self.mul(&f_factor, f), &self.mul(&g_factor, g))
    }

    /// Compute a Gröbner basis of the ideal generated by `generators` with
    /// respect to `ord`, using Buchberger's algorithm.
    ///
    /// The algorithm keeps a working basis and a queue of pairs. For each pair
    /// it forms the S-polynomial and reduces it modulo the current basis; a
    /// non-zero remainder exposes a leading monomial not yet captured, so it is
    /// appended to the basis and paired with every existing element. When the
    /// queue empties, every S-polynomial reduces to zero and the basis is a
    /// Gröbner basis. Buchberger's first criterion is applied to skip pairs with
    /// coprime leading monomials, whose S-polynomials are guaranteed to reduce
    /// to zero.
    ///
    /// The returned basis is correct but not canonical; use
    /// [`Self::reduced_groebner_basis`] for the unique reduced form.
    pub fn groebner_basis(
        &self,
        ord: &MonomialOrdering,
        generators: &[MultiPolynomial<FS::Elem>],
    ) -> Vec<MultiPolynomial<FS::Elem>> {
        let mut basis: Vec<MultiPolynomial<FS::Elem>> = generators
            .iter()
            .map(|g| self.reduce(g.clone()))
            .filter(|g| !self.is_zero(g))
            .collect();

        let mut pairs: Vec<(usize, usize)> = (0..basis.len())
            .flat_map(|i| ((i + 1)..basis.len()).map(move |j| (i, j)))
            .collect();

        while let Some((i, j)) = pairs.pop() {
            let lm_i = self.leading_monomial(ord, &basis[i]).unwrap();
            let lm_j = self.leading_monomial(ord, &basis[j]).unwrap();
            if Monomial::is_coprime(&lm_i, &lm_j) {
                continue;
            }
            let s = self.s_polynomial(ord, &basis[i], &basis[j]);
            let r = self.reduce_modulo(ord, &s, &basis);
            if !self.is_zero(&r) {
                let new_index = basis.len();
                for k in 0..new_index {
                    pairs.push((k, new_index));
                }
                basis.push(r);
            }
        }

        basis
    }

    /// Compute the unique *reduced* Gröbner basis of the ideal generated by
    /// `generators` with respect to `ord`.
    ///
    /// A reduced Gröbner basis is the canonical fingerprint of an ideal for a
    /// fixed ordering: two ideals are equal iff their reduced Gröbner bases
    /// coincide. It is obtained from any Gröbner basis by
    /// 1. making every element monic,
    /// 2. discarding every element whose leading monomial is divisible by the
    ///    leading monomial of another (minimalisation), and
    /// 3. replacing each remaining element by its normal form modulo the others,
    ///    so that no term of any element is divisible by the leading monomial of
    ///    another (reduction).
    ///
    /// The result is sorted by leading monomial in descending order so that the
    /// output is deterministic.
    pub fn reduced_groebner_basis(
        &self,
        ord: &MonomialOrdering,
        generators: &[MultiPolynomial<FS::Elem>],
    ) -> Vec<MultiPolynomial<FS::Elem>> {
        let basis = self.groebner_basis(ord, generators);

        // Make monic and keep only leading-monomial-minimal elements.
        let mut minimal: Vec<MultiPolynomial<FS::Elem>> = Vec::new();
        for p in basis {
            let p = self.make_monic(ord, &p);
            if self.is_zero(&p) {
                continue;
            }
            let lm_p = self.leading_monomial(ord, &p).unwrap();
            if minimal
                .iter()
                .any(|q| self.leading_monomial(ord, q).unwrap().divides(&lm_p))
            {
                continue;
            }
            minimal.retain(|q| !lm_p.divides(&self.leading_monomial(ord, q).unwrap()));
            minimal.push(p);
        }

        // Reduce each element against all the others.
        let mut reduced: Vec<MultiPolynomial<FS::Elem>> = minimal.clone();
        for i in 0..reduced.len() {
            let others: Vec<MultiPolynomial<FS::Elem>> = reduced
                .iter()
                .enumerate()
                .filter(|(k, _)| *k != i)
                .map(|(_, q)| q.clone())
                .collect();
            let r = self.reduce_modulo(ord, &reduced[i], &others);
            reduced[i] = self.make_monic(ord, &r);
        }

        reduced.sort_by(|a, b| {
            ord.cmp(
                &self.leading_monomial(ord, b).unwrap(),
                &self.leading_monomial(ord, a).unwrap(),
            )
        });
        reduced
    }

    /// Test whether `basis` is a Gröbner basis (with respect to `ord`) of the
    /// ideal it generates. By Buchberger's criterion this holds iff the
    /// S-polynomial of every pair of basis elements reduces to zero modulo the
    /// basis, which is exactly what this routine checks.
    pub fn is_groebner_basis(
        &self,
        ord: &MonomialOrdering,
        basis: &[MultiPolynomial<FS::Elem>],
    ) -> bool {
        let basis: Vec<MultiPolynomial<FS::Elem>> = basis
            .iter()
            .map(|g| self.reduce(g.clone()))
            .filter(|g| !self.is_zero(g))
            .collect();
        for i in 0..basis.len() {
            for j in (i + 1)..basis.len() {
                let s = self.s_polynomial(ord, &basis[i], &basis[j]);
                if !self.is_zero(&self.reduce_modulo(ord, &s, &basis)) {
                    return false;
                }
            }
        }
        true
    }

    /// Test whether the polynomial `f` lies in the ideal generated by
    /// `generators`. Membership is decided by computing a Gröbner basis of the
    /// ideal and checking that the normal form of `f` modulo that basis is zero
    /// (the "ideal membership problem", solved by Gröbner bases).
    pub fn ideal_contains_element(
        &self,
        ord: &MonomialOrdering,
        generators: &[MultiPolynomial<FS::Elem>],
        f: &MultiPolynomial<FS::Elem>,
    ) -> bool {
        let gb = self.groebner_basis(ord, generators);
        self.is_zero(&self.reduce_modulo(ord, f, &gb))
    }

    /// Test whether the two given generating sets generate the same ideal. Two
    /// ideals are equal iff they have the same reduced Gröbner basis with
    /// respect to any single fixed ordering.
    pub fn ideals_equal(
        &self,
        ord: &MonomialOrdering,
        generators_a: &[MultiPolynomial<FS::Elem>],
        generators_b: &[MultiPolynomial<FS::Elem>],
    ) -> bool {
        let a = self.reduced_groebner_basis(ord, generators_a);
        let b = self.reduced_groebner_basis(ord, generators_b);
        a.len() == b.len() && a.iter().zip(b.iter()).all(|(p, q)| self.equal(p, q))
    }

    /// The elimination ideal that removes the variables in `eliminate`.
    ///
    /// Given an ideal `I` of `k[x_1, ..., x_n]`, eliminating a subset of the
    /// variables yields `I ∩ k[remaining variables]`: the polynomial
    /// consequences of the generators that involve only the variables we keep.
    /// By the Elimination Theorem this is computed with a Gröbner basis for an
    /// elimination order in which any monomial containing an eliminated variable
    /// outweighs any monomial that does not (here lexicographic order with the
    /// eliminated variables most significant); the basis elements free of the
    /// eliminated variables generate the elimination ideal.
    pub fn eliminate_variables(
        &self,
        generators: &[MultiPolynomial<FS::Elem>],
        eliminate: &[Variable],
    ) -> Vec<MultiPolynomial<FS::Elem>> {
        let eliminate_set: HashSet<Variable> = eliminate.iter().cloned().collect();
        let mut other_vars: Vec<Variable> = generators
            .iter()
            .flat_map(|g| g.free_vars())
            .filter(|v| !eliminate_set.contains(v))
            .collect::<HashSet<Variable>>()
            .into_iter()
            .collect();
        other_vars.sort_by_key(Variable::ident);

        let mut order_vars = eliminate.to_vec();
        order_vars.extend(other_vars);
        let ord = MonomialOrdering::lex(order_vars);

        self.reduced_groebner_basis(&ord, generators)
            .into_iter()
            .filter(|g| g.free_vars().iter().all(|v| !eliminate_set.contains(v)))
            .collect()
    }

    /// The intersection `I ∩ J` of two ideals given by their generators.
    ///
    /// The intersection is computed with the classical elimination trick:
    /// introducing a fresh variable `t`, the ideal
    /// `t·I + (1 - t)·J = (t f_1, ..., t f_r, (1 - t) g_1, ..., (1 - t) g_s)`
    /// of `k[t, x_1, ..., x_n]` satisfies `(t·I + (1 - t)·J) ∩ k[x] = I ∩ J`.
    /// Eliminating `t` therefore returns generators of the intersection.
    pub fn ideal_intersection(
        &self,
        generators_a: &[MultiPolynomial<FS::Elem>],
        generators_b: &[MultiPolynomial<FS::Elem>],
    ) -> Vec<MultiPolynomial<FS::Elem>> {
        let t = Variable::new("t");
        let t_poly = self.var(t.clone());
        let one_minus_t = self.sub(&self.one(), &t_poly);

        let mut combined: Vec<MultiPolynomial<FS::Elem>> = Vec::new();
        for f in generators_a {
            combined.push(self.mul(&t_poly, f));
        }
        for g in generators_b {
            combined.push(self.mul(&one_minus_t, g));
        }

        self.eliminate_variables(&combined, &[t])
    }

    /// The ideal quotient (colon ideal) `(I : J) = { h : h·J ⊆ I }`.
    ///
    /// It is built from the identities `I : J = ⋂_j (I : (g_j))` over the
    /// generators `g_j` of `J`, and `I : (g) = (1/g)·(I ∩ (g))`: intersecting
    /// `I` with the principal ideal `(g)` produces generators all divisible by
    /// `g`, and dividing each by `g` gives `I : (g)`. Intersections of the
    /// resulting ideals are taken pairwise.
    pub fn ideal_quotient(
        &self,
        ord: &MonomialOrdering,
        generators_i: &[MultiPolynomial<FS::Elem>],
        generators_j: &[MultiPolynomial<FS::Elem>],
    ) -> Vec<MultiPolynomial<FS::Elem>> {
        let non_zero_j: Vec<MultiPolynomial<FS::Elem>> = generators_j
            .iter()
            .map(|g| self.reduce(g.clone()))
            .filter(|g| !self.is_zero(g))
            .collect();

        // (I : (0)) is the whole ring, and the quotient by the zero ideal is the
        // whole ring as well.
        if non_zero_j.is_empty() {
            return vec![self.one()];
        }

        let mut result: Option<Vec<MultiPolynomial<FS::Elem>>> = None;
        for g in &non_zero_j {
            let intersection = self.ideal_intersection(generators_i, std::slice::from_ref(g));
            let colon_g: Vec<MultiPolynomial<FS::Elem>> = intersection
                .iter()
                .map(|h| self.try_divide(h, g).unwrap())
                .collect();
            result = Some(match result {
                None => colon_g,
                Some(prev) => self.ideal_intersection(&prev, &colon_g),
            });
        }

        self.reduced_groebner_basis(ord, &result.unwrap())
    }

    /// Whether the ideal with the given Gröbner basis is zero-dimensional, i.e.
    /// its quotient ring is a finite-dimensional vector space over the field.
    /// This holds iff for every variable some leading monomial of the basis is a
    /// pure power of that variable (the standard finiteness criterion).
    fn is_zero_dimensional(
        &self,
        ord: &MonomialOrdering,
        gb: &[MultiPolynomial<FS::Elem>],
    ) -> bool {
        let variables: HashSet<Variable> = gb.iter().flat_map(|g| g.free_vars()).collect();
        let leading_monomials: Vec<Monomial> = gb
            .iter()
            .filter_map(|g| self.leading_monomial(ord, g))
            .collect();
        variables.iter().all(|v| {
            leading_monomials
                .iter()
                .any(|lm| lm.degree() > 0 && lm.free_vars().len() == 1 && lm.get_var_pow(v) > 0)
        })
    }

    /// Convert a Gröbner basis of a *zero-dimensional* ideal from one monomial
    /// ordering to another using the FGLM algorithm.
    ///
    /// Recomputing a Gröbner basis from scratch under a hard ordering (such as
    /// lexicographic) can be very expensive. When the ideal is zero-dimensional
    /// the quotient ring `k[x]/I` is a finite-dimensional vector space, and FGLM
    /// exploits this: it walks through monomials in increasing `target_order`,
    /// computes the normal form of each modulo the source basis as a vector in
    /// the (finite) basis of standard monomials, and uses linear algebra over
    /// the field to detect when a monomial becomes a leading term of the target
    /// basis (a linear dependence) versus a new standard monomial (independence).
    /// The cost is polynomial in the dimension of the quotient ring rather than
    /// in the doubly-exponential worst case of Buchberger's algorithm.
    ///
    /// `generators` need not already be a Gröbner basis for `source_order`; a
    /// reduced Gröbner basis is computed internally. Returns `None` when the
    /// ideal is not zero-dimensional, in which case fall back to
    /// [`Self::reduced_groebner_basis`] with the target ordering.
    pub fn convert_groebner_basis(
        &self,
        source_order: &MonomialOrdering,
        generators: &[MultiPolynomial<FS::Elem>],
        target_order: &MonomialOrdering,
    ) -> Option<Vec<MultiPolynomial<FS::Elem>>> {
        let source_gb = self.reduced_groebner_basis(source_order, generators);

        // The unit ideal (1) and the zero ideal are handled directly.
        if source_gb
            .iter()
            .any(|g| self.as_constant(g).is_some() && !self.is_zero(g))
        {
            return Some(vec![self.one()]);
        }
        if source_gb.is_empty() {
            return Some(vec![]);
        }
        if !self.is_zero_dimensional(source_order, &source_gb) {
            return None;
        }

        let variables: Vec<Variable> = {
            let mut vs: Vec<Variable> = source_gb
                .iter()
                .flat_map(|g| g.free_vars())
                .collect::<HashSet<Variable>>()
                .into_iter()
                .collect();
            vs.sort_by_key(Variable::ident);
            vs
        };

        // Normal form of a monomial modulo the source basis, as a vector
        // (map from monomial to field coefficient).
        let normal_form = |monomial: &Monomial| -> HashMap<Monomial, FS::Elem> {
            let nf = self.reduce_modulo(
                source_order,
                &self.term_poly(self.coeff_ring().one(), monomial.clone()),
                &source_gb,
            );
            nf.terms
                .into_iter()
                .map(|t| (t.monomial, t.coeff))
                .collect()
        };

        // A row of the running Gaussian elimination over the field. `residual`
        // is the part of the normal form not yet explained by earlier basis
        // monomials; `combination[k]` records how `residual` is built from the
        // normal forms of the standard monomials `staircase[k]`.
        struct EchelonRow<E> {
            residual: HashMap<Monomial, E>,
            combination: HashMap<usize, E>,
            pivot: Monomial,
        }

        let mut echelon: Vec<EchelonRow<FS::Elem>> = Vec::new();
        // Standard monomials of the target basis, in discovery order.
        let mut staircase: Vec<Monomial> = Vec::new();
        let mut target_basis: Vec<MultiPolynomial<FS::Elem>> = Vec::new();
        let mut target_leading: Vec<Monomial> = Vec::new();

        // Candidate monomials still to process, kept deduplicated.
        let mut candidates: Vec<Monomial> = vec![Monomial::one()];
        let mut seen: HashSet<Monomial> = HashSet::from([Monomial::one()]);

        // The leading monomial of a vector under the source order (used purely
        // as a pivot for elimination).
        let vector_leading = |vec: &HashMap<Monomial, FS::Elem>| -> Option<Monomial> {
            vec.keys().max_by(|&a, &b| source_order.cmp(a, b)).cloned()
        };

        while !candidates.is_empty() {
            // Pick the smallest candidate in the target order.
            let idx = (0..candidates.len())
                .min_by(|&i, &j| target_order.cmp(&candidates[i], &candidates[j]))
                .unwrap();
            let monomial = candidates.swap_remove(idx);

            // Skip monomials already divisible by a target leading monomial:
            // they are not standard and cannot be new leading monomials.
            if target_leading.iter().any(|lm| lm.divides(&monomial)) {
                continue;
            }

            let mut residual = normal_form(&monomial);
            let mut combination: HashMap<usize, FS::Elem> = HashMap::new();

            // Reduce the residual against the current echelon rows.
            while let Some(lead) = vector_leading(&residual) {
                let Some(row) = echelon.iter().find(|r| r.pivot == lead) else {
                    break;
                };
                let factor = self
                    .coeff_ring()
                    .try_divide(&residual[&lead], &row.residual[&lead])
                    .unwrap();
                for (mon, coeff) in &row.residual {
                    let entry = residual
                        .entry(mon.clone())
                        .or_insert_with(|| self.coeff_ring().zero());
                    *entry = self
                        .coeff_ring()
                        .sub(entry, &self.coeff_ring().mul(&factor, coeff));
                }
                for (k, coeff) in &row.combination {
                    let entry = combination
                        .entry(*k)
                        .or_insert_with(|| self.coeff_ring().zero());
                    *entry = self
                        .coeff_ring()
                        .add(entry, &self.coeff_ring().mul(&factor, coeff));
                }
                residual.retain(|_, c| !self.coeff_ring().is_zero(c));
            }
            residual.retain(|_, c| !self.coeff_ring().is_zero(c));

            if residual.is_empty() {
                // Linear dependence: `monomial` equals a combination of standard
                // monomials, so `monomial - Σ c_k · staircase[k]` is in the
                // ideal and has leading monomial `monomial` in the target order.
                let mut g = self.term_poly(self.coeff_ring().one(), monomial.clone());
                for (k, coeff) in &combination {
                    let term = self.term_poly(coeff.clone(), staircase[*k].clone());
                    g = self.sub(&g, &term);
                }
                target_leading.push(monomial);
                target_basis.push(self.make_monic(target_order, &g));
            } else {
                // Linear independence: a new standard monomial.
                let new_index = staircase.len();
                echelon.push(EchelonRow {
                    pivot: vector_leading(&residual).unwrap(),
                    residual,
                    combination: {
                        let mut c = combination
                            .into_iter()
                            .map(|(k, v)| (k, self.coeff_ring().neg(&v)))
                            .collect::<HashMap<usize, FS::Elem>>();
                        c.insert(new_index, self.coeff_ring().one());
                        c
                    },
                });
                staircase.push(monomial.clone());
                for v in &variables {
                    let next = Monomial::mul(
                        &monomial,
                        &Monomial::new(vec![VariablePower {
                            var: v.clone(),
                            pow: 1,
                        }]),
                    );
                    if seen.insert(next.clone()) {
                        candidates.push(next);
                    }
                }
            }
        }

        target_basis.sort_by(|a, b| {
            target_order.cmp(
                &self.leading_monomial(target_order, b).unwrap(),
                &self.leading_monomial(target_order, a).unwrap(),
            )
        });
        Some(target_basis)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomic;

    /// Build a monomial from `(variable, exponent)` pairs.
    fn mon(pairs: &[(&Variable, usize)]) -> Monomial {
        Monomial::new(
            pairs
                .iter()
                .map(|(v, p)| VariablePower {
                    var: (*v).clone(),
                    pow: *p,
                })
                .collect(),
        )
    }

    #[test]
    fn test_leading_term_per_ordering() {
        // The textbook example f = 4xy^2z + 4z^2 - 5x^3 + 7x^2z^2 has a
        // different leading monomial under each of the three orderings.
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let zv = Variable::new("z");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let z = MultiPolynomial::<Rational>::var(zv.clone()).into_ergonomic();

        let f = (4 * &x * y.pow(2) * &z + 4 * z.pow(2) - 5 * x.pow(3) + 7 * x.pow(2) * z.pow(2))
            .into_verbose();
        let mp = MultiPolynomial::<Rational>::structure();

        let vars = vec![xv.clone(), yv.clone(), zv.clone()];
        let lex = MonomialOrdering::lex(vars.clone());
        let grlex = MonomialOrdering::grlex(vars.clone());
        let grevlex = MonomialOrdering::grevlex(vars.clone());

        assert_eq!(
            mp.leading_monomial(&lex, &f).unwrap(),
            mon(&[(&xv, 3)]),
            "lex leading monomial should be x^3"
        );
        assert_eq!(
            mp.leading_monomial(&grlex, &f).unwrap(),
            mon(&[(&xv, 2), (&zv, 2)]),
            "grlex leading monomial should be x^2 z^2"
        );
        assert_eq!(
            mp.leading_monomial(&grevlex, &f).unwrap(),
            mon(&[(&xv, 1), (&yv, 2), (&zv, 1)]),
            "grevlex leading monomial should be x y^2 z"
        );
    }

    #[test]
    fn test_monomial_helpers() {
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x2y = mon(&[(&xv, 2), (&yv, 1)]);
        let xy = mon(&[(&xv, 1), (&yv, 1)]);
        let xz = mon(&[(&xv, 1), (&Variable::new("z"), 1)]);

        assert!(xy.divides(&x2y));
        assert!(!x2y.divides(&xy));
        assert_eq!(x2y.try_div(&xy).unwrap(), mon(&[(&xv, 1)]));
        assert!(x2y.try_div(&xz).is_none());
        assert_eq!(Monomial::lcm(&x2y, &xy), x2y);
        assert!(!Monomial::is_coprime(&x2y, &xy));
        assert!(Monomial::is_coprime(&mon(&[(&xv, 1)]), &mon(&[(&yv, 1)])));
    }

    #[test]
    fn test_buchberger_and_verification() {
        // Cox, Little, O'Shea: with grlex (x > y) the set {f1, f2} below is NOT
        // a Gröbner basis; its reduced Gröbner basis is {x^2, xy, y^2 - x/2}.
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();

        let f1 = (x.pow(3) - 2 * &x * &y).into_verbose();
        let f2 = (x.pow(2) * &y - 2 * y.pow(2) + &x).into_verbose();
        let mp = MultiPolynomial::<Rational>::structure();
        let grlex = MonomialOrdering::grlex(vec![xv.clone(), yv.clone()]);

        let gens = [f1.clone(), f2.clone()];
        assert!(
            !mp.is_groebner_basis(&grlex, &gens),
            "the raw generators are not a Gröbner basis"
        );

        let gb = mp.reduced_groebner_basis(&grlex, &gens);
        assert!(
            mp.is_groebner_basis(&grlex, &gb),
            "the computed reduced basis must be a Gröbner basis"
        );
        // The generators lie in the ideal, and the basis generates the same ideal.
        assert!(mp.ideal_contains_element(&grlex, &gb, &f1));
        assert!(mp.ideal_contains_element(&grlex, &gb, &f2));
        assert!(mp.ideals_equal(&grlex, &gens, &gb));

        // The reduced basis has the expected leading monomials x^2, xy, y^2.
        let mut leading: Vec<Monomial> = gb
            .iter()
            .map(|g| mp.leading_monomial(&grlex, g).unwrap())
            .collect();
        leading.sort_by(|a, b| grlex.cmp(b, a));
        assert_eq!(
            leading,
            vec![
                mon(&[(&xv, 2)]),
                mon(&[(&xv, 1), (&yv, 1)]),
                mon(&[(&yv, 2)]),
            ]
        );
    }

    #[test]
    fn test_membership_uses_normal_form() {
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let mp = MultiPolynomial::<Rational>::structure();
        let grevlex = MonomialOrdering::grevlex(vec![xv.clone(), yv.clone()]);

        let gens = [(x.pow(2) - &y).into_verbose(), (&y - 1).into_verbose()];
        // x^2 - 1 = (x^2 - y) + (y - 1) lies in the ideal.
        assert!(mp.ideal_contains_element(&grevlex, &gens, &(x.pow(2) - 1).into_verbose()));
        // x is not in the ideal.
        assert!(!mp.ideal_contains_element(&grevlex, &gens, &x.clone().into_verbose()));
    }

    #[test]
    fn test_fglm_conversion_matches_direct() {
        // A zero-dimensional ideal: x^2 + y^2 = 4 and xy = 1.
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let mp = MultiPolynomial::<Rational>::structure();

        let gens = [
            (x.pow(2) + y.pow(2) - 4).into_verbose(),
            (&x * &y - 1).into_verbose(),
        ];
        let grevlex = MonomialOrdering::grevlex(vec![xv.clone(), yv.clone()]);
        let lex = MonomialOrdering::lex(vec![xv.clone(), yv.clone()]);

        let converted = mp
            .convert_groebner_basis(&grevlex, &gens, &lex)
            .expect("ideal is zero-dimensional");
        let direct = mp.reduced_groebner_basis(&lex, &gens);

        assert!(mp.is_groebner_basis(&lex, &converted));
        assert_eq!(converted.len(), direct.len());
        assert!(
            converted
                .iter()
                .zip(direct.iter())
                .all(|(a, b)| mp.equal(a, b)),
            "FGLM conversion must match the directly computed lex Gröbner basis"
        );
    }

    #[test]
    fn test_fglm_rejects_positive_dimensional() {
        // The ideal (xy) is one-dimensional, so FGLM is not applicable.
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let mp = MultiPolynomial::<Rational>::structure();

        let gens = [(&x * &y).into_verbose()];
        let grevlex = MonomialOrdering::grevlex(vec![xv.clone(), yv.clone()]);
        let lex = MonomialOrdering::lex(vec![xv.clone(), yv.clone()]);
        assert!(mp.convert_groebner_basis(&grevlex, &gens, &lex).is_none());
    }

    #[test]
    fn test_ideal_intersection() {
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let mp = MultiPolynomial::<Rational>::structure();
        let grevlex = MonomialOrdering::grevlex(vec![xv.clone(), yv.clone()]);

        // (x) ∩ (y) = (xy).
        let inter = mp.ideal_intersection(&[x.clone().into_verbose()], &[y.clone().into_verbose()]);
        assert!(mp.ideals_equal(&grevlex, &inter, &[(&x * &y).into_verbose()]));

        // (x) ∩ (x + y) = (x^2 + xy) since x and x + y are coprime.
        let inter2 =
            mp.ideal_intersection(&[x.clone().into_verbose()], &[(&x + &y).into_verbose()]);
        assert!(mp.ideals_equal(&grevlex, &inter2, &[(x.pow(2) + &x * &y).into_verbose()]));
    }

    #[test]
    fn test_ideal_quotient() {
        let xv = Variable::new("x");
        let yv = Variable::new("y");
        let x = MultiPolynomial::<Rational>::var(xv.clone()).into_ergonomic();
        let y = MultiPolynomial::<Rational>::var(yv.clone()).into_ergonomic();
        let mp = MultiPolynomial::<Rational>::structure();
        let grevlex = MonomialOrdering::grevlex(vec![xv.clone(), yv.clone()]);

        // (x^2 y) : (x) = (xy).
        let quotient = mp.ideal_quotient(
            &grevlex,
            &[(x.pow(2) * &y).into_verbose()],
            &[x.clone().into_verbose()],
        );
        assert!(mp.ideals_equal(&grevlex, &quotient, &[(&x * &y).into_verbose()]));
    }
}
