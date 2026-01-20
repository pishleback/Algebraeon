# Project Ideas

This document contains a list of project ideas for anyone wishing to take on a project to contribute a rewarding chunk of work to Algebraeon with mentoring.


## Finding Four Squares of Lagrange's Theorem
For every non-negative integer $n$ there exists $4$ square numbers $a^2$, $b^2$, $c^2$ and $d^2$ such that
```math
n = a^2 + b^2 + c^2 + d^2
```
In other words, every non-negative integer can be expressed a sum of four square numbers. One of the more enlightening proofs of this fact makes use of integer quaternions and some modular arithmetic. For this project you will implement this proof, including a new framework for working with quaternions, in Algebraeon.

Expected outcomes:
 - A framework for working with quaternions and integer quaternions.
 - Algorithms to efficiently solve the necessary problems in modular arithmetic.
 - Using the above, an algorithm to express any non-negative integer, even extremely large ones, as a sum of four squares.

Relative difficulty: Easy

Expected duration: ~90 hours

You must be familiar with:
 - Modular arithmetic.
 - The basics of ring theory.
 - Writing complex programs in a high-level language, such as Python.

It is preferable if you are familiar with:
 - The algebraic structure of quaternions.
 - Euclidean domains and Euclidean division in $\mathbb{Z}[i]$, as it provides a nice analogue to the integer quaternions.
 - Programming in Rust.

Possible mentors:
 - Pishleback
 

## Generalized Pell Equations
For a fixed non-square positive integer $d$ and non-zero integer $n$, the generalized Pell equation is
```math
x^2 - dy^2 = n
```
where $x$ and $y$ are integers to be found. Pell's (ungeneralized) equation is the special case $n=1$. 
 - An (ungeneralized) Pell equation is the special case $n=1$. A Pell equation always has infinitely many solutions, and there is a so-called fundamental solution from which all other solutions are generated.
 - For arbitrary $n$ there is always a finite set of solutions from which all others can be generated once a fundamental solution to the corresponding Pell's equation is known.

Expected outcomes:
 - An algorithm to find the fundamental solution to Pell's equation using the existing capability in Algebraeon for computing continued fractions.
 - An algorithm to find the finite set of solutions to a generalized Pell equation, from which all others are generated.

Relative difficulty: Easy

Expected duration: ~90 hours

You must be familiar with:
 - Continued fractions.
 - Writing complex programs in a high-level language, such as Python.

It is preferable if you are familiar with:
 - The eventual periodicity of continued fractions on quadratic values.
 - Programming in Rust.

Possible mentors:
 - Pishleback
 

## Class Groups in Quadratic Number Fields (TODO)

Reduction theory for binary integral quadratic forms.

TODO: Compute the class group of an quadratic number field.

Relative difficulty: Medium

Expected duration: ~175 hours

You must be familiar with:
 - Group theory.
 - Quadratic number fields including the unique factorization of ideals and finiteness of the class group.
 - Writing complex programs in a high-level language, such as Python.

It is preferable if you are familiar with:
 - Quadratic forms.
 - Programming in Rust.

Possible mentors:
 - Pishleback
 


## Factoring Integers using the Quadratic Number Field Seive
basic proof-of-concept algorithm

get it faster than ECM using optimizations like MPQNFS


## Computations in the Algebraic Closure of Finite Fields
open ended explore ways to compute in the algebraic closure of finite fields


## Faster Factoring Integer Polynomials
The Berlekamp–Zassenhaus algorithm is currently used by Algebraeon for factoring integer polynomials. The essential steps used to factor a polynomial $f(x) \in \mathbb{Z}$ are as follows:
 - $f(x)$ is factored into irreducibles modulo some prime $p$, say $f(x) = f_1(x) \cdots f_n(x) \pmod p$ with $f_1(x), \dots, f_n(x)$ the irreducible factors of $f(x)$ modulo $p$. This can already be done very quickly in Algebraeon using the Cantor–Zassenhaus algorithm, whose inner workings need not be understood for this project. The polynomials $f_1(x), \dots, f_n(x)$ are called the modular factors of $f(x)$.
 - The factorization $f(x) = f_1(x) \cdots f_n(x) \pmod p$ mod $p$ is "lifted" via Hensel's lemma to factorizations modulo $p^k$ for increasingly large $k$. In a $p$-adic sense, these lifted factorizations give increasingly good approximations to the desired factorization of $f(x)$ in $\mathbb{Z}[x]$.
 - There is not a 1-to-1 correspondence between the modular factors of $f(x)$ and its true factors. Rather, each true factor corresponds to some subset of the modular factors, and the set of all true factors of $f(x)$ corresponds to a partition of the modular factors. To finish the factorization of $f(x)$, the Berlekamp–Zassenhaus algorithm loops over all subsets of the modular factors and checks whether they correspond to a true irreducible factor of $f(x)$.

This algorithm performs quite well in practice. However, in the last step, the number of subsets of modular factors to check grows like $2^n$. 

Mark van Hoeij's Knapsack Algorithm is a state-of-the-art improvement to the Berlekamp–Zassenhaus algorithm turns the exponential-time task of checking all subsets of modular factors into a knapsack type problem which can be solved in polynomial type using LLL. A good implementation of the Knapsack algorithm can vastly speed up the factorization of integer polynomials compared to Berlekamp–Zassenhaus. The algorithm is described here https://www.math.fsu.edu/~hoeij/knapsack/paper/knapsack.pdf.

Outcomes:
 - You will implement the Knapsack algorithm for factoring integer polynomials in Algebraeon.

Stretch outcomes:
 - Paul Zimmermann has a list of benchmark polynomials here https://homepages.loria.fr/PZimmermann/mupad/. It would be delightful if Algebraeon can factor all 8 of these polynomials in a reasonable time. The current implementation of the Berlekamp–Zassenhaus algorithm can only factor the first two.

Relative difficulty: Hard

Expected duration: ~350 hours

You must be familiar with:
 - Linear algebra.
 - Ring theory including UFDs and polynomial rings.
 - Finite fields.
 - $p$-adic numbers.

It is preferable if you are familiar with:
 - The Berlekamp–Zassenhaus algorithm for factoring integer polynomials.
 - Hensel lifting.
 - LLL lattice basis reduction.
 - Programming in Rust.

Possible mentors:
 - Pishleback


## <TITLE/DESCRIPTION>
A more detailed description of the project (2-5 sentences).

Expected outcomes.

Estimated time and difficulty (easy / medium / hard). ~90 hours, ~175 hours or ~350 hours.

Required skills:
 - x
 - y

Preferred skills:
 - a
 - b

Possible mentors:
 - Pishleback