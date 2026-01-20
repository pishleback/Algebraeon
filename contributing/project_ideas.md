# Project Ideas

This document contains a list of project ideas for those wishing to take on a project to contribute a rewarding chunk of work to Algebraeon, with mentoring.

### General Requirements
To take on a projects with us there are some general requirements:
 - You should be comfortable writing complex programs in a high level programming language.
 - You need to be comfortable with the mathematical concepts underpinning a particular project. Generally this will mean some combination of: Ring Theory, Linear Algebra, Modular Arithmetic, Number Theory, Group Theory, Galois Theory, and Algebraic Number Theory.

Bonus points if you are already familiar with programming in Rust.

## Finding The Four Squares of Lagrange's Theorem
A theorem of Lagrange states that for every non-negative integer $n \in \mathbb{Z}$ there exists $4$ square numbers $a^2$, $b^2$, $c^2$ and $d^2$ such that
```math
n = a^2 + b^2 + c^2 + d^2
```
In other words, every non-negative integer can be expressed a sum of four square numbers. One enlightening proof of this fact makes use of integer quaternions and some modular arithmetic. For this project you will use this proof to implement an algorithm which efficiently expresses any non-negative integer a sum of four square numbers, including a new framework for working with quaternion integers in Algebraeon.

### Expected Outcomes
 - A new framework for working with quaternions and integer quaternions.
 - Algorithms to efficiently solve the necessary problems in modular arithmetic.
 - Using the above, an algorithm to express any non-negative integer, even extremely large ones, as a sum of four squares.

### Difficulty
Easy

### Duration
~175 hours

### Requirements
You must be familiar with:
 - Modular arithmetic.
 - Ring theory.

It is preferable if you are familiar with:
 - The algebraic structure of quaternions.
 - Euclidean domains and Euclidean division in $\mathbb{Z}[i]$, as it provides a nice analogue to the integer quaternions.

### Possible Mentors
 - Pishleback
 

## Generalized Pell Equations
For a fixed non-square positive integer $d \in \mathbb{Z}$ and non-zero integer $n \in \mathbb{Z}$, the generalized Pell equation is
```math
x^2 - dy^2 = n
```
where $x$ and $y$ are integers to be found.
 - An (ungeneralized) Pell equation is the special case $n=1$. The Pell equation $x^2 - dy^2 = 1$ always has infinitely many solutions, and there is a so-called fundamental solution from which all other solutions are generated in a certain sense.
 - For the generalized Pell equation $x^2 - dy^2 = n$ there is always a finite set of solutions from which all others can be generated, once a fundamental solution to the corresponding Pell $x^2 - dy^2 = 1$ equation is known.

### Expected Outcomes
Main outcome:
 - An algorithm to find the fundamental solution to Pell's equation $x^2 - dy^2 = 1$ using the existing capability in Algebraeon for computing continued fractions.

Stretch outcome:
 - An algorithm to find the finite set of solutions to a generalized Pell equation $x^2 - dy^2 = n$ from which all others are generated.

### Difficulty
Easy

### Duration
~90 hours

### Requirements
You must be familiar with:
 - Number theory.

It is preferable if you are familiar with:
 - Continued fractions including the eventual periodicity of continued fractions on quadratic values.

### Possible Mentors
 - Pishleback
 

## Computing Class Groups in Quadratic Number Fields
A quadratic number field is the ring of integers of an algebraic number field of the form $\mathbb{Q}[\sqrt{d}]$ where $\sqrt{d}$ is a square-free integer other than $1$. The ring of integers of the quadratic number field is given by
```math
\mathcal{O}_d = \begin{cases} \mathbb{Z}[\frac{1 + \sqrt{d}}{2}] & d \equiv 1 \pmod 4 \\ \mathbb{Z}[\sqrt{d}] & d \equiv 2, 3 \pmod 4  \end{cases}
```
The fractional ideals of $\mathcal{O}_d$ form an abelian group and the subgroup of principal fractional ideals has finite index. The finite quotient group of the fractional ideals by the principal fractional ideals is called the class group.

In this project you will use the theory of binary integral quadratic forms as a means to compute the size and structure of the class group. The approach splits naturally into two cases; the imaginary case when $d < 0$ and the real case when $d > 0$.

### Expected Outcomes
Main outcome:
 - An algorithm to compute the size and structure of the class group of imaginary quadratic number fields using reduction theory for binary integral quadratic forms.

Likely secondary outcomes:
 - An algorithm to compute the size and structure of the class group of real quadratic number fields using reduction theory for binary integral quadratic forms.
 - Solving the principal ideal problem for quadratic number fields - an algorithm which determines whether an ideal in $\mathcal{O}_d$ is principal, and return a generator if it is.

Stretch outcomes:
 - Look further into the theory of binary integral quadratic forms and exposing some more of their structure, for example Genus theory.
 - Faster algorithms exist for computing the class group of quadratic number fields, and could also be implemented.

### Difficulty
Medium

### Duration
~175 hours

### Requirements
You must be familiar with:
 - Group theory.
 - Quadratic number fields including ideals, fractional ideals, and the class group.

It is preferable if you are familiar with:
 - Binary integral quadratic forms.


### Possible mentors:
 - Pishleback
 

<!-- ## Factoring Integers using the Quadratic Number Field Seive
basic proof-of-concept algorithm

get it faster than ECM using optimizations like MPQNFS


## Computations in the Algebraic Closure of Finite Fields
open ended explore ways to compute in the algebraic closure of finite fields -->


## Faster Factoring of Integer Polynomials
The Berlekamp–Zassenhaus algorithm is currently used by Algebraeon for factoring integer polynomials. The essential steps used to factor a polynomial $f(x) \in \mathbb{Z}[x]$ are as follows:
 - $f(x)$ is factored into irreducibles modulo some prime $p$, say $f(x) = f_1(x) \cdots f_n(x) \pmod p$ with $f_1(x), \dots, f_n(x)$ the irreducible factors of $f(x)$ modulo $p$. This can already be done very quickly in Algebraeon using the Cantor–Zassenhaus algorithm, whose inner workings need not be understood for this project. The polynomials $f_1(x), \dots, f_n(x)$ are called the modular factors of $f(x)$.
 - The factorization $f(x) = f_1(x) \cdots f_n(x) \pmod p$ mod $p$ is "lifted" via Hensel's lemma to factorizations modulo $p^k$ for increasingly large $k$. In a $p$-adic sense, these lifted factorizations give increasingly good approximations to the desired factorization of $f(x)$ in $\mathbb{Z}[x]$.
 - There is not a 1-to-1 correspondence between the modular factors of $f(x)$ and its true factors. Rather, each true factor corresponds to some subset of the modular factors, and the set of all true factors of $f(x)$ corresponds to a partition of the modular factors. To finish the factorization of $f(x)$, the Berlekamp–Zassenhaus algorithm loops over all subsets of the modular factors and checks whether they correspond to a true irreducible factor of $f(x)$.

This algorithm performs quite well in practice. However, in the last step, the number of subsets of modular factors to check grows like $2^n$. 

Mark van Hoeij's Knapsack Algorithm is a state-of-the-art improvement to the Berlekamp–Zassenhaus algorithm turns the exponential-time task of checking all subsets of modular factors into a knapsack type problem which can be solved in polynomial type using LLL. A good implementation of the Knapsack algorithm can vastly speed up the factorization of integer polynomials compared to Berlekamp–Zassenhaus. A description of the algorithm can be found here https://www.math.fsu.edu/~hoeij/knapsack/paper/knapsack.pdf.

### Expected Outcomes
Main outcome:
 - You will implement the Knapsack algorithm for factoring integer polynomials in Algebraeon.

Stretch outcome:
 - Paul Zimmermann has a list of benchmark polynomials here https://homepages.loria.fr/PZimmermann/mupad/. It would be delightful if Algebraeon can factor all 8 of these polynomials in a reasonable time. The current implementation of the Berlekamp–Zassenhaus algorithm can only factor the first two.

### Difficulty
Hard

### Duration
~350 hours

### Requirements
You must be familiar with:
 - Linear algebra.
 - Ring theory including UFDs and polynomial rings.
 - Finite fields.

It is preferable if you are familiar with:
 - $p$-adic numbers.
 - Hensel lifting.
 - LLL lattice basis reduction.
 - Details of the Berlekamp–Zassenhaus algorithm for factoring integer polynomials.

### Possible Mentors
 - Pishleback

## Your Own Project Proposal
You can come to us with your own project proposal. Perhaps there is an interesting algorithm you've always wanted to implement? Or something you find useful in another CAS which Algebraeon is not yet capable of? Get in touch.