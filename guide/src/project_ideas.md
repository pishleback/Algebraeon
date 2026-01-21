# Project Ideas

This document contains a list of project ideas for those wishing to take on a project to contribute a chunk of work to Algebraeon, with mentoring.

**General Requirements**
To take on a project with mentoring with us there are some general requirements:
 - You should be comfortable writing complex programs in a high level programming language.
 - You need to be comfortable with the mathematical concepts underpinning a particular project. Generally this will mean some combination of: Ring Theory, Linear Algebra, Modular Arithmetic, Number Theory, Group Theory, Galois Theory, and Algebraic Number Theory.

Bonus points if you are already familiar with programming in Rust.

## Finding The Four Squares of Lagrange's Theorem
A theorem of Lagrange states that for every non-negative integer \\(n \in \mathbb{Z}\\) there exist \\(4\\) square numbers \\(a^2\\), \\(b^2\\), \\(c^2\\) and \\(d^2\\) such that
\\[n = a^2 + b^2 + c^2 + d^2\\]
In other words, every non-negative integer can be expressed a sum of four square numbers. One enlightening proof of this fact makes use of integer quaternions and some modular arithmetic. For this project you will use this proof to implement an algorithm which efficiently expresses any non-negative integer a sum of four square numbers, including a framework for working with quaternion integers in Algebraeon.

**Expected Outcomes**
 - A framework for working with quaternions and integer quaternions.
 - Algorithms to efficiently solve the necessary problems in modular arithmetic.
 - Using the above, an algorithm to express non-negative integers, even extremely large ones, as a sum of four squares.

**Difficulty:** Easy

**Duration:** ~175 hours

**Requirements**

You must be familiar with:
 - Modular arithmetic.
 - Ring theory.

It is preferable if you are familiar with:
 - The algebraic structure of quaternions.
 - Euclidean domains and Euclidean division in $\mathbb{Z}[i]$, as it provides a nice analogue to the integer quaternions.

**Possible Mentors**
 - Pishleback
 

## Generalized Pell Equations
For a fixed non-square positive integer \\(d \in \mathbb{Z}\\) and non-zero integer \\(n \in \mathbb{Z}\\), the generalized Pell equation is the diophantine equation
\\[x^2 - dy^2 = n\\]
where \\(x\\) and \\(y\\) are integers to be found.
 - An (ungeneralized) Pell equation is the special case \\(n=1\\). The Pell equation \\(x^2 - dy^2 = 1\\) always has infinitely many solutions, and there is a so-called fundamental solution from which all other solutions are generated in a certain sense.
 - For the generalized Pell equation \\(x^2 - dy^2 = n\\) there is always a finite set of solutions from which all others can be generated, once a fundamental solution to the corresponding Pell equation \\(x^2 - dy^2 = 1\\) is known.

**Expected Outcomes**

Main outcomes:
 - Find the fundamental solution to Pell's equation \\(x^2 - dy^2 = 1\\) using the existing capability in Algebraeon for computing continued fractions.
 - Generate all other solutions to Pell's equation once the fundamental solution is found.

Stretch outcomes:
 - Find the finite set of solutions to a generalized Pell equation \\(x^2 - dy^2 = n\\) from which all others are generated.
 - Generate all solutions to a generalized Pell's equation.

**Difficulty:** Easy

**Duration:** ~90 hours

**Requirements**

You must be familiar with:
 - Number theory.

It is preferable if you are familiar with:
 - Continued fractions including the eventual periodicity of continued fractions on quadratic values.

**Possible Mentors**
 - Pishleback
 

## More Efficient Polynomials over Finite Fields
The simplest way to represent a polynomial \\(c_0 + c_1 x + c_2 x^2 + \dots + c_d x^d \in R[x]\\) over a ring \\(R\\) is by storing a `Vec` of coefficients \\(c_0, c_1, c_2, \dots, c_d \in R\\). This is how Algebraeon currently does it, but for certain rings such as \\(R=\frac{\mathbb{Z}}{n \mathbb{Z}}\\) it is possible to do better. For example, over the ring \\(R = \frac{\mathbb{Z}}{2 \mathbb{Z}}\\) polynomials can be more efficiently represented as the binary expansion of a natural number - polynomial addition becomes a bitwise XOR and polynomial multiplication becomes multiplication of natural numbers. Optimizations of a similar kind can be made for \\(R=\frac{\mathbb{Z}}{n \mathbb{Z}}\\) when \\(n > 2\\) too.

**Expected Outcomes**

Main outcome:
 - Efficient representations of polynomials over \\(\frac{\mathbb{Z}}{2 \mathbb{Z}}\\) with coefficients represented internally as the bits in the binary expantion of a natural number.

Secondary outcomes:
 - Efficient representations of polynomials over \\(\frac{\mathbb{Z}}{n \mathbb{Z}}\\) for \\(n > 2\\) with coefficients represented internally as blocks of bits in the binary expantion of a natural number.
 - Refine Algebraeon's polynomials structure system to allow general-purpose polynomial algorithms to operate on polynomials represented in this way.

Stretch outcomes:
 - Make use of the new polynomial representations in the existing implementation of the Berlekamp–Zassenhaus algorithm for factoring polynomials over the integers. Do we see the expected performance boost in the benchmarks?

**Difficulty:** Medium

**Duration:** ~175 hours

**Requirements**

You must be familiar with:
 - Modular arithmetic.
 - Binary expansions.

It is preferable if you are familiar with:
 - Writing optimized low-level code exploiting the binary representations of data structures. 

**Possible Mentors**

 - Pishleback
 

## Computing The Class Group of Quadratic Number Fields
A quadratic number field is the ring of integers of an algebraic number field of the form \\(\mathbb{Q}[\sqrt{d}]\\) where \\(\sqrt{d}\\) is a square-free integer (excluding \\(1\\)). The ring of integers of the quadratic number field is given by
\\[\mathcal{O}_d = \begin{cases} \mathbb{Z}[\frac{1 + \sqrt{d}}{2}] & d \equiv 1 \pmod 4 \\\\ \mathbb{Z}[\sqrt{d}] & d \equiv 2, 3 \pmod 4  \end{cases}\\]
The class group of \\(\mathcal{O}_d\\) is the quotient of the fractional ideals by the principal fractional ideals. It is a theorem of algebraic number theory that the class group is always finite.

For this project you will use the reduction theory of binary integral quadratic forms as a means to compute the size and structure of the class group.

**Expected Outcomes**

Main outcome:
 - An algorithm to compute the size and structure of the class group of imaginary quadratic number fields using reduction theory for binary integral quadratic forms.

Secondary outcomes:
 - An algorithm to compute the size and structure of the class group of real quadratic number fields using reduction theory for binary integral quadratic forms.
 - Solving the principal ideal problem for quadratic number fields - an algorithm which determines whether an ideal in \\(\mathcal{O}_d\\) is principal, and return a generator if it is.

Stretch outcomes:
 - Look further into the theory of binary integral quadratic forms and exposing some more of their structure, for example Genus theory.
 - Faster algorithms exist for computing the class group of quadratic number fields, and could also be implemented.

**Difficulty:** Medium

**Duration:** ~175 hours

**Requirements**

You must be familiar with:
 - Group theory.
 - Quadratic number fields including ideals, fractional ideals, and the class group.

It is preferable if you are familiar with:
 - Binary integral quadratic forms.

**Possible mentors**
 - Pishleback



## Working in The Algebraic Closure of Finite Fields
Let \\(p \in \mathbb{N}\\) be a prime number. For any positive integer \\(k\\) there is a unique finite field \\(\mathbb{F}\_{p^k}\\) of order \\(p^k\\), and \\(\mathbb{F}\_{p^a}\\) is a subfield of \\(\mathbb{F}\_{p^b}\\) if and only if \\(a\\) divides \\(b\\). The algebraic closure of \\(\mathbb{F}\_p\\) can be viewed as the inverse limit with respect to the subfield inclusions \\(\mathbb{F}\_{p^a} \subseteq \mathbb{F}\_{p^b}\\) of all the finite fields \\(\mathbb{F}\_{p^k}\\)
\\[\overline{\mathbb{F}}_p = \lim\_{\underset{k \ge 1}{\longleftarrow}} \mathbb{F}\_{p^k}\\]
There are a variety of techniques for doing computations in \\(\overline{\mathbb{F}}_p\\) with their own advantages. The use of Conway polynomials is one way, which Algebraeon already implements.

In this project you will explore and implement other methods for computing in \\(\overline{\mathbb{F}}_p\\).

**Difficulty:** Medium/Hard

**Duration:** ~175 hours

**Requirements**

You must be familiar with:
 - Group theory.
 - Field theory including the Galois theory of finite fields.
 - The construction of an algebraic closure as an inverse limit.

**Possible mentors**

 - Pishleback
 


## Faster Factoring of Integer Polynomials
The Berlekamp–Zassenhaus algorithm is currently used by Algebraeon for factoring integer polynomials. At a high level it works by factoring the polynomial modulo some prime \\(p\\), lifting that factorization so that it holds modulo higher powers of \\(p\\), and finally finding which subsets of the irreducible factors combine to produce irreducible factors over the integers. The algorithm is very fast when the number of factors modulo \\(p\\) is small. However, due to the search over all subsets of factors modulo \\(p\\), the running time is exponential in the number of factors modulo \\(p\\). Mark van Hoeij's Knapsack Algorithm is a state-of-the-art improvement to the Berlekamp–Zassenhaus algorithm which turns the exponential-time task of checking all subsets of modular factors into a knapsack type problem which can be solved in polynomial type using LLL.

**Expected Outcomes**

Main outcome:
 - You will implement version of the Knapsack algorithm for factoring integer polynomials in Algebraeon.

Stretch outcome:
 - Paul Zimmermann has a list of benchmark polynomials [here](https://homepages.loria.fr/PZimmermann/mupad/). It would be delightful if Algebraeon could factor all 8 of these polynomials in a reasonable time. The current implementation of the Berlekamp–Zassenhaus algorithm can factor the first two.

**Difficulty:** Hard

**Duration:** ~350 hours

**Requirements**

You must be familiar with:
 - Linear algebra.
 - Ring theory including UFDs and polynomial rings.
 - Finite fields.

It is preferable if you are familiar with:
 - \\(p\\)-adic numbers.
 - Hensel lifting.
 - LLL lattice basis reduction.
 - Details of the Berlekamp–Zassenhaus algorithm for factoring integer polynomials.

**Possible Mentors**
 - Pishleback

## Your Own Project Proposal
You can come to us with your own project proposal. Perhaps there is an interesting algorithm you've always wanted to implement? Or something you find useful in another CAS which Algebraeon is not yet capable of? Get in touch.