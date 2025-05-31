# Free Modules

This section explains how to work with free modules, submodules and cosets of submodules.

If \\(R\\) is a commutative ring then the set \\(R^n\\) of all length \\(n\\) tuples of elements from \\(R\\) is a free \\(R\\)-module with standard basis \\(e_1, \dots, e_n \in R^n\\).
\\[e_1 = (1, 0, \dots, 0)\\]
\\[e_2 = (0, 1, \dots, 0)\\]
\\[\vdots\\]
\\[e_n = (0, 0, \dots, 1)\\]

Algebraeon supports working with free modules over the following rings \\(R\\):
 - The integers \\(\mathbb{Z}\\)
 - The rationals \\(\mathbb{Q}\\)
 - Any field
 - Any Euclidean domain

The examples in this section primarily illustrate how to use Algebraeon in the case \\(R = \mathbb{Z}\\). 

## Free Modules

The free module \\(R^n\\) is represented by objects of type `FinitelyFreeModuleStructure`. A free module structure can be obtained from the ring of scalars by calling `.free_module(n)` (the module will take the scalar ring structure by reference) or `.into_free_module(n)` (the module will take ownership of the scalar ring structure).

For example, to obtain \\(\mathbb{Z}^3\\)

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::RingToFinitelyFreeModuleSignature;
# use algebraeon::sets::structure::MetaType;
#
let module = Integer::structure().into_free_module(3);
```

Elements of \\(\mathbb{Z}^3\\) are represented by objects of type `Vec<Integer>` and basic operations with the elements are provided by the module structure.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::rings::structure::*;
# use algebraeon::sets::structure::*;
# 
# let module = Integer::structure().into_free_module(3);
# 
let a = vec![1.into(), 2.into(), 3.into()];
let b = vec![(-1).into(), 2.into(), (-2).into()];

assert!(module.equal(
    &module.neg(&a),
    &vec![(-1).into(), (-2).into(), (-3).into()]
));

assert!(
    module.equal(&module.add(&a, &b), 
    &vec![0.into(), 4.into(), 1.into()]
));

assert!(
    module.equal(&module.sub(&a, &b),
    &vec![2.into(), 0.into(), 5.into()]
));

assert!(module.equal(
    &module.scalar_mul(&5.into(), &a),
    &vec![5.into(), 10.into(), 15.into()]
));
```

The scalar ring structure can be obtained from a module by calling `.ring()`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let module = Integer::structure().into_free_module(3);
# 
let ring = module.ring();
assert_eq!(ring, &Integer::structure());
```

## Submodules

The set of submodules of the free module \\(R^n\\) is represented by objects of type `FinitelyFreeSubmoduleStructure`. This structure can be obtained from a module by calling `.submodules()` (the structure will take the module structure by reference) or `.into_submodules()` (the structure will take ownership of the module structure).

For example, to obtain the set of all submodules of \\(\mathbb{Z}^3\\)

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
let submodules = Integer::structure().into_free_module(3).into_submodules();
```

Submodules of \\(R^n\\) are represented by objects of type `FinitelyFreeSubmodule`.

### Constructing Submodules

The zero submodule \\(\{0\} \subseteq R^n\\) can be constructed using `.zero_submodule()` and the full submodule \\(R^n \subseteq R^n\\) can be constructed using `.full_submodule()`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
assert_eq!(submodules.zero_submodule().rank(), 0);
assert_eq!(submodules.full_submodule().rank(), 3);
```

The submodule given by the span of some elements of the module can be constructed using `.span(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
assert_eq!(
    submodules
        .span(vec![
            &vec![1.into(), 2.into(), 2.into()],
            &vec![2.into(), 1.into(), 1.into()],
            &vec![3.into(), 3.into(), 3.into()]
        ])
        .rank(),
    2
);
```

The submodule given by the kernel of some elements can be constructed using `.kernel(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
assert!(submodules.equal(
    &submodules.kernel(vec![
        &vec![1.into(), 2.into()],
        &vec![2.into(), 1.into()],
        &vec![3.into(), 3.into()],
    ]),
    &submodules.span(vec![&vec![1.into(), 1.into(), (-1).into()]])
));
```

### Basic Operations

Test submodules for equality using `.equal(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
assert!(submodules.equal(
    &submodules.span(vec![
        &vec![2.into(), 2.into(), 0.into()],
        &vec![2.into(), (-2).into(), 0.into()],
    ]),
    &submodules.span(vec![
        &vec![4.into(), 0.into(), 0.into()],
        &vec![0.into(), 4.into(), 0.into()],
        &vec![2.into(), 2.into(), 0.into()],
    ])
));

assert!(!submodules.equal(
    &submodules.span(vec![
        &vec![1.into(), 1.into(), 0.into()],
        &vec![2.into(), 3.into(), 0.into()],
    ]),
    &submodules.span(vec![
        &vec![1.into(), 2.into(), 3.into()],
        &vec![1.into(), 1.into(), 0.into()],
    ])
));
```

Check whether a submodule contains an element using `.contains_element(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
let a = submodules.span(vec![
    &vec![2.into(), 2.into(), 0.into()],
    &vec![2.into(), (-2).into(), 0.into()],
]);

assert!(submodules.contains_element(&a, &vec![4.into(), 4.into(), 0.into()]));
assert!(!submodules.contains_element(&a, &vec![3.into(), 4.into(), 0.into()]));
assert!(!submodules.contains_element(&a, &vec![4.into(), 4.into(), 1.into()]));
```

Check whether a submodule is a subset of another submodule using `.contains(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
let a = submodules.span(vec![&vec![3.into(), 3.into(), 3.into()]]);
let b = submodules.span(vec![&vec![6.into(), 6.into(), 6.into()]]);

assert!(submodules.contains(&a, &b));
assert!(!submodules.contains(&b, &a));
```

Compute the sum of two submodules using `.sum(..)` and compute the intersection of two submodules using `.intersection(..)`.

```rust
# use algebraeon::nzq::Integer;
# use algebraeon::rings::linear::finitely_free_module::*;
# use algebraeon::sets::structure::*;
# 
# let submodules = Integer::structure().into_free_module(3).into_submodules();
# 
let a = submodules.span(vec![&vec![4.into(), 4.into(), 4.into()]]);
let b = submodules.span(vec![&vec![6.into(), 6.into(), 6.into()]]);

let sum_ab = submodules.span(vec![&vec![2.into(), 2.into(), 2.into()]]);
assert!(submodules.equal(&submodules.sum(&a, &b), &sum_ab));

let intersect_ab = submodules.span(vec![&vec![12.into(), 12.into(), 12.into()]]);
assert!(submodules.equal(&submodules.intersect(&a, &b), &intersect_ab));
```

<!-- ### Other Operations
 - Reducing an element (sometimes unique)
 - `extension_basis` -->

<!-- ## Cosets

### Constructing Cosets

 - From an element
 - Full coset
 - From a submodule
 - From a submodule and an offset

### Operations

 - Equality of cosets
 - Contains a point
 - Contains a cosets
 - Sum
 - Intersection

## Affine Subsets

### Constructing Affine Subsets

 - Empty
 - From an element
 - Full coset
 - From a submodule
 - From a submodule and an offset

### Operations

 - Equality
 - Contains a point
 - Contains an affine subset
 - Sum
 - Intersection
 -->

