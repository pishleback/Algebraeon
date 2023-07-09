#![allow(dead_code)]

use rayon::prelude::*;

use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fmt::Debug,
    hash::Hash,
    borrow::Borrow,
};

use crate::todd_coxeter;

include!("group.rs");
include!("subset.rs");
include!("subgroup.rs");
include!("normal_subgroup.rs");
include!("generating_set.rs");
include!("partition.rs");
include!("homomorphism.rs");
include!("iso_rep.rs");