//! Chain complexes.
//! 
//! This module provide tools for chain complexes. 
//! 
//! - [Examples](#example-homology-calculation)
//! - [Tools](#tools)
//! - [Background](#background)
//! 
//! # Example homology calculations
//! 
//! See the [Background](#background) section below for a description of the overall approach we use to compute (persistent) homology in OAT. This isn't the only approach -- not even the only approach available in OAT!  But it's one of the most convenient.
//! 
//! - [Calculate the homology of a simplicial complex, with generators](crate::topology::simplicial::from::relation) 
//! - [Calculate the barcode of a filtered Vietoris-Rips complex, with generators](crate::topology::simplicial::from::graph_weighted::ChainComplexVrFiltered)
//! 
//! # Tools
//! 
//! This module provides the following submodules. See [Background](#background) for an explanation of terms.
//! 
//! - [factored]
//!   - tools for the U-match decomposition JM = DJ
//!   - a [data structure called FactoredBoundaryMatrix](factored::FactoredBoundaryMatrix) that stores the U-match decomposition. This data structure was a wide range of methods for enumerating subsets of indices and basis vectors (birth simplices, death simplices, etc).
//!   - a [function](factored::factor_boundary_matrix), which will factor the differential matrix D and return a [FactoredBoundaryMatrix](factored::FactoredBoundaryMatrix)
//!   - additional tools for working with the matrix J. 
//! 
//! - [jordan]
//!   - tools for the matrix J
//!   - a data structure for the matrix J, [jordan::JordanBasisMatrix]
//! 
//! - [barcode]
//!   - handy utilities for barcodes
//!   - a function to [extract the barcode](barcode::barcode) from a U-match decomposition
//!   - data structures for bars and barcodes, with a variety of associated methods, e.g. to compute
//!     - betti curves
//!     - left- and right-endpoints
//!     - grouping by dimension
//!     - iterators
//!     - etc.
//! 
//! 
//! 
//! # Background
//! 
//! ## The differential matrix
//! 
//! The differential matrix of a chain complex with a fixed basis is a matrix D with one row and column for each basis vector in the complex. Column j of D represents the boundary of the jth basis vector. This matrix is often called the boundary matrix (this documentation does that too), though the term boundary matrix often refers to just a part of the differential matrix indexed by chains of a given dimension, e.g. the boundaries of just the triangles in a simplicial complex.
//! 
//! Often we are interested in the chain complex of a simplicial complex; in that case the differential matrix has a row and column for each simplex of every dimension.
//! 
//! 
//! ## Homology, persistent homology, and U-match
//! 
//! To compute the (persistent) homology of a chain complex C, OAT uses [U-match decomposition](https://arxiv.org/abs/2108.08831).
//! 
//! Specifically, we obtain an upper triangular matrix J with 1's on the diagonal, and a generalized matching matrix M, such that 
//! 
//! JM = DJ
//! 
//! where D is the differential matrix of the chain complex.  
//! 
//! Provided that the rows and colums of D are sorted in ascending order, we can then use the matrix J to compute the homology of C, and the *persistent* homology of C, if C is a filtered chain complex.  
//! 
//! **If C is a filtered complex**
//! 
//! * **Finite bars** The barcode has one interval of form `[birth(s), birth(t))` for every `d`-simplex `s` 
//! and `d+1`-simplex `t` such that `M[s,t]` is nonzero.  Here `birth(s)` denotes the time when
//! `s` enters the filtration.
//! The corresponding basis vector is column `s` of the matrix `R`.  We call `s` and `t` birth
//! and death simplices.  Usually we exclude empty intervals.
//! 
//! * **Infinite bars**  The barcode has one interval of form `[birth(s), inf)` for every `d`-simplex `s`
//! such that `M[s,:]=0` and `M[:,s]=0`.  The corresponding basis vector is column `s` of `C`.
//! 
//! **If C is an unfiltered complex**
//! 
//! In this case we can still treat C as a filtered complex, where the filtration is trivial (it consists of C and no subcomplexes). In that case we can still compute the barcode, and the cycle representatives for the infinite bars provide a basis for the homology of C.

pub mod factored;
pub mod jordan;
pub mod barcode;