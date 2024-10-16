//! [Open applied topology (OAT)](https://openappliedtopology.github.io) is a library for fast, user-friendly algebra and topology. OAT has 
//! 
//! - a user-friendly frontend for Python users, called [oat_python](https://crates.io/crates/oat_python)
//! - a fast backend written in Rust, called [oat_rust](https://crates.io/crates/oat_rust) 
//! - a variety of tutorials published as [jupyter notebooks](https://openappliedtopology.github.io)
//! 
//! This package contains the source code for [oat_rust](https://crates.io/crates/oat_rust).
//! 
//! # Contents
//! 
//! 
//! - [Welcome!](#welcome)
//! - [Community](#community)
//! - [Values](#values)
//! - [Mission](#mission)
//! - [Get started](#get-started)
//! 
//! # Welcome!
//! 
//! Welcome!  This package is [oat_rust](https://crates.io/crates/oat_rust), part of the Open Applied Topology ecosystem.  It provides powerful tools for applied topology, including
//! 
//! - [Algebra](crate::algebra)
//!   - [Rings](crate::algebra::rings)
//!   - [Vectors](crate::algebra::vectors)
//!   - [Matrices](crate::algebra::matrices)
//!   - [Chain complexes](crate::algebra::chains)
//!   - [Barcodes and persistence diagrams](crate::algebra::chains::barcode)
//!   - [(Persistent) (relative) (co)homology](crate::algebra::chains)
//!   - [Optimal cycle representatives](crate::algebra::chains::factored)
//! - [Topology](crate::topology)
//!   - [Simplicial complexes](crate::topology::simplicial)
//!   - [Boundary matrices](crate::topology::simplicial)
//!   - [Point clouds](crate::topology::point_cloud)
//! - [Utilities](crate::utilities)
//!   - [Iterators](crate::utilities::iterators)
//!   - [Optimization](crate::utilities::optimization)
//!   - [Order](crate::utilities::order)
//!   - [Sets](crate::utilities::sets)
//!   - [Indexing](crate::utilities::indexing)
//!   - and more!
//! - See our sister library, [oat_python](https://crates.io/crates/oat_python), for additional tools on
//!   - Persistent homology
//!   - Cycle optimization
//!   - Interactive 2d and 3d visualization
//! 
//! 
//! # Community
//! 
//! OAT is by and for the open source community.  Reach out to the developers at [openappliedtopology@gmail.com](mailto:openappliedtopology@gmail.com) if you
//! - Need help getting started
//! - Wish for a missing feature
//! - Want to try coding
//! 
//! A collaboration of 20 research centers at colleges, universities, private, and public organizations support OAT's
//! development. The founding developers are Princton University, Macalester College, and the University of Delaware
//! The National Science Foundation granted seed funding for OAT in
//! 2019 under the [ExHACT]((https://www.nsf.gov/awardsearch/showAward?AWD_ID=1854748&HistoricalAwards=false))
//! project, and [Pacific Northwest National Laboratory](https://www.pnnl.gov/) (PNNL) provides continuing financial
//! support.  PNNL now coordinates development.  See <span style="color: SteelBlue;">[here](./ATTRIBUTIONS.md)</span>
//! for further details.
//! 
//! # Values
//! 
//! Our [shared values](https://github.com/OpenAppliedTopology/oat_python/blob/main/CODE_OF_CONDUCT.md) are
//! 
//! - Inclusion
//! - Respect, and 
//! - A shared passion to expand human knowledge, through algebraic topology
//! 
//! 
//! # Mission
//! 
//! **Performance**
//!     
//! OAT is a first-class solver for cutting-edge applications.  It is ideally suited to large, sparse data sets.
//!     The core library is written in Rust, a low-level systems programming language designed for safety and performance.
//!     High-level wrappers are available in Python. 
//! 
//! 
//! **Reliability**
//! 
//! OAT has more unit tests than type definitions and function definitions, combined.
//!   Its modular design enables end users to write their own checks for correctness, with ease.
//!   The library inherits strong safety guarantees from the the Rust compiler.
//! 
//! 
//! **Transparency**
//! 
//! OAT documentation emphasizes clarity and accessibility for users with all backgrounds.  It includes more than 180 working examples, and describes both code and underlying mathematical concepts in detail. 
//! [Online Jupyter notebook tutorials](crate::tutorials) illustrate how to combine multiple tools into larger applications.
//! The platform's modular design breaks large solvers into basic components, which it exposes to the user for inspection.  In addition, the library provides powerful methods to inspect and analyze objects, consistent with the way humans naturally think about problems; for example, you can look up rows and columns of boundary matrices using *cubes*, *simplices*, or *cells* as keys.
//!   
//! 
//! **Modularity**
//! 
//!   OAT reduces complex problems to the same basic building blocks that topologists use when writing on a chalk board. Users can mix and match those blocks with confidence, using a simple, streamlined interface.  They can even create new components that work seemlessly with the rest of the library, including coefficient rings, sparse matrix data structures, and customized filtrations on simplicial complexes.
//! 
//! 
//! # Get Started
//! 
//! 
//! **Python users**
//! 
//! If you'd like to use OAT for a project coded in Python, check out the [oat-python package](https://pypi.org/project/oat-python/)! You can also find [Jupyter notebook tutorials](https://openappliedtopology.github.io) on the OAT homepage.
//! 
//! **First time Rust users** 
//! 
//! If you're new to Rust, there are a few good resources:
//! 
//! - [The Rust Book](https://doc.rust-lang.org/book/ch00-00-introduction.html) is written with begginers in mind
//! - [Rust By Example](https://doc.rust-lang.org/stable/rust-by-example/) complements the explanations in the Rust book with examples -- it's often much easier to grasp the concepts.
//! - [Tutorials](crate::tutorials) is part of the documenation  for oat_rust. Because the Rust Book and Rust By Example are fairly long, it can be hard to get a birds-eye-view of the development process. The tutorials page is designed to fill in some gaps, and also clarify a few of the finer points of coding in Rust.
//! - [introduction_to_rust](https://github.com/pnnl/introduction_to_rust_pnnl) is a small Rust package that shows how to perform several common tasks, like creating a new package from scratch, building and viewing the documentation, and getting dependencies.
//! 
//! Finally, we recommend
//! 
//! - [Visual Studio Code](https://code.visualstudio.com) with the [Rust Analyzer Extension](https://rust-analyzer.github.io) for writing your code. They are simple to use, and highly convenient.








pub mod algebra;
pub mod topology;
pub mod utilities;
pub mod tutorials;
