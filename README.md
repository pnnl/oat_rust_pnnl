# About this fork

This is a fork of the [oat_rust](https://github.com/OpenAppliedTopology/oat_rust) repository. The fork is maintained by Pacific Northwest National Laboratory (PNNL), which is the lead contributor to the Open Applied Topology project.

# Open Applied Topology

[Open applied topology (OAT)](https://openappliedtopology.github.io) is a library for fast, user-friendly algebra and topology. OAT has 

- a user-friendly frontend for Python users, called [oat_python](https://github.com/OpenAppliedTopology/oat_python)
- a fast backend written in Rust, called [oat_rust](https://github.com/OpenAppliedTopology/oat_rust) 
- a variety of tutorials published as [jupyter notebooks](https://openappliedtopology.github.io)

This package contains the source code for [oat_rust](https://github.com/OpenAppliedTopology/oat_rust).


# Caution: breaking changes

OAT is in early stages of develpoment, and it's evolving quickly. Code that you write today may not work tomorrow, due to these changes. We will do our very best to make sure that if/when this happens, you will only need to make small changes to your code to fix the problem (e.g., updating the name of a function). However, please do bear this in mind as you write your code!

**These efforts notwithstanding, expect major changes to oat_rust in late 2024/early 2025.**

# Install and explore

**Python users** 

If you'd like to use OAT as a Python user, you'll want our sister library, oat_python, which is available on [PyPI](https://pypi.org/project/oat_python/). 

<!-- You don't need to have Rust installed to use oat_python; just use your favorite package manager, for example `pip install oat_python`, `conda_install oat_python` etc.!  Explore the [Jupyter notebook tutorials on GitHub](https://github.com/OpenAppliedTopology/oat)! -->

**Rust users** 

The oat_rust package is a pure Rust library. Rust libraries aren't installed in the same was as Python packages. Instead, users typically use Rust packages in one of two ways: (1) As a dependency for other Rust packages. In this case you'll have a folder containing code for another package, including a file named `Cargo.toml`. Add oat_rust to the list of dependencies in that file, and Rust will automatically download and use oat_rust when you compile your project. Here's an [example Cargo.toml](https://github.com/OpenAppliedTopology/oat_python/blob/main/Cargo.toml). The relevant line in that file is `oat_rust = "X.Y.Z"`.  Here X.Y.Z refers to a specific version of oat_rust. You'll want to check [Crates.io](https://crates.io/crates/oat_rust) to find the version number for the most recent version of oat. [Rust By Example](https://doc.rust-lang.org/rust-by-example/cargo/deps.html) gives a great explanation of dependencies! (2) To compile executable files. This allows you to write programs which can be executed even if you don't have Rust installed!  You'll need Rust to compile those programs first, though. Check out this example from [Rust By Example](https://doc.rust-lang.org/stable/rust-by-example/hello.html?highlight=executable#hello-world) to get started!



# Documentation

There are two ways to view the documenation for this package.

**Web browser** 

Browse the documentation on [Crates.io](https://docs.rs/oat_rust).

**Build from source** 

Most users won't need this option, but if you are modifying and extending oat_rust, then you can use it to access the documenation for your extension. First [install Rust](https://www.rust-lang.org/tools/install).  If you already have Rust, make sure you have the most recent version.  Next, obtain a copy of the the [oat_rust git repository](https://github.com/OpenAppliedTopology/oat_rust). Open a shell, CD into the repository, and run `cargo doc --no-deps --open`.  The documentaiton homepage will then open in a browser.

# Contributing

For information on **contributing**, see [`CONTRIBUTING.md`](https://github.com/OpenAppliedTopology/oat_python/blob/main/CONTRIBUTING.md).

# License

For inforamtion on copyright and licensing, see [`LICENSE`](https://github.com/OpenAppliedTopology/oat_python/blob/main/LICENSE).

# Attributions

OAT is an extension of the ExHACT library. See [`ATTRIBUTIONS.md`](https://github.com/OpenAppliedTopology/oat_python/blob/main/ATTRIBUTIONS.md) for details.

<!-- # Python Installation from source

1. Download and install the most recent version of [Rust](https://www.rust-lang.org/).  Make sure your installation is up to date by running `rustup update` in a command shell.

2. Create a virtual Python environment, e.g. using Anaconda, or open one that you already have.  In this example, let's assume the environment name is `myenv`.  Activate `myenv`, and run

    ```bash
    pip install maturin
    ```

    A number of warning messages may appear; this is normal. 

    **If you already have maturin installed, make sure it is up to date!**

3. [Clone](https://github.com/OpenAppliedTopology/oat_python) a copy of oat_python. Open a shell and CD into the oat_python folder.  Activate `myenv` and run

    ```bash
    maturin develop --release
    ```
    
5. oat_python should now be installed.  Try running the Jupyter notebooks with `myenv`! -->




<!-- # Open Applied Topology

Open applied topology (OAT) is a library for fast, user-friendly algebra and topology.

## Documentation

There are two ways to view the documenation for this package.

**Web browser** 

Select `Documenation` from the menu on the righthand side of the page at [Crates.io](https://crates.io/crates/oat_rust)

**Build from source** 

First [install Rust](https://www.rust-lang.org/tools/install).  If you already have Rust, make sure you have the most recent version.  Next, obtain a copy of the the [oat_rust git repository](https://github.com/OpenAppliedTopology/oat_rust). Open a shell, CD into the repository, and run `cargo doc --no-deps --open`.  The documentaiton homepage will then open in a browser.  Documentation includes



## Installation

Installation of this package is managed by the Rust feature `Cargo`.  

- See the [`Cargo`](https://doc.rust-lang.org/cargo/) documentation for details on compiling and running the code.
- See the [`introduction_to_rust`](https://github.com/OpenAppliedTopology/introduction_to_rust) crate for a gentle introduction to programming in Rust.


## Contributing

For information on contributing, see [`CONTRIBUTING.md`](https://github.com/OpenAppliedTopology/oat_rust/blob/main/CONTRIBUTING).

## License and attributions

For inforamtion on copyright and licensing, see [`LICENSE`](https://github.com/OpenAppliedTopology/oat_rust/blob/main/LICENSE). 

OAT is an extension of the ExHACT project, with major contributions from Pacific Northwest National Laboratory. For details see [`ATTRIBUTIONS.md`](https://github.com/OpenAppliedTopology/oat_rust/blob/main/ATTRIBUTIONS.md). -->
