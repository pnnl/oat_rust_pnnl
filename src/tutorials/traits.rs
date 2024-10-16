//! # A simple introduction to traits in Rust
//! 
//! Also check out the  [introduction to traits in the Rust book!](https://doc.rust-lang.org/book/ch10-02-traits.html)
//! 
//! Let’s start with a tiny example that demonstrates how traits allow us to write flexible and reusable functions. Below is a function that accepts any type T that can be converted into a usize. We express this constraint using the From trait:
//! 
//! ```
//! fn eat_and_print<T>(item: T)
//! where
//!     usize: From<T>,
//! {
//!     let converted: usize = usize::from(item);
//!     println!("{}", converted);
//! }
//! ```
//! 
//! Explanation:
//! 
//! - `T` is called a generic type. It's just a symbol that acts as a placeholder for other types, the same way that we use `x` as a placeholder in the statement "solve 2x = 1".
//! - `where usize: From<T>` specifies that there must be an implementation of the From trait for usize from the type T. This means that usize can be created from a value of type T.
//! - The function converts the input to usize using `usize::from(item)` and prints the result.
//! 
//! **Example usage**
//! 
//! ```
//! fn eat_and_print<T>(item: T)
//! where
//!     usize: From<T>,
//! {
//!     let converted: usize = usize::from(item);
//!     println!("{}", converted);
//! }
//! 
//! fn main() {
//!     eat_and_print(42u8);  // u8 implements From<u8> for usize
//!     eat_and_print(5u16);  // u32 implements From<u32> for usize
//! }
//! ```
//! 
//! # A more complex example
//! 
//! Next, let’s create another function, `add_one``, which works with any type T that implements both `From<usize>`` and the `Add`` trait (for adding values). This function will consume the input, add one to it, and return the result.
//! 
//! **Function Example: `add_one`**
//! 
//! 
//! ```
//! fn add_one<T>(item: T) -> T
//! where
//!     T: From<usize> + std::ops::Add<usize, Output = T>,
//! {
//!     item + 1
//! }
//! ```
//! 
//! **Explanation:**
//! 
//! - `T: From<usize> + std::ops::Add<usize>` means that T must implement both the `From<usize>` trait (so that it can be created from a usize) and the `Add` trait (so that we can add a `usize` to it).
//! - The function consumes the input item, adds one to it, and returns the result.
//! 
//! **Example usage:**
//! 
//! ```
//! fn add_one<T>(item: T) -> T
//! where
//!     T: From<usize> + std::ops::Add<usize, Output = T>,
//! {
//!     item + 1
//! }
//! 
//! fn main() {
//!     let result = add_one(41usize);
//!     println!("{}", result);  // Outputs: 42
//! }
//! ```
//! 
//! Here, `add_one`` works with any type that supports both being constructed from a `usize` and adding a `usize` to it.
//! 
//! # Make a trait available for a new struct
//! 
//! Let's create a simple struct `CookieJar` that holds the number of cookies. We'll implement the Add trait for it so that we can easily combine two `CookieJar` instances by adding their cookie counts together.
//! 
//! ```
//! struct CookieJar {
//!     number_of_cookies: usize,
//! }
//! ```
//! 
//! Next, we’ll implement the `Add` trait for `CookieJar`:
//! 
//! 
//! ```
//! use std::ops::Add;
//! 
//! struct CookieJar {
//!     number_of_cookies: usize,
//! }
//! 
//! impl Add for CookieJar {
//!     type Output = CookieJar;
//! 
//!     fn add(self, other: CookieJar) -> CookieJar {
//!         CookieJar {
//!             number_of_cookies: self.number_of_cookies + other.number_of_cookies,
//!         }
//!     }
//! }
//! ```
//! 
//! **Example usage:**
//! 
//! ```
//! use std::ops::Add;
//! 
//! struct CookieJar {
//!     number_of_cookies: usize,
//! }
//! 
//! impl Add for CookieJar {
//!     type Output = CookieJar;
//! 
//!     fn add(self, other: CookieJar) -> CookieJar {
//!         CookieJar {
//!             number_of_cookies: self.number_of_cookies + other.number_of_cookies,
//!         }
//!     }
//! }
//! 
//! fn main() {
//!     let jar1 = CookieJar { number_of_cookies: 5 };
//!     let jar2 = CookieJar { number_of_cookies: 3 };
//!     
//!     let combined = jar1 + jar2;
//!     println!("Total cookies: {}", combined.number_of_cookies);  // Outputs: 8
//! }
//! ```
//! 
//! # Conclusion
//! - Traits allow us to define shared behavior and write reusable, flexible code.
//! - In the first example, we constrained a function to work with types that can be converted to usize using the From trait.
//! - In the second example, we built a function (add_one) that works on any type implementing both From<usize> and Add.
//! - Finally, we defined a custom struct (CookieJar) and implemented the Add trait for it to allow cookie jars to be combined.
//! By leveraging traits, we gain the ability to generalize behavior across types, making Rust code both powerful and highly reusable!









//! Gentle introduction to traits
//! 
//! One way to sort a list `L = [1,2,3]` in Python is to call `L.sort()`.
//! The function `.sort()` is an example of a *method*, which many
//! programming languages support.
//! 
//! Most languages treat the methods defined on two different objects as
//! different, even if they have the same name.  For example the `.sort()`
//! method used on lists is technically different from the `.sort()` method 
//! used on Numpy arrays.  However, sometimes it can be useful to write code
//! that will work with *any* object that has a `.sort()` method.
//! 
//! This is where traits come in.  A `trait` in Rust essentially consists of
//! - a function name, e.g. `my_method`, and 
//! - a function signature, which specifies the inputs that the function `my_method` should take
//! and the outputs that it should return.
//!
//! After defining a trait, you can define a method with that trait name
//! on multiple different objects -- as long as the method has the correct
//! function signature.  *Then you can write a function `F` that takes any
//! input `I` such that `I` has a method called `my_method`.*
//! 
//! This becomes important in OAT because the most basic operations on matrices, vectors, and
//! vector entries are all governed by trait methods.  For example, looking up a 
//! row or column of matrix is a trait method.
//! 
