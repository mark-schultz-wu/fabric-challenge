//! TODO: Fill in documentation

#![warn(clippy::all)]
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(unused_imports)]
#![warn(unused_variables)]
#![warn(dead_code)]
#![warn(missing_debug_implementations)]

mod ff;
pub use ff::*;
mod poly;

fn main() {
    println!("Hello, world!");
}
