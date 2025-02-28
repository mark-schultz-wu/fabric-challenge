#![doc = include_str!("../README.md")]
#![warn(clippy::all)]
#![warn(missing_docs)]
#![warn(rust_2018_idioms)]
#![warn(unused_imports)]
#![warn(unused_variables)]
#![warn(dead_code)]
#![warn(missing_debug_implementations)]

pub use crate::ff::{Field, MontgomeryFp};
pub use crate::poly::{
    GeneralMultivariatePolynomial, MultilinearPolynomial, MultivariatePolynomial,
    UnivariatePolynomial,
};
pub use crate::sumcheck::{ProtocolError, Prover, ProverMessage, Verifier, VerifierMessage};

mod ff;
mod poly;
mod sumcheck;
