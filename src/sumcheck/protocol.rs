//! The prover in the sumcheck protocol

use crate::Field;
use crate::MultivariatePolynomial;
use crate::UnivariatePolynomial;
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::marker::PhantomData;

/// The errors of the Sumcheck protocol
#[derive(Debug, Clone, Copy)]
pub enum ProtocolError {
    ProtocolAlreadyCompleted,
    ProcessMessageInWrongOrder,
    AttemptToShrinkConstantPoly,
}

/// The prover in the sumcheck protocol
#[derive(Debug, Clone)]
pub struct Prover<F: Field, P: MultivariatePolynomial<F>> {
    /// A multivariate polynomial, which initially is `g`, the polynomial sumcheck is being applied to.
    /// Prover modifies this polynomial during each round though
    shrunk_g: P,
    _field_data: PhantomData<F>,
}

#[derive(Debug, Clone)]
pub enum ProverMessage<F: Field> {
    InitialMessage(F, UnivariatePolynomial<F>),
    OtherMessage(UnivariatePolynomial<F>),
}

impl<F: Field, Poly: MultivariatePolynomial<F>> Prover<F, Poly> {
    /// Constructs a new prover
    pub fn new(g: Poly) -> Self {
        Self {
            shrunk_g: g,
            _field_data: PhantomData,
        }
    }
    /// A message the prover sends during the protocol.
    pub fn message(
        &mut self,
        message: VerifierMessage<F>,
    ) -> Result<ProverMessage<F>, ProtocolError> {
        match message {
            // Protocol execution is over, prover should not be processing this message
            VerifierMessage::Reject | VerifierMessage::Accept => {
                Err(ProtocolError::ProtocolAlreadyCompleted)
            }
            VerifierMessage::Initial => {
                // At the start of the protocol
                // 1. compute the sum of g over the boolean hypercube, and
                // 2. take a univariate slice of g, and
                // transmit these to the verifier
                //
                let sum = self.shrunk_g.sum_over_boolean_hypercube();
                let slice = self.shrunk_g.univariate_slice_last();
                Ok(ProverMessage::InitialMessage(sum, slice))
            }
            VerifierMessage::Challenge(f) => {
                // Map the ShrinkError to a ProtocolError
                self.shrunk_g
                    .shrink_last(&f)
                    .map_err(|_| ProtocolError::AttemptToShrinkConstantPoly)?;

                // If we get here, shrinking succeeded, so return the slice
                Ok(ProverMessage::OtherMessage(
                    self.shrunk_g.univariate_slice_last(),
                ))
            }
        }
    }
}

/// The verifier in the sumcheck protocol
#[derive(Debug, Clone)]
pub struct Verifier<F: Field, P: MultivariatePolynomial<F>> {
    /// The polynomial that the protocol is being performed on
    g: P,
    /// The claimed value of \sum_{0,1} g that the prover is attempting to prove
    sum_val: Option<F>,
    /// The list of challenges that the verifier is performing
    challenges: Vec<F>,
    /// The RNG that `Verifier` uses
    rng: StdRng,
}

#[derive(Debug, Clone)]
/// The messages the verifier sends during sumcheck
pub enum VerifierMessage<F: Field> {
    Challenge(F),
    Initial,
    Accept,
    Reject,
}

impl<F: Field, Poly: MultivariatePolynomial<F>> Verifier<F, Poly> {
    /// Constructs a new verifier
    pub fn new(g: Poly, seed: Option<u64>) -> Self {
        let rng = if let Some(s) = seed {
            StdRng::seed_from_u64(s)
        } else {
            StdRng::from_os_rng()
        };
        let challenges = Vec::with_capacity(g.num_variables());
        Self {
            g,
            sum_val: None,
            challenges,
            rng,
        }
    }

    /// The verifier's message in the protocol.
    /// They take as input a univariate slice,
    /// check that
    ///     1. s(0) + s(1) = C, the claimed value of the sum, and
    ///     2. deg(s) <= deg_i(g), where i is the round #.
    /// They then output a randomly chosen field element to the prover,
    /// or Accept/Reject in the final round.
    pub fn message(
        &mut self,
        message: ProverMessage<F>,
    ) -> Result<VerifierMessage<F>, ProtocolError> {
        let slice = match message {
            ProverMessage::InitialMessage(sum_val, slice) => {
                self.sum_val = Some(sum_val);
                slice
            }
            ProverMessage::OtherMessage(slice) => slice,
        };

        // The index of the variable we are processing
        let variable_number = self.g.num_variables() - self.challenges.len();
        // self.sum_val is only `None` if we have not processed
        // the initial message
        let sum_val: &F = self
            .sum_val
            .as_ref()
            .ok_or(ProtocolError::ProcessMessageInWrongOrder)?;

        // Check sum constraint
        if slice.evaluate(F::zero()) + slice.evaluate(F::one()) != *sum_val {
            // The sum check fails
            return Ok(VerifierMessage::Reject);
        }

        // Get the degree of the univariate slice
        // If it's None, that means we have a zero polynomial, which should have degree 0
        let slice_degree = slice.degree().unwrap_or(0);

        // Get the degree bound for this round. Again, it is `None` for the zero polynomial.
        let degree_bound = self.g.degree(variable_number).unwrap_or(0);

        // Check degree constraint
        if slice_degree > degree_bound {
            return Ok(VerifierMessage::Reject);
        }

        // Check if this is the final round
        if variable_number == 1 {
            // Last variable
            // Generate the final challenge
            let final_challenge = F::random(&mut self.rng);

            // Store the final challenge
            self.challenges.push(final_challenge.clone());

            // Evaluate the original polynomial at all challenge points
            let expected = self.g.evaluate(&self.challenges);

            // Evaluate the final univariate slice at the final challenge
            let actual = slice.evaluate(final_challenge);

            // Check if the values match
            if expected == actual {
                return Ok(VerifierMessage::Accept);
            } else {
                return Ok(VerifierMessage::Reject);
            }
        }

        // Generate and return challenge for non-final rounds
        let challenge = F::random(&mut self.rng);
        self.challenges.push(challenge.clone());
        Ok(VerifierMessage::Challenge(challenge))
    }
}
