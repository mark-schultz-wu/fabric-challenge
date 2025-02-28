//! The prover in the sumcheck protocol
#![allow(dead_code)]

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

/// The message variants the prover sends during the sumcheck protoocl
#[derive(Debug, Clone, PartialEq, Eq)]
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
            VerifierMessage::Reject(_) | VerifierMessage::Accept => {
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
    hypercube_sum: Option<F>,
    /// The value of the previous-round's challenge evaluated on the previous-round's univariate slice
    slice_challenge_sum: Option<F>,
    /// The list of challenges that the verifier is performing
    challenges: Vec<F>,
    /// The RNG that `Verifier` uses
    rng: StdRng,
}

#[derive(Debug, Clone, PartialEq, Eq)]
/// The messages the verifier sends during sumcheck
pub enum VerifierMessage<F: Field> {
    Challenge(F),
    Initial,
    Accept,
    Reject(String),
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
            hypercube_sum: None,
            slice_challenge_sum: None,
            challenges,
            rng,
        }
    }

    /// The verifier's message in the protocol.
    /// They take as input a univariate slice,
    /// check that
    ///     1. g_j(0) + g_j(1) = g_{j-1}(r_{j-1}), where r_{j-1} is the previous challenge
    ///     2. deg(g_j) <= deg_j(g), where j is the round #.
    /// They then output a randomly chosen field element to the prover,
    /// or Accept/Reject in the final round.
    pub fn message(
        &mut self,
        message: ProverMessage<F>,
    ) -> Result<VerifierMessage<F>, ProtocolError> {
        let slice = match message {
            ProverMessage::InitialMessage(hypercube_sum, slice) => {
                self.hypercube_sum = Some(hypercube_sum.clone());
                self.slice_challenge_sum = Some(hypercube_sum);
                slice
            }
            ProverMessage::OtherMessage(slice) => slice,
        };
        // The index of the variable we are processing
        let variable_number = self.g.num_variables() - self.challenges.len() - 1;
        // self.slice_challenge_sum is only `None` in the first round
        let slice_challenge_sum: &F = self
            .slice_challenge_sum
            .as_ref()
            .ok_or(ProtocolError::ProcessMessageInWrongOrder)?;

        // Get the degree of the univariate slice
        // If it's None, that means we have a zero polynomial, which should have degree 0
        let slice_degree = slice.degree().unwrap_or(0);

        // Get the degree bound for this round. Again, it is `None` for the zero polynomial.
        let degree_bound = self.g.degree(variable_number);

        // Check sum constraint
        if slice.evaluate(&F::zero()) + slice.evaluate(&F::one()) != *slice_challenge_sum {
            return Ok(VerifierMessage::Reject(String::from("Sum Check Violation")));
        }

        // Check degree constraint
        if slice_degree > degree_bound {
            return Ok(VerifierMessage::Reject(String::from(
                "Degree Bound Violation",
            )));
        }

        // Check if this is the final round
        if variable_number == 0 {
            // Last variable
            // Generate the final challenge
            let final_challenge = F::random(&mut self.rng);

            // Store the final challenge
            self.challenges.push(final_challenge.clone());

            // Reverses the vector so that the challenge associated with the first variable has the smallest index
            self.challenges.reverse();
            // Evaluate the original polynomial at all challenge points
            let expected = self.g.evaluate(&self.challenges);

            // Evaluate the final univariate slice at the final challenge
            let actual = slice.evaluate(&final_challenge);

            // Check if the values match
            if expected == actual {
                return Ok(VerifierMessage::Accept);
            } else {
                return Ok(VerifierMessage::Reject(String::from(
                    "Final Evaluation Failed",
                )));
            }
        }

        // Generate and return challenge for non-final rounds
        let challenge = F::random(&mut self.rng);
        self.challenges.push(challenge.clone());
        self.slice_challenge_sum = Some(slice.evaluate(&challenge));
        Ok(VerifierMessage::Challenge(challenge))
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{GeneralMultivariatePolynomial, MontgomeryFp, MultilinearPolynomial};
    use bitvec::prelude::*;
    use std::collections::HashMap;

    // Using a small prime (251) for testing
    type F = MontgomeryFp<251>;

    // Helper function to create a polynomial for testing
    fn create_test_polynomial() -> GeneralMultivariatePolynomial<F> {
        // Create a simple polynomial: f(x,y,z) = 5 + 3xy + 2yz^2
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term

        GeneralMultivariatePolynomial::from_coefficients(coeffs)
    }

    // Helper function to create a multilinear polynomial for testing
    fn create_multilinear_test_polynomial() -> MultilinearPolynomial<F> {
        // Create a simple multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz
        let mut coeffs = HashMap::new();
        let mut bv = BitVec::new();

        // [0, 0, 0] - Constant term
        bv.clear();
        bv.push(false);
        bv.push(false);
        bv.push(false);
        coeffs.insert(bv.clone(), F::from(5));

        // [1, 1, 0] - xy term
        bv.clear();
        bv.push(true);
        bv.push(true);
        bv.push(false);
        coeffs.insert(bv.clone(), F::from(3));

        // [0, 1, 1] - yz term
        bv.clear();
        bv.push(false);
        bv.push(true);
        bv.push(true);
        coeffs.insert(bv.clone(), F::from(2));

        MultilinearPolynomial::from_coefficients(coeffs)
    }

    #[test]
    fn test_successful_protocol_execution() {
        // Test with a general multivariate polynomial
        let poly = create_test_polynomial();

        // Create prover and verifier
        let mut prover = Prover::new(poly.clone());
        let mut verifier = Verifier::new(poly, Some(42)); // Use a fixed seed for deterministic testing

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");
        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // Run the protocol rounds
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            let prover_msg = prover
                .message(VerifierMessage::Challenge(challenge))
                .expect("Failed to get prover message");
            current_verifier_msg = verifier
                .message(prover_msg)
                .expect("Failed to process prover message");
        }

        // Check that the verifier accepted
        assert_eq!(current_verifier_msg, VerifierMessage::Accept);
    }

    #[test]
    fn test_multilinear_protocol_execution() {
        // Test with a multilinear polynomial
        let poly = create_multilinear_test_polynomial();

        // Create prover and verifier
        let mut prover = Prover::new(poly.clone());
        let mut verifier = Verifier::new(poly, Some(42)); // Use a fixed seed for deterministic testing

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");
        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // Run the protocol rounds
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            let prover_msg = prover
                .message(VerifierMessage::Challenge(challenge))
                .expect("Failed to get prover message");
            current_verifier_msg = verifier
                .message(prover_msg)
                .expect("Failed to process prover message");
        }

        // Check that the verifier accepted
        assert_eq!(current_verifier_msg, VerifierMessage::Accept);
    }

    #[test]
    fn test_protocol_with_dishonest_prover() {
        let poly = create_test_polynomial();

        // Create a dishonest prover by modifying the polynomial
        let mut dishonest_poly = poly.clone();
        // Add a term that wasn't in the original polynomial
        dishonest_poly.set_coefficient(vec![2, 0, 0], F::from(10));

        let mut prover = Prover::new(dishonest_poly);
        let mut verifier = Verifier::new(poly, Some(42)); // Verifier uses the original polynomial

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");
        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // Run the protocol rounds
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            match prover.message(VerifierMessage::Challenge(challenge)) {
                Ok(prover_msg) => {
                    current_verifier_msg = verifier
                        .message(prover_msg)
                        .expect("Failed to process prover message");

                    // If the verifier rejected, we're done
                    if let VerifierMessage::Reject(_) = current_verifier_msg {
                        break;
                    }
                }
                Err(_) => {
                    // If the prover errors, we'll count that as a rejection
                    current_verifier_msg =
                        VerifierMessage::Reject(String::from("The prover has errored"));
                    break;
                }
            }
        }
        // Check that the verifier rejected
        match current_verifier_msg {
            VerifierMessage::Reject(_) => {
                // Test passes - we got a rejection as expected
            }
            _ => {
                panic!(
                    "Expected a Reject message but got {:?}",
                    current_verifier_msg
                );
            }
        }
    }

    #[test]
    fn test_protocol_error_handling() {
        let poly = create_test_polynomial();
        let mut prover = Prover::new(poly.clone());

        // Test handling of completed protocol
        let result = prover.message(VerifierMessage::Accept);
        assert!(matches!(
            result,
            Err(ProtocolError::ProtocolAlreadyCompleted)
        ));
        let result = prover.message(VerifierMessage::Reject(String::from("Test rejection")));
        assert!(matches!(
            result,
            Err(ProtocolError::ProtocolAlreadyCompleted)
        ));

        // Test handling of wrong message order in verifier
        let mut verifier = Verifier::new(poly, Some(42));
        let result = verifier.message(ProverMessage::OtherMessage(UnivariatePolynomial::new(
            vec![F::from(1)],
        )));
        assert!(matches!(
            result,
            Err(ProtocolError::ProcessMessageInWrongOrder)
        ));
    }

    #[test]
    fn test_shrink_constant_polynomial() {
        // Create a constant polynomial with no variables
        let constant_poly = GeneralMultivariatePolynomial::<F>::zero(0);
        let mut prover = Prover::new(constant_poly);

        // Try to shrink a constant polynomial
        let challenge = F::from(5);
        let result = prover.message(VerifierMessage::Challenge(challenge));

        // Should error with AttemptToShrinkConstantPoly
        assert!(matches!(
            result,
            Err(ProtocolError::AttemptToShrinkConstantPoly)
        ));
    }

    #[test]
    fn test_protocol_with_zero_polynomial() {
        // Test with a zero polynomial
        let zero_poly = GeneralMultivariatePolynomial::<F>::zero(3);
        let sum = zero_poly.sum_over_boolean_hypercube();
        // Sum should be zero
        assert_eq!(&sum, &F::zero());

        let mut prover = Prover::new(zero_poly.clone());
        let mut verifier = Verifier::new(zero_poly, Some(42));

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");

        // Extract the sum value and check
        if let ProverMessage::InitialMessage(claimed_sum, _) = prover_msg {
            assert_eq!(claimed_sum, F::zero());
        } else {
            panic!("Unexpected prover message");
        }

        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // Run the protocol rounds
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            let prover_msg = prover
                .message(VerifierMessage::Challenge(challenge))
                .expect("Failed to get prover message");
            current_verifier_msg = verifier
                .message(prover_msg)
                .expect("Failed to process prover message");
        }

        // Check that the verifier accepted
        assert_eq!(current_verifier_msg, VerifierMessage::Accept);
    }

    #[test]
    fn test_protocol_with_variable_number_mismatch() {
        // Create polynomials with different numbers of variables
        let poly1 = create_test_polynomial(); // 3 variables

        // Create a polynomial with 4 variables
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0, 0], F::from(5));
        coeffs.insert(vec![1, 1, 0, 0], F::from(3));
        let poly2 = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Create prover with poly1 and verifier with poly2
        let mut prover = Prover::new(poly1);
        let mut verifier = Verifier::new(poly2, Some(42));

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");
        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // Run the protocol rounds until error or rejection
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            match prover.message(VerifierMessage::Challenge(challenge)) {
                Ok(prover_msg) => {
                    current_verifier_msg = verifier
                        .message(prover_msg)
                        .expect("Failed to process prover message");

                    // If the verifier rejected, we're done
                    if let VerifierMessage::Reject(_) = current_verifier_msg {
                        break;
                    }
                }
                Err(_) => {
                    // If the prover errors, we'll count that as a rejection
                    current_verifier_msg = VerifierMessage::Reject(String::from("Prover Errored"));
                    break;
                }
            }
        }

        // Check that the protocol did not succeed
        // Either rejected or errored due to the mismatch
        assert_ne!(current_verifier_msg, VerifierMessage::Accept);
    }

    #[test]
    fn test_check_degree_constraint() {
        // Create a polynomial with degree higher than claimed
        let honest_poly = create_test_polynomial();

        // Create a modified copy for the prover that has a higher degree
        let mut dishonest_poly = honest_poly.clone();
        dishonest_poly.set_coefficient(vec![5, 0, 0], F::from(1)); // Add x^5 term

        let mut prover = Prover::new(dishonest_poly);
        let mut verifier = Verifier::new(honest_poly, Some(42));

        // Initialize the protocol
        let prover_msg = prover
            .message(VerifierMessage::Initial)
            .expect("Failed to get initial message");
        let verifier_msg = verifier
            .message(prover_msg)
            .expect("Failed to process initial message");

        // The verifier should reject due to degree constraint violation
        let mut current_verifier_msg = verifier_msg;
        while let VerifierMessage::Challenge(challenge) = current_verifier_msg {
            match prover.message(VerifierMessage::Challenge(challenge)) {
                Ok(prover_msg) => {
                    current_verifier_msg = verifier
                        .message(prover_msg)
                        .expect("Failed to process prover message");

                    // If the verifier rejected, we're done
                    if let VerifierMessage::Reject(_) = current_verifier_msg {
                        break;
                    }
                }
                Err(_) => {
                    // If the prover errors, we'll count that as a rejection
                    current_verifier_msg = VerifierMessage::Reject(String::from("Prover Errored"));
                    break;
                }
            }
        }

        // Check that the verifier rejected
        match current_verifier_msg {
            VerifierMessage::Reject(_) => {
                // Test passes - we got a rejection as expected
            }
            _ => {
                panic!(
                    "Expected a Reject message but got {:?}",
                    current_verifier_msg
                );
            }
        }
    }
}
*/
