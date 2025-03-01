//! This example shows a step-by-step execution of the Sumcheck protocol
//! for a simple multilinear polynomial with 3 variables.

use fabric_challenge::{
    Field, GeneralMultivariatePolynomial, MontgomeryFp, MultivariatePolynomial, Prover, Verifier,
    VerifierMessage,
};
use std::collections::HashMap;
use std::time::Instant;

fn main() {
    println!("=================================================");
    println!("       SUMCHECK PROTOCOL DEMONSTRATION");
    println!("        (General Multivariate Backend)");
    println!("=================================================");

    // Define our field - using a small prime (13) for readability
    const PRIME: u32 = 13;
    type F = MontgomeryFp<PRIME>;

    let start_time = Instant::now();

    println!("\n1. SETUP PHASE");
    println!("---------------");
    println!("Field: GF({})", PRIME);

    // Create our example multilinear polynomial:
    // p(x_1, x_2, x_3) = 2x_1x_2x_3 + 3x_1x_2 + x_2x_3 + 4x_1 + 2
    let mut coeffs = HashMap::new();
    coeffs.insert(vec![1, 1, 1], F::from(2)); // 2x_1x_2x_3
    coeffs.insert(vec![1, 1, 0], F::from(3)); // 3x_1x_2
    coeffs.insert(vec![0, 1, 1], F::from(1)); // x_2x_3
    coeffs.insert(vec![1, 0, 0], F::from(4)); // 4x_1
    coeffs.insert(vec![0, 0, 0], F::from(2)); // 2 (constant term)

    let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

    println!("\nPolynomial: p(x_1, x_2, x_3) = 2x_1x_2x_3 + 3x_1x_2 + x_2x_3 + 4x_1 + 2");

    let setup_time = start_time.elapsed();
    println!("Setup completed in: {:?}", setup_time);

    // Calculate the sum over the boolean hypercube
    let true_sum = poly.sum_over_boolean_hypercube();

    let sum_calc_time = start_time.elapsed() - setup_time;
    println!("Sum calculation completed in: {:?}", sum_calc_time);

    // Calculate and display sum over all boolean inputs
    println!("\nCalculating sum over all boolean inputs (2^3 = 8 evaluations):");

    let mut evaluations = Vec::new();

    // Calculate all evaluations
    for x1 in 0..2 {
        for x2 in 0..2 {
            for x3 in 0..2 {
                let point = vec![F::from(x1), F::from(x2), F::from(x3)];
                let value = poly.evaluate(&point);

                // Calculate each term contribution
                let term1 = 2 * x1 * x2 * x3; // 2x1x2x3
                let term2 = 3 * x1 * x2; // 3x1x2
                #[allow(clippy::identity_op)]
                let term3 = 1 * x2 * x3; // x2x3
                let term4 = 4 * x1; // 4x1
                let term5 = 2; // constant term

                // Sum should match the polynomial evaluation
                let sum = term1 + term2 + term3 + term4 + term5;

                println!(
                    "  p({}, {}, {}) = {} + {} + {} + {} + {} = {}",
                    x1, x2, x3, term1, term2, term3, term4, term5, sum
                );

                evaluations.push(value);
            }
        }
    }

    // Show the simple sum verification
    println!(
        "\n{} + {} + {} + {} + {} + {} + {} + {} ≡ {} (mod {})",
        evaluations[0],
        evaluations[1],
        evaluations[2],
        evaluations[3],
        evaluations[4],
        evaluations[5],
        evaluations[6],
        evaluations[7],
        true_sum,
        PRIME
    );

    println!("\nTrue sum: H = Sum p(x_1, x_2, x_3) = {}", true_sum);

    // Initialize the Sumcheck protocol
    println!("\n2. PROTOCOL EXECUTION");
    println!("---------------------");

    println!("\nInitializing Prover with the polynomial");
    let mut prover = Prover::new(poly.clone());

    println!(
        "Initializing Verifier with the polynomial and claimed sum: {}",
        true_sum
    );
    let mut verifier = Verifier::new(poly.clone(), Some(42)); // Fixed seed for reproducibility

    println!("\n[ROUND 1] - Variable x_3");
    println!("------------------------");

    println!("Prover's task: Compute g_3(X_3) = Sum_{{x_1,x_2 in {{0,1}}}} p(x_1, x_2, X_3)");

    let round1_time = Instant::now();

    // Round 1: Initial message exchange
    let prover_msg1 = prover.message(VerifierMessage::Initial).unwrap();

    // Display information from the prover's message
    if let Some(claimed_sum) = prover_msg1.get_claimed_sum() {
        println!("\nProver claims the sum is: {}", claimed_sum);
    }
    println!(
        "Prover sends univariate polynomial g_3(X_3) = {}",
        prover_msg1.get_polynomial()
    );

    // For demonstration, we'll calculate the values g_3(0) and g_3(1)
    let msg1 = prover_msg1.get_polynomial().clone();
    let g3_at_0 = msg1.evaluate(&F::zero());
    let g3_at_1 = msg1.evaluate(&F::one());

    println!("\nVerifier checks:");
    println!("  g_3(0) = {}", g3_at_0);
    println!("  g_3(1) = {}", g3_at_1);
    println!(
        "  g_3(0) + g_3(1) = {} + {} = {}",
        g3_at_0,
        g3_at_1,
        g3_at_0 + g3_at_1
    );
    println!("  Claimed sum H = {}", true_sum);
    println!(
        "  Check: g_3(0) + g_3(1) = H? {}",
        g3_at_0 + g3_at_1 == true_sum
    );

    // Verifier processes the prover's message
    let verifier_msg1 = verifier.message(prover_msg1).unwrap();

    // Get the challenge from the verifier's response
    if let Some(r3) = verifier_msg1.get_challenge() {
        println!("\nVerifier sends random challenge r_3 = {}", r3);

        println!(
            "\nBoth parties now focus on p(x_1, x_2, r_3) = p(x_1, x_2, {})",
            r3
        );

        println!("\n[ROUND 2] - Variable x_2");
        println!("------------------------");

        println!("Prover's task: Compute g_2(X_2) = Sum_{{x_1 in {{0,1}}}} p(x_1, X_2, r_3)");

        println!("Round 1 completed in: {:?}", round1_time.elapsed());

        // Round 2: Prover processes challenge and sends next message
        let round2_time = Instant::now();
        #[allow(clippy::clone_on_copy)]
        let r3_owned = r3.clone(); // Clone for ownership
        let prover_msg2 = prover
            .message(VerifierMessage::Challenge(r3_owned))
            .unwrap();

        // Get the polynomial from the prover's message
        let msg2 = prover_msg2.get_polynomial().clone();
        println!("\nProver sends univariate polynomial g_2(X_2) = {}", msg2);

        // For demonstration, calculate g_2(0), g_2(1), and g_3(r_3)
        let g2_at_0 = msg2.evaluate(&F::zero());
        let g2_at_1 = msg2.evaluate(&F::one());
        let g3_at_r3 = msg1.evaluate(r3);

        println!("\nVerifier checks:");
        println!("  g_2(0) = {}", g2_at_0);
        println!("  g_2(1) = {}", g2_at_1);
        println!(
            "  g_2(0) + g_2(1) = {} + {} = {}",
            g2_at_0,
            g2_at_1,
            g2_at_0 + g2_at_1
        );
        println!("  g_3(r_3) = g_3({}) = {}", r3, g3_at_r3);
        println!(
            "  Check: g_2(0) + g_2(1) = g_3(r_3)? {}",
            g2_at_0 + g2_at_1 == g3_at_r3
        );

        // Verifier processes the prover's message
        let verifier_msg2 = verifier.message(prover_msg2).unwrap();

        // Get the challenge from the verifier's response
        if let Some(r2) = verifier_msg2.get_challenge() {
            println!("\nVerifier sends random challenge r_2 = {}", r2);

            println!(
                "\nBoth parties now focus on p(x_1, r_2, r_3) = p(x_1, {}, {})",
                r2, r3
            );
            println!("Round 2 completed in: {:?}", round2_time.elapsed());

            println!("\n[ROUND 3] - Variable x_1 (Final Round)");
            println!("-------------------------------------");

            println!("Prover's task: Compute g_1(X_1) = p(X_1, r_2, r_3)");

            let round3_time = Instant::now();
            // Round 3: Prover processes challenge and sends final message
            #[allow(clippy::clone_on_copy)]
            let r2_owned = r2.clone(); // Clone for ownership
            let prover_msg3 = prover
                .message(VerifierMessage::Challenge(r2_owned))
                .unwrap();

            // Get the polynomial from the prover's message
            let msg3 = prover_msg3.get_polynomial();
            println!("\nProver sends univariate polynomial g_1(X_1) = {}", msg3);

            // For demonstration, calculate g_1(0), g_1(1), and g_2(r_2)
            let g1_at_0 = msg3.evaluate(&F::zero());
            let g1_at_1 = msg3.evaluate(&F::one());
            let g2_at_r2 = msg2.evaluate(r2);

            println!("\nVerifier checks:");
            println!("  g_1(0) = {}", g1_at_0);
            println!("  g_1(1) = {}", g1_at_1);
            println!(
                "  g_1(0) + g_1(1) = {} + {} = {}",
                g1_at_0,
                g1_at_1,
                g1_at_0 + g1_at_1
            );
            println!("  g_2(r_2) = g_2({}) = {}", r2, g2_at_r2);
            println!(
                "  Check: g_1(0) + g_1(1) = g_2(r_2)? {}",
                g1_at_0 + g1_at_1 == g2_at_r2
            );

            // Verifier processes the prover's message
            let verifier_final = verifier.message(prover_msg3).unwrap();

            println!("Round 3 completed in: {:?}", round3_time.elapsed());

            println!("\n3. FINAL VERIFICATION");
            println!("---------------------");

            // Display the result
            if verifier_final.is_accept() {
                println!("Verifier ACCEPTS the proof!");
                println!("\nThe Sumcheck protocol succeeds because:");
                println!("1. The sum H was correctly claimed by the prover");
                println!("2. All consistency checks passed at each round");
                println!("3. The final evaluation matched the direct computation");
            } else if let Some(reason) = verifier_final.rejection_reason() {
                println!("Verifier REJECTS the proof!");
                println!("Reason: {}", reason);
            } else {
                println!("Unexpected protocol outcome");
            }
            println!("Total execution time: {:?}", start_time.elapsed());
        }
    }

    println!("\nThis demonstrates that the prover does indeed know a polynomial");
    println!("whose sum over the Boolean hypercube equals the claimed value H,");
    println!("without revealing the full polynomial to the verifier.");
}
