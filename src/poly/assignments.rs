use crate::Field;

/// Represents how a variable should be treated during polynomial operations,
/// such as taking a univariate slice or summing over the boolean hypercube.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VariableAssignment<F: Field> {
    /// Variable is fixed to a specific value
    Fixed(F),
    /// Variable is the "free" symbolic variable (not summed over)
    Free,
    /// Variable should be summed over boolean domain {0,1}
    Sum,
}

/// A validated set of variable assignments with exactly one free variable.
/// This is a common precondition for operations on `MultivariatePolynomial`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedAssignments<F: Field> {
    /// The assignments for each variable
    assignments: Vec<VariableAssignment<F>>,
    /// The index of the single free variable
    free_variable_index: usize,
}

impl<F: Field> ValidatedAssignments<F> {
    /// Creates a new validated set of assignments, ensuring exactly one variable is marked as Free.
    pub fn new(assignments: Vec<VariableAssignment<F>>) -> Result<Self, &'static str> {
        let mut free_index = None;
        for (i, assignment) in assignments.iter().enumerate() {
            if *assignment == VariableAssignment::Free {
                if free_index.is_some() {
                    return Err("Multiple free variables found in assignments");
                } else {
                    free_index = Some(i);
                }
            }
        }
        if let Some(i) = free_index {
            Ok(Self {
                assignments,
                free_variable_index: i,
            })
        } else {
            Err("No free variable found in assignments")
        }
    }

    /// Returns a reference to the underlying assignments.
    pub fn assignments(&self) -> &[VariableAssignment<F>] {
        &self.assignments
    }

    /// Returns the index of the free variable.
    pub fn free_variable_index(&self) -> usize {
        self.free_variable_index
    }

    /// Creates assignments for a specific round of the sumcheck protocol.
    ///
    /// # Arguments
    /// * `round` - Current round index (0-based)
    /// * `challenges` - Previous challenge values from completed rounds
    /// * `num_vars` - Total number of variables in the polynomial
    ///
    /// # Returns
    /// Validated assignments for this round, or an error if parameters are invalid
    pub fn new_for_round(
        round: usize,
        challenges: &[F],
        num_vars: usize,
    ) -> Result<Self, &'static str> {
        if round >= num_vars {
            return Err("Round index out of bounds");
        }

        if challenges.len() != round {
            return Err("Incorrect number of challenge values for this round");
        }

        let mut assignments = Vec::with_capacity(num_vars);

        // Set previous variables to Fixed with challenge values
        for i in 0..round {
            assignments.push(VariableAssignment::Fixed(challenges[i].clone()));
        }

        // Set current variable to Free
        assignments.push(VariableAssignment::Free);

        // Set remaining variables to Sum
        for _ in (round + 1)..num_vars {
            assignments.push(VariableAssignment::Sum);
        }

        Self::new(assignments)
    }
}

// Implement indexing for ValidatedAssignments
impl<F: Field> std::ops::Index<usize> for ValidatedAssignments<F> {
    type Output = VariableAssignment<F>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.assignments[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ff::MontgomeryFp;

    // Use a small prime for testing
    const TEST_PRIME: u32 = 13;
    // Testing only a fixed (small) field as nothing depends on
    // field arithmetic
    type TestField = MontgomeryFp<TEST_PRIME>;

    #[test]
    fn test_validated_assignments_new_success() {
        let fixed_val = TestField::from_u32(5);
        // Test with exactly one free variable
        let assignments = vec![
            VariableAssignment::Fixed(fixed_val),
            VariableAssignment::Free,
            VariableAssignment::Sum,
        ];

        let validated = ValidatedAssignments::new(assignments).unwrap();
        assert_eq!(validated.free_variable_index(), 1);
        assert_eq!(validated.assignments().len(), 3);

        assert_eq!(
            validated.assignments[0],
            VariableAssignment::Fixed(fixed_val)
        );
        assert_eq!(validated.assignments[1], VariableAssignment::Free);
        assert_eq!(validated.assignments[2], VariableAssignment::Sum);
    }

    #[test]
    fn test_validated_assignments_new_no_free() {
        // Test with no free variables
        let assignments = vec![
            VariableAssignment::Fixed(TestField::from_u32(1)),
            VariableAssignment::Fixed(TestField::from_u32(2)),
            VariableAssignment::Sum,
        ];

        let result = ValidatedAssignments::new(assignments);
        assert!(result.is_err());
    }

    #[test]
    fn test_validated_assignments_new_multiple_free() {
        // Test with multiple free variables
        let assignments = vec![
            VariableAssignment::Free,
            VariableAssignment::Fixed(TestField::from_u32(2)),
            VariableAssignment::Free,
        ];

        let result = ValidatedAssignments::new(assignments);
        assert!(result.is_err());
    }

    #[test]
    fn test_for_round_success() {
        // Test creating assignments for round 2 with 5 variables
        let mut challenges = vec![TestField::from_u32(10), TestField::from_u32(7)];
        let validated = ValidatedAssignments::new_for_round(2, &challenges, 5).unwrap();

        assert_eq!(validated.free_variable_index(), 2);
        assert_eq!(validated.assignments().len(), 5);

        let handmade_assignments = vec![
            VariableAssignment::Fixed(TestField::from_u32(10)),
            VariableAssignment::Fixed(TestField::from_u32(7)),
            VariableAssignment::Free,
            VariableAssignment::Sum,
            VariableAssignment::Sum,
        ];

        assert_eq!(
            validated,
            ValidatedAssignments::new(handmade_assignments).unwrap()
        );
    }

    #[test]
    fn test_for_round_round_out_of_bounds() {
        // Test with round >= num_vars
        let challenges = vec![TestField::from_u32(10), TestField::from_u32(7)];
        let result = ValidatedAssignments::new_for_round(5, &challenges, 5);

        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Round index out of bounds");
    }

    #[test]
    fn test_for_round_incorrect_challenges() {
        // Test with incorrect number of challenges
        let challenges = vec![
            TestField::from_u32(10),
            TestField::from_u32(7),
            TestField::from_u32(3),
        ];
        let result = ValidatedAssignments::new_for_round(2, &challenges, 5);

        assert!(result.is_err());
    }

    #[test]
    fn test_for_round_first_round() {
        // Test the first round (round 0) - should have no fixed variables
        let empty_challenges: Vec<TestField> = vec![];
        let validated = ValidatedAssignments::new_for_round(0, &empty_challenges, 3).unwrap();

        assert_eq!(validated.free_variable_index(), 0);
        assert_eq!(validated.assignments().len(), 3);
        for assignment in validated.assignments() {
            if let VariableAssignment::Fixed(_) = assignment {
                panic!("Found a fixed assignment, shouldn't be any in round 0");
            }
        }
    }

    #[test]
    fn test_for_round_last_round() {
        // Test the last round (all previous variables fixed)
        let challenges = vec![
            TestField::from_u32(4),
            TestField::from_u32(7),
            TestField::from_u32(2),
        ];
        let validated = ValidatedAssignments::new_for_round(3, &challenges, 4).unwrap();

        assert_eq!(validated.free_variable_index(), 3);
        assert_eq!(validated.assignments().len(), 4);

        for (i, assignment) in validated.assignments.iter().enumerate() {
            if i < 3 {
                assert!(matches!(assignment, VariableAssignment::Fixed(_)));
            } else {
                assert_eq!(assignment, &VariableAssignment::Free);
            }
        }
    }

    #[test]
    fn test_indexing() {
        // Test indexing capability
        let assignments = vec![
            VariableAssignment::Fixed(TestField::from_u32(5)),
            VariableAssignment::Free,
            VariableAssignment::Sum,
        ];

        let validated = ValidatedAssignments::new(assignments).unwrap();

        assert_eq!(
            validated[0],
            VariableAssignment::Fixed(TestField::from_u32(5))
        );
        assert_eq!(validated[1], VariableAssignment::Free);
        assert_eq!(validated[2], VariableAssignment::Sum);
    }
}
