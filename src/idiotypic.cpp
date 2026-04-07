// [[Rcpp::depends(RcppArmadillo)]]
#include "affinity_distance.h"
#include <cmath>

// Idiotypic Network Dynamics -- C++ backend
// Implements Varela & Coutinho's second-generation immune network model
// with bell-shaped (double-threshold) activation function

// Bell-shaped activation: too little = death, moderate = activation, too much = suppression
// Returns: -1 (suppressed), 0 (neutral), +1 (activated)
static double bell_activation(double stimulation, double theta_low, double theta_high) {
  if (stimulation < theta_low) return -1.0;  // insufficient stimulation -> death
  if (stimulation > theta_high) return -1.0; // over-stimulation -> suppression
  // In the activation window: scale linearly from 0 to 1 at peak, back to 0
  const double mid = (theta_low + theta_high) / 2.0;
  const double half_width = (theta_high - theta_low) / 2.0;
  if (half_width <= 0.0) return 0.0;  // degenerate case: theta_low == theta_high
  return 1.0 - std::abs(stimulation - mid) / half_width;
}


// Simulate idiotypic network dynamics
// Population dynamics: dA_i/dt = source - decay*A_i + A_i * sum(activation) - A_i * sum(suppression)
//
// Returns:
//   - population: final population levels per antibody
//   - keep: logical vector of which antibodies survive
//   - Ab_Ab_affinity: the antibody-antibody affinity matrix
//
// [[Rcpp::export]]
Rcpp::List idiotypic_dynamics_cpp(const arma::mat& A,
                                  const std::string& affinity_type,
                                  double aff_alpha,
                                  double aff_c,
                                  double aff_p,
                                  double theta_low,
                                  double theta_high,
                                  double source_rate,
                                  double decay_rate,
                                  double dt,
                                  int timeSteps,
                                  double survival_threshold) {
  arma::uword m = A.n_rows;

  // Compute Ab-Ab affinity matrix
  AffinityType at = parse_affinity_type(affinity_type);
  arma::mat Ab_Ab = affinity_matrix_cpp(A, A, at, aff_alpha, aff_c, aff_p);

  // Zero out diagonal (no self-interaction)
  Ab_Ab.diag().zeros();

  // Initialize population levels (all start at 1.0)
  arma::vec population = arma::ones<arma::vec>(m);

  // Simulate dynamics
  for (int t = 0; t < timeSteps; ++t) {
    arma::vec dp = arma::zeros<arma::vec>(m);

    for (arma::uword i = 0; i < m; ++i) {
      double activation_sum = 0.0;
      double suppression_sum = 0.0;

      for (arma::uword j = 0; j < m; ++j) {
        if (i == j) continue;

        // Total stimulation from j to i = affinity * population_j
        double stimulation = Ab_Ab(i, j) * population(j);
        double response = bell_activation(stimulation, theta_low, theta_high);

        if (response > 0) {
          activation_sum += response;
        } else if (response < 0) {
          suppression_sum += std::abs(response);
        }
      }

      // Population dynamics ODE
      dp(i) = source_rate
             - decay_rate * population(i)
             + population(i) * activation_sum
             - population(i) * suppression_sum;
    }

    // Euler step
    population += dt * dp;

    // Clamp to non-negative
    population.clamp(0.0, arma::datum::inf);
  }

  // Determine survivors
  Rcpp::LogicalVector keep(m);
  for (arma::uword i = 0; i < m; ++i) {
    keep[i] = population(i) > survival_threshold;
  }

  return Rcpp::List::create(
    Rcpp::Named("population") = population,
    Rcpp::Named("keep") = keep,
    Rcpp::Named("Ab_Ab_affinity") = Ab_Ab
  );
}
