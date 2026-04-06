// [[Rcpp::depends(RcppArmadillo)]]
#include "affinity_distance.h"
#include <cmath>

// Somatic Hypermutation Engine -- C++ backend
// Five mutation strategies inspired by real immunological mechanisms

// Strategy 1: UNIFORM -- classic random Gaussian noise (current bHIVE behavior)
// mutation_rate = max((1 - affinity) * decay^(iter-1), mutationMin)
// mutated = antibody + N(0, mutation_rate)
static arma::rowvec mutate_uniform(const arma::rowvec& antibody,
                                   double affinity, int iter,
                                   double decay, double mutationMin) {
  arma::uword d = antibody.n_elem;
  double rate = std::max((1.0 - affinity) * std::pow(decay, iter - 1), mutationMin);
  arma::rowvec noise(d);
  for (arma::uword i = 0; i < d; ++i) {
    noise(i) = R::rnorm(0.0, rate);
  }
  return antibody + noise;
}

// Strategy 2: AIRS -- affinity-proportional mutation (Watkins & Timmis, AIRS2)
// rate = c * exp(-affinity / temperature)
// High affinity -> fine-tuning, low affinity -> broad exploration
static arma::rowvec mutate_airs(const arma::rowvec& antibody,
                                double affinity, double c_rate,
                                double temperature) {
  arma::uword d = antibody.n_elem;
  double rate = c_rate * std::exp(-affinity / temperature);
  arma::rowvec noise(d);
  for (arma::uword i = 0; i < d; ++i) {
    noise(i) = R::rnorm(0.0, rate);
  }
  return antibody + noise;
}

// Strategy 3: HOTSPOT -- feature-importance-weighted mutation
// Features with higher gradient magnitude mutate more (analogous to AID targeting)
// gradient = data_point - antibody (direction of improvement)
static arma::rowvec mutate_hotspot(const arma::rowvec& antibody,
                                   const arma::rowvec& data_point,
                                   double affinity, double base_rate) {
  arma::uword d = antibody.n_elem;
  arma::rowvec gradient = data_point - antibody;
  arma::rowvec grad_mag = arma::abs(gradient);

  // Normalize gradient magnitudes to get per-feature weights
  double total = arma::accu(grad_mag) + 1e-15;
  arma::rowvec weights = grad_mag / total;

  // Scale mutation rate by feature importance
  double overall_rate = (1.0 - affinity) * base_rate;
  arma::rowvec noise(d);
  for (arma::uword i = 0; i < d; ++i) {
    double feature_rate = overall_rate * (1.0 + d * weights(i));
    noise(i) = R::rnorm(0.0, feature_rate);
  }
  return antibody + noise;
}

// Strategy 4: ENERGY -- mutation budget constrained by E = E_0 * (1-affinity)^2
// Inspired by Kleinstein's E_SHM ~ N_Mut^2 model
// Total mutation magnitude is bounded, distributed across features
static arma::rowvec mutate_energy(const arma::rowvec& antibody,
                                  double affinity, double E_0) {
  arma::uword d = antibody.n_elem;

  // Energy budget
  double E_budget = E_0 * std::pow(1.0 - affinity, 2);

  // Generate random direction, then scale to energy budget
  arma::rowvec direction(d);
  double norm_sq = 0.0;
  for (arma::uword i = 0; i < d; ++i) {
    direction(i) = R::rnorm(0.0, 1.0);
    norm_sq += direction(i) * direction(i);
  }
  double norm = std::sqrt(norm_sq + 1e-15);

  // Magnitude sampled from [0, sqrt(E_budget)]
  double max_magnitude = std::sqrt(E_budget);
  double magnitude = R::runif(0.0, max_magnitude);

  return antibody + (magnitude / norm) * direction;
}

// Strategy 5: ADAPTIVE -- per-feature adaptive rate with moment tracking
// Each feature maintains running mean (m1) and variance (m2) of past gradients
// Mutation rate per feature adapts like Adam optimizer
// [[Rcpp::export]]
Rcpp::List shm_mutate_cpp(const arma::mat& A,
                          const arma::mat& X,
                          const arma::vec& affinities,
                          const arma::uvec& top_k_indices,
                          const arma::uvec& data_indices,
                          const std::string& method,
                          int iter,
                          double decay,
                          double mutationMin,
                          double c_rate,
                          double temperature,
                          double E_0,
                          double base_rate,
                          double beta1,
                          double beta2,
                          double adam_epsilon,
                          arma::mat m1_state,
                          arma::mat m2_state,
                          const std::string& affinity_type,
                          double aff_alpha,
                          double aff_c,
                          double aff_p) {
  arma::uword d = A.n_cols;
  arma::mat A_out = A;

  for (arma::uword idx = 0; idx < data_indices.n_elem; ++idx) {
    arma::uword i = data_indices(idx);     // data point index
    arma::rowvec x_i = X.row(i);

    for (arma::uword ki = 0; ki < top_k_indices.n_elem; ++ki) {
      arma::uword j = top_k_indices(ki);   // antibody index
      double aff = affinities(ki);

      arma::rowvec mutated;

      if (method == "uniform") {
        mutated = mutate_uniform(A_out.row(j), aff, iter, decay, mutationMin);
      } else if (method == "airs") {
        mutated = mutate_airs(A_out.row(j), aff, c_rate, temperature);
      } else if (method == "hotspot") {
        mutated = mutate_hotspot(A_out.row(j), x_i, aff, base_rate);
      } else if (method == "energy") {
        mutated = mutate_energy(A_out.row(j), aff, E_0);
      } else if (method == "adaptive") {
        // Adam-inspired adaptive per-feature mutation
        arma::rowvec gradient = x_i - A_out.row(j);

        // Update moments
        m1_state.row(j) = beta1 * m1_state.row(j) + (1.0 - beta1) * gradient;
        m2_state.row(j) = beta2 * m2_state.row(j) + (1.0 - beta2) * (gradient % gradient);

        // Bias correction
        double bc1 = 1.0 - std::pow(beta1, iter);
        double bc2 = 1.0 - std::pow(beta2, iter);
        arma::rowvec m1_hat = m1_state.row(j) / bc1;
        arma::rowvec m2_hat = m2_state.row(j) / bc2;

        // Adaptive step size (like Adam learning rate)
        double lr = (1.0 - aff) * base_rate;
        arma::rowvec step = lr * m1_hat / (arma::sqrt(m2_hat) + adam_epsilon);

        // Add stochastic noise scaled by adaptive rate
        arma::rowvec noise(d);
        for (arma::uword dd = 0; dd < d; ++dd) {
          double feat_rate = std::abs(step(dd)) + mutationMin;
          noise(dd) = R::rnorm(0.0, feat_rate);
        }
        mutated = A_out.row(j) + step + 0.1 * noise;
      } else {
        Rcpp::stop("Unknown SHM method: " + method);
      }

      // Evaluate mutated antibody using scalar function (no matrix allocation)
      const AffinityType at = parse_affinity_type(affinity_type);
      const double f_mutated = affinity_scalar_cpp(
        x_i, mutated, at, aff_alpha, aff_c, aff_p
      );

      if (std::isfinite(f_mutated) && f_mutated > aff) {
        A_out.row(j) = mutated;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("A") = A_out,
    Rcpp::Named("m1_state") = m1_state,
    Rcpp::Named("m2_state") = m2_state
  );
}
