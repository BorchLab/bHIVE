// [[Rcpp::depends(RcppArmadillo)]]
#include "affinity_distance.h"

// ============================================================
// AFFINITY FUNCTIONS (bulk matrix computation)
// ============================================================

// Gaussian/RBF: exp(-alpha * ||x-y||^2)
// Uses BLAS trick: ||x-y||^2 = ||x||^2 + ||y||^2 - 2*x.y
static arma::mat affinity_gaussian(const arma::mat& X, const arma::mat& A,
                                   double alpha) {
  arma::uword n = X.n_rows;
  arma::uword m = A.n_rows;

  // Squared norms
  arma::vec x_sq = arma::sum(X % X, 1);  // n x 1
  arma::vec a_sq = arma::sum(A % A, 1);  // m x 1

  // Squared distance matrix: n x m
  // D_sq(i,j) = ||x_i||^2 + ||a_j||^2 - 2 * x_i . a_j
  arma::mat D_sq = arma::repmat(x_sq, 1, m) +
                   arma::repmat(a_sq.t(), n, 1) -
                   2.0 * X * A.t();

  // Clamp to avoid numerical issues with negative values
  D_sq.clamp(0.0, arma::datum::inf);

  return arma::exp(-alpha * D_sq);
}

// Laplace: exp(-alpha * ||x-y||_1)
static arma::mat affinity_laplace(const arma::mat& X, const arma::mat& A,
                                  double alpha) {
  arma::uword n = X.n_rows;
  arma::uword m = A.n_rows;
  arma::mat result(n, m);

  for (arma::uword j = 0; j < m; ++j) {
    arma::rowvec a_j = A.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      double l1 = arma::accu(arma::abs(X.row(i) - a_j));
      result(i, j) = std::exp(-alpha * l1);
    }
  }
  return result;
}

// Polynomial: (x.y + c)^p
static arma::mat affinity_polynomial(const arma::mat& X, const arma::mat& A,
                                     double c, double p) {
  arma::mat dot = X * A.t();  // n x m via BLAS
  return arma::pow(dot + c, p);
}

// Cosine similarity: (x.y) / (||x|| * ||y||)
static arma::mat affinity_cosine(const arma::mat& X, const arma::mat& A) {
  arma::vec x_norm = arma::sqrt(arma::sum(X % X, 1));  // n x 1
  arma::vec a_norm = arma::sqrt(arma::sum(A % A, 1));  // m x 1

  arma::mat dot = X * A.t();  // n x m via BLAS

  // Denominator matrix
  arma::mat denom = x_norm * a_norm.t();

  // Avoid division by zero (1e-12 is safer than 1e-15 across platforms)
  arma::mat result = dot / (denom + 1e-12);
  result.clamp(-1.0, 1.0);

  return result;
}

// Hamming similarity: proportion of matching elements
static arma::mat affinity_hamming(const arma::mat& X, const arma::mat& A) {
  const arma::uword n = X.n_rows;
  const arma::uword m = A.n_rows;
  const arma::uword d = X.n_cols;
  const double inv_d = 1.0 / static_cast<double>(d);
  arma::mat result(n, m);

  // Pre-convert both matrices to integer once (avoids per-pair conversion)
  const arma::imat X_int = arma::conv_to<arma::imat>::from(X);
  const arma::imat A_int = arma::conv_to<arma::imat>::from(A);

  for (arma::uword j = 0; j < m; ++j) {
    const arma::irowvec a_int = A_int.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      result(i, j) = static_cast<double>(arma::accu(X_int.row(i) == a_int)) * inv_d;
    }
  }
  return result;
}

// Dispatcher (bulk matrix)
arma::mat affinity_matrix_cpp(const arma::mat& X, const arma::mat& A,
                              AffinityType type, double alpha, double c,
                              double p) {
  switch (type) {
    case AffinityType::GAUSSIAN:
      return affinity_gaussian(X, A, alpha);
    case AffinityType::LAPLACE:
      return affinity_laplace(X, A, alpha);
    case AffinityType::POLYNOMIAL:
      return affinity_polynomial(X, A, c, p);
    case AffinityType::COSINE:
      return affinity_cosine(X, A);
    case AffinityType::HAMMING:
      return affinity_hamming(X, A);
    default:
      Rcpp::stop("Unknown affinity type");
      return arma::mat();  // unreachable
  }
}


// ============================================================
// SCALAR AFFINITY (single point-to-point, no matrix allocation)
// Critical for hot-path clone/mutate evaluations
// ============================================================

double affinity_scalar_cpp(const arma::rowvec& x, const arma::rowvec& a,
                           AffinityType type, double alpha, double c,
                           double p) {
  switch (type) {
    case AffinityType::GAUSSIAN: {
      const double dist2 = arma::accu(arma::square(x - a));
      return std::exp(-alpha * dist2);
    }
    case AffinityType::LAPLACE: {
      const double l1 = arma::accu(arma::abs(x - a));
      return std::exp(-alpha * l1);
    }
    case AffinityType::POLYNOMIAL: {
      return std::pow(arma::dot(x, a) + c, p);
    }
    case AffinityType::COSINE: {
      const double denom = arma::norm(x) * arma::norm(a);
      if (denom < 1e-12) return 0.0;
      double cs = arma::dot(x, a) / denom;
      return std::max(-1.0, std::min(1.0, cs));
    }
    case AffinityType::HAMMING: {
      const arma::uword d = x.n_elem;
      const arma::irowvec xi = arma::conv_to<arma::irowvec>::from(x);
      const arma::irowvec ai = arma::conv_to<arma::irowvec>::from(a);
      return static_cast<double>(arma::accu(xi == ai)) / static_cast<double>(d);
    }
    default:
      return 0.0;
  }
}


// ============================================================
// DISTANCE FUNCTIONS (bulk matrix computation)
// ============================================================

// Euclidean distance: sqrt(||x-y||^2)
static arma::mat distance_euclidean(const arma::mat& X, const arma::mat& A) {
  arma::vec x_sq = arma::sum(X % X, 1);
  arma::vec a_sq = arma::sum(A % A, 1);
  arma::uword n = X.n_rows;
  arma::uword m = A.n_rows;

  arma::mat D_sq = arma::repmat(x_sq, 1, m) +
                   arma::repmat(a_sq.t(), n, 1) -
                   2.0 * X * A.t();
  D_sq.clamp(0.0, arma::datum::inf);
  return arma::sqrt(D_sq);
}

// Manhattan distance: ||x-y||_1
static arma::mat distance_manhattan(const arma::mat& X, const arma::mat& A) {
  arma::uword n = X.n_rows;
  arma::uword m = A.n_rows;
  arma::mat result(n, m);

  for (arma::uword j = 0; j < m; ++j) {
    arma::rowvec a_j = A.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      result(i, j) = arma::accu(arma::abs(X.row(i) - a_j));
    }
  }
  return result;
}

// Minkowski distance: (sum |x-y|^p)^(1/p)
static arma::mat distance_minkowski(const arma::mat& X, const arma::mat& A,
                                    double p) {
  // For p=2, delegate to Euclidean (BLAS-optimized)
  if (std::abs(p - 2.0) < 1e-10) return distance_euclidean(X, A);

  const arma::uword n = X.n_rows;
  const arma::uword m = A.n_rows;
  arma::mat result(n, m);
  const double inv_p = 1.0 / p;

  for (arma::uword j = 0; j < m; ++j) {
    const arma::rowvec a_j = A.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      // Clamp to >= 0 before pow to prevent NaN from numerical error
      double val = arma::accu(arma::pow(arma::abs(X.row(i) - a_j), p));
      result(i, j) = std::pow(std::max(0.0, val), inv_p);
    }
  }
  return result;
}

// Cosine distance: 1 - cosine_similarity
static arma::mat distance_cosine(const arma::mat& X, const arma::mat& A) {
  return 1.0 - affinity_cosine(X, A);
}

// Mahalanobis distance: sqrt((x-y)^T Sigma^{-1} (x-y))
// Sigma_inv is pre-computed inverse covariance matrix
static arma::mat distance_mahalanobis(const arma::mat& X, const arma::mat& A,
                                      const arma::mat& Sigma_inv) {
  arma::uword n = X.n_rows;
  arma::uword m = A.n_rows;
  arma::mat result(n, m);

  for (arma::uword j = 0; j < m; ++j) {
    arma::rowvec a_j = A.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      arma::rowvec diff = X.row(i) - a_j;
      arma::rowvec tmp = diff * Sigma_inv;
      result(i, j) = std::sqrt(arma::accu(tmp % diff));
    }
  }
  return result;
}

// Hamming distance: count of mismatches
static arma::mat distance_hamming(const arma::mat& X, const arma::mat& A) {
  const arma::uword n = X.n_rows;
  const arma::uword m = A.n_rows;
  arma::mat result(n, m);

  // Pre-convert both matrices to integer once
  const arma::imat X_int = arma::conv_to<arma::imat>::from(X);
  const arma::imat A_int = arma::conv_to<arma::imat>::from(A);

  for (arma::uword j = 0; j < m; ++j) {
    const arma::irowvec a_int = A_int.row(j);
    for (arma::uword i = 0; i < n; ++i) {
      result(i, j) = static_cast<double>(arma::accu(X_int.row(i) != a_int));
    }
  }
  return result;
}

// Dispatcher
arma::mat distance_matrix_cpp(const arma::mat& X, const arma::mat& A,
                              DistanceType type, double p,
                              const arma::mat& Sigma_inv) {
  switch (type) {
    case DistanceType::EUCLIDEAN:
      return distance_euclidean(X, A);
    case DistanceType::MANHATTAN:
      return distance_manhattan(X, A);
    case DistanceType::MINKOWSKI:
      return distance_minkowski(X, A, p);
    case DistanceType::COSINE:
      return distance_cosine(X, A);
    case DistanceType::MAHALANOBIS:
      return distance_mahalanobis(X, A, Sigma_inv);
    case DistanceType::HAMMING:
      return distance_hamming(X, A);
    default:
      Rcpp::stop("Unknown distance type");
      return arma::mat();  // unreachable
  }
}


// ============================================================
// PAIRWISE DISTANCE (antibody-antibody, m x m)
// ============================================================

arma::mat pairwise_distance_cpp(const arma::mat& A, DistanceType type,
                                double p, const arma::mat& Sigma_inv) {
  return distance_matrix_cpp(A, A, type, p, Sigma_inv);
}


// ============================================================
// R-EXPOSED FUNCTIONS
// ============================================================

// [[Rcpp::export]]
arma::mat compute_affinity_matrix(const arma::mat& X, const arma::mat& A,
                                  const std::string& affinity_type,
                                  double alpha = 1.0, double c = 1.0,
                                  double p = 2.0) {
  AffinityType type = parse_affinity_type(affinity_type);
  return affinity_matrix_cpp(X, A, type, alpha, c, p);
}

// [[Rcpp::export]]
arma::mat compute_distance_matrix(const arma::mat& X, const arma::mat& A,
                                  const std::string& dist_type,
                                  double p,
                                  const arma::mat& Sigma_inv) {
  DistanceType type = parse_distance_type(dist_type);
  return distance_matrix_cpp(X, A, type, p, Sigma_inv);
}

// [[Rcpp::export]]
arma::mat compute_pairwise_distance(const arma::mat& A,
                                    const std::string& dist_type,
                                    double p,
                                    const arma::mat& Sigma_inv) {
  DistanceType type = parse_distance_type(dist_type);
  return pairwise_distance_cpp(A, type, p, Sigma_inv);
}
