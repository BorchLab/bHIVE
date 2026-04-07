// [[Rcpp::depends(RcppArmadillo)]]
#include "affinity_distance.h"

// Network suppression: remove antibodies within epsilon distance of each other
// Returns logical vector indicating which antibodies to keep
// [[Rcpp::export]]
Rcpp::LogicalVector network_suppression_cpp(const arma::mat& A,
                                            const std::string& dist_type,
                                            double epsilon,
                                            double p,
                                            const arma::mat& Sigma_inv) {
  arma::uword m = A.n_rows;
  DistanceType dtype = parse_distance_type(dist_type);

  std::vector<bool> keep(m, true);

  // Compute pairwise distances only as needed (upper triangle)
  for (arma::uword u = 0; u < m; ++u) {
    if (!keep[u]) continue;
    for (arma::uword v = u + 1; v < m; ++v) {
      if (!keep[v]) continue;

      double dist_uv;
      arma::rowvec diff = A.row(u) - A.row(v);

      switch (dtype) {
        case DistanceType::EUCLIDEAN:
          dist_uv = std::sqrt(arma::accu(diff % diff));
          break;
        case DistanceType::MANHATTAN:
          dist_uv = arma::accu(arma::abs(diff));
          break;
        case DistanceType::MINKOWSKI:
          dist_uv = std::pow(arma::accu(arma::pow(arma::abs(diff), p)), 1.0 / p);
          break;
        case DistanceType::COSINE: {
          double dot = arma::accu(A.row(u) % A.row(v));
          double nu = std::sqrt(arma::accu(A.row(u) % A.row(u)));
          double nv = std::sqrt(arma::accu(A.row(v) % A.row(v)));
          double denom = nu * nv;
          dist_uv = (denom < 1e-15) ? 1.0 : (1.0 - dot / denom);
          break;
        }
        case DistanceType::MAHALANOBIS: {
          arma::rowvec tmp = diff * Sigma_inv;
          dist_uv = std::sqrt(arma::accu(tmp % diff));
          break;
        }
        case DistanceType::HAMMING: {
          arma::irowvec u_int = arma::conv_to<arma::irowvec>::from(A.row(u));
          arma::irowvec v_int = arma::conv_to<arma::irowvec>::from(A.row(v));
          dist_uv = static_cast<double>(arma::accu(u_int != v_int));
          break;
        }
        default:
          Rcpp::stop("Unknown distance type in suppression");
      }

      if (dist_uv < epsilon) {
        keep[v] = false;
      }
    }
  }

  Rcpp::LogicalVector result(m);
  for (arma::uword i = 0; i < m; ++i) {
    result[i] = keep[i];
  }
  return result;
}
