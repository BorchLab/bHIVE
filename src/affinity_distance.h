#ifndef BHIVE_AFFINITY_DISTANCE_H
#define BHIVE_AFFINITY_DISTANCE_H

#include "bHIVE_types.h"

// ============================================================
// Bulk affinity matrix: returns n x m matrix of affinities
// X is n x d (data), A is m x d (antibodies)
// ============================================================

arma::mat affinity_matrix_cpp(const arma::mat& X, const arma::mat& A,
                              AffinityType type, double alpha, double c,
                              double p);

// ============================================================
// Scalar affinity: returns single affinity value between two vectors
// Critical for hot-path single-point evaluations in clone/mutate
// ============================================================

double affinity_scalar_cpp(const arma::rowvec& x, const arma::rowvec& a,
                           AffinityType type, double alpha, double c,
                           double p);

// ============================================================
// Bulk distance matrix: returns n x m matrix of distances
// X is n x d (data), A is m x d (antibodies)
// ============================================================

arma::mat distance_matrix_cpp(const arma::mat& X, const arma::mat& A,
                              DistanceType type, double p,
                              const arma::mat& Sigma_inv);

// ============================================================
// Pairwise distance matrix: returns m x m matrix (antibody-antibody)
// Used for network suppression and idiotypic interactions
// ============================================================

arma::mat pairwise_distance_cpp(const arma::mat& A, DistanceType type,
                                double p, const arma::mat& Sigma_inv);

#endif // BHIVE_AFFINITY_DISTANCE_H
