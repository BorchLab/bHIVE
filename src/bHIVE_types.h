#ifndef BHIVE_TYPES_H
#define BHIVE_TYPES_H

#include <RcppArmadillo.h>

enum class AffinityType {
  GAUSSIAN = 0,
  LAPLACE = 1,
  POLYNOMIAL = 2,
  COSINE = 3,
  HAMMING = 4
};

enum class DistanceType {
  EUCLIDEAN = 0,
  MANHATTAN = 1,
  MINKOWSKI = 2,
  COSINE = 3,
  MAHALANOBIS = 4,
  HAMMING = 5
};

enum class TaskType {
  CLUSTERING = 0,
  CLASSIFICATION = 1,
  REGRESSION = 2
};

// Convert R string to enum
inline AffinityType parse_affinity_type(const std::string& s) {
  if (s == "gaussian")   return AffinityType::GAUSSIAN;
  if (s == "laplace")    return AffinityType::LAPLACE;
  if (s == "polynomial") return AffinityType::POLYNOMIAL;
  if (s == "cosine")     return AffinityType::COSINE;
  if (s == "hamming")    return AffinityType::HAMMING;
  Rcpp::stop("Invalid affinity type: " + s);
  return AffinityType::GAUSSIAN; // unreachable
}

inline DistanceType parse_distance_type(const std::string& s) {
  if (s == "euclidean")   return DistanceType::EUCLIDEAN;
  if (s == "manhattan")   return DistanceType::MANHATTAN;
  if (s == "minkowski")   return DistanceType::MINKOWSKI;
  if (s == "cosine")      return DistanceType::COSINE;
  if (s == "mahalanobis") return DistanceType::MAHALANOBIS;
  if (s == "hamming")     return DistanceType::HAMMING;
  Rcpp::stop("Invalid distance type: " + s);
  return DistanceType::EUCLIDEAN; // unreachable
}

inline TaskType parse_task_type(const std::string& s) {
  if (s == "clustering")     return TaskType::CLUSTERING;
  if (s == "classification") return TaskType::CLASSIFICATION;
  if (s == "regression")     return TaskType::REGRESSION;
  Rcpp::stop("Invalid task type: " + s);
  return TaskType::CLUSTERING; // unreachable
}

#endif // BHIVE_TYPES_H
