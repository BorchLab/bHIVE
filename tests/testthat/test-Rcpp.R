# Tests for C++ backend functions (via Rcpp)

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Affinity Matrix ====

test_that("compute_affinity_matrix returns correct dimensions", {
  A <- X[1:10, ]
  aff <- compute_affinity_matrix(X, A, "gaussian", alpha = 1, c = 1, p = 2)
  expect_equal(dim(aff), c(150, 10))
})

test_that("compute_affinity_matrix gaussian values in [0, 1]", {
  A <- X[1:5, ]
  aff <- compute_affinity_matrix(X, A, "gaussian", alpha = 1, c = 1, p = 2)
  expect_true(all(aff >= 0))
  expect_true(all(aff <= 1))
})

test_that("compute_affinity_matrix self-affinity is 1 for gaussian", {
  A <- X[1:5, ]
  aff <- compute_affinity_matrix(A, A, "gaussian", alpha = 1, c = 1, p = 2)
  diag_vals <- diag(aff)
  expect_equal(diag_vals, rep(1, 5), tolerance = 1e-10)
})

test_that("compute_affinity_matrix works for all types", {
  A <- X[1:5, ]
  for (aff_type in c("gaussian", "laplace", "polynomial", "cosine")) {
    result <- compute_affinity_matrix(X, A, aff_type, alpha = 1, c = 1, p = 2)
    expect_equal(dim(result), c(150, 5))
    expect_true(all(is.finite(result)))
  }
})

# ==== Distance Matrix ====

test_that("compute_distance_matrix returns correct dimensions", {
  A <- X[1:10, ]
  dist <- compute_distance_matrix(X, A, "euclidean", p = 2.0,
                                   Sigma_inv = matrix(0, 0, 0))
  expect_equal(dim(dist), c(150, 10))
})

test_that("compute_distance_matrix self-distance is 0", {
  A <- X[1:5, ]
  dist <- compute_distance_matrix(A, A, "euclidean", p = 2.0,
                                   Sigma_inv = matrix(0, 0, 0))
  diag_vals <- diag(dist)
  expect_equal(diag_vals, rep(0, 5), tolerance = 1e-10)
})

test_that("compute_distance_matrix works for different metrics", {
  A <- X[1:5, ]
  Sigma_inv <- matrix(0, 0, 0)
  for (dist_type in c("euclidean", "manhattan", "minkowski")) {
    result <- compute_distance_matrix(X, A, dist_type, p = 2.0,
                                       Sigma_inv = Sigma_inv)
    expect_equal(dim(result), c(150, 5))
    expect_true(all(result >= 0))
    expect_true(all(is.finite(result)))
  }
})

# ==== Pairwise Distance ====

test_that("compute_pairwise_distance returns m x m matrix", {
  A <- X[1:10, ]
  pw <- compute_pairwise_distance(A, "euclidean", p = 2.0,
                                   Sigma_inv = matrix(0, 0, 0))
  expect_equal(dim(pw), c(10, 10))
  expect_true(all(diag(pw) < 0.1))  # self-distance should be small
  # Symmetric
  expect_equal(pw, t(pw), tolerance = 1e-10)
})

# ==== Network Suppression ====

test_that("network_suppression_cpp keeps distinct antibodies", {
  A <- X[1:10, ]
  keep <- network_suppression_cpp(A, "euclidean", epsilon = 0.001,
                                   p = 2.0, Sigma_inv = matrix(0, 0, 0))
  expect_true(is.logical(keep))
  expect_equal(length(keep), 10)
  expect_true(all(keep))  # very small epsilon => all kept
})

test_that("network_suppression_cpp removes duplicates", {
  A <- rbind(X[1, ], X[1, ] + 0.001, X[2, ], X[2, ] + 0.001)
  keep <- network_suppression_cpp(A, "euclidean", epsilon = 0.1,
                                   p = 2.0, Sigma_inv = matrix(0, 0, 0))
  expect_true(sum(keep) < 4)  # some duplicates removed
})

# ==== K-means++ Initialization ====

test_that("init_kmeanspp_cpp returns correct dimensions", {
  centers <- init_kmeanspp_cpp(X, 10L)
  expect_equal(dim(centers), c(10, 4))
  expect_true(all(is.finite(centers)))
})

test_that("init_kmeanspp_cpp returns unique centers", {
  centers <- init_kmeanspp_cpp(X, 5L)
  expect_equal(nrow(unique(centers)), 5)
})

# ==== Clonal Selection Iteration ====

test_that("clonal_selection_iteration_cpp runs for clustering", {
  A <- X[sample(150, 10), ]
  y_num <- rep(0, 150)
  result <- clonal_selection_iteration_cpp(
    A, X, y_num, 0L, 3L, 5.0, Inf, 1.0, 0.01, 1L,
    "gaussian", 1.0, 1.0, 2.0, 0L
  )
  expect_true(is.matrix(result$A))
  expect_equal(ncol(result$A), 4)
  expect_true(nrow(result$A) >= 10)  # clones added
})

test_that("clonal_selection_iteration_cpp runs for classification", {
  A <- X[sample(150, 10), ]
  y_num <- as.numeric(iris$Species) - 1
  result <- clonal_selection_iteration_cpp(
    A, X, y_num, 1L, 3L, 5.0, Inf, 1.0, 0.01, 1L,
    "gaussian", 1.0, 1.0, 2.0, 3L
  )
  expect_true(is.matrix(result$A))
  expect_true(!is.null(result$class_counts))
})

# ==== Final Assignment ====

test_that("final_assignment_cpp runs for clustering", {
  A <- X[sample(150, 10), ]
  fa <- final_assignment_cpp(X, A, "gaussian", "euclidean", 0L,
                              1.0, 1.0, 2.0, matrix(0, 0, 0))
  expect_true(!is.null(fa$assignments))
  expect_equal(length(fa$assignments), 150)
})

test_that("final_assignment_cpp runs for classification", {
  A <- X[sample(150, 10), ]
  fa <- final_assignment_cpp(X, A, "gaussian", "euclidean", 1L,
                              1.0, 1.0, 2.0, matrix(0, 0, 0))
  expect_true(!is.null(fa$best_antibody_idx))
  expect_equal(length(fa$best_antibody_idx), 150)
})

# ==== Idiotypic Dynamics ====

test_that("idiotypic_dynamics_cpp runs and returns expected structure", {
  A <- matrix(rnorm(30), nrow = 6, ncol = 5)
  result <- idiotypic_dynamics_cpp(
    A, "gaussian", 1.0, 1.0, 2.0,
    0.01, 0.5, 0.5, 0.1, 0.1, 20L, 0.5
  )
  expect_true(is.list(result))
  expect_true("Ab_Ab_affinity" %in% names(result))
  expect_true("population" %in% names(result))
  expect_true("keep" %in% names(result))
  expect_equal(length(result$population), 6)
  expect_equal(length(result$keep), 6)
  expect_equal(dim(result$Ab_Ab_affinity), c(6, 6))
})
