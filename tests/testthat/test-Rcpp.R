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
    "gaussian", 1.0, 1.0, 2.0, 0L,
    "uniform", 1.0, 0.5, 1.0, 0.1, 0.9, 0.999, 1e-8,
    matrix(0, 0, 0), matrix(0, 0, 0)
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
    "gaussian", 1.0, 1.0, 2.0, 3L,
    "uniform", 1.0, 0.5, 1.0, 0.1, 0.9, 0.999, 1e-8,
    matrix(0, 0, 0), matrix(0, 0, 0)
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

# ==== Values against R reference implementations ====

set.seed(7)
Xv <- matrix(rnorm(24), nrow = 6)
Av <- matrix(rnorm(12), nrow = 3)
Hx <- matrix(sample(0:2, 24, replace = TRUE), nrow = 6)  # integer-coded
Ha <- matrix(sample(0:2, 12, replace = TRUE), nrow = 3)
Z  <- matrix(0, 0, 0)
S  <- crossprod(matrix(rnorm(40), 10)) / 10 + diag(4)

# n x m matrix of f(X[i, ], A[j, ])
ref_pairs <- function(X, A, f) {
  outer(seq_len(nrow(X)), seq_len(nrow(A)),
        Vectorize(function(i, j) f(X[i, ], A[j, ])))
}
cos_sim <- function(x, a) sum(x * a) / sqrt(sum(x^2) * sum(a^2))

test_that("affinity kernels match their formulas", {
  expect_equal(compute_affinity_matrix(Xv, Av, "gaussian", alpha = 0.7),
               ref_pairs(Xv, Av, function(x, a) exp(-0.7 * sum((x - a)^2))))
  expect_equal(compute_affinity_matrix(Xv, Av, "laplace", alpha = 0.7),
               ref_pairs(Xv, Av, function(x, a) exp(-0.7 * sum(abs(x - a)))))
  expect_equal(compute_affinity_matrix(Xv, Av, "polynomial", c = 0.5, p = 3),
               ref_pairs(Xv, Av, function(x, a)
                 min(1, max(0, (sum(x * a) + 0.5) /
                              max(sum(x^2) + 0.5, sum(a^2) + 0.5)))^3))
  expect_equal(compute_affinity_matrix(Xv, Av, "cosine"),
               ref_pairs(Xv, Av, cos_sim))
  expect_equal(compute_affinity_matrix(Hx, Ha, "hamming"),
               ref_pairs(Hx, Ha, function(x, a) mean(x == a)))
})

test_that("cosine affinity of a zero vector is 0, not NaN", {
  expect_equal(compute_affinity_matrix(rbind(c(0, 0, 0)), rbind(c(1, 2, 3)),
                                       "cosine"), matrix(0))
})

test_that("distances match stats::dist() and stats::mahalanobis()", {
  between <- function(method, ...) {
    as.matrix(dist(rbind(Xv, Av), method = method, ...))[1:6, 7:9,
                                                         drop = FALSE]
  }
  expect_equal(compute_distance_matrix(Xv, Av, "euclidean", 2, Z),
               unname(between("euclidean")))
  expect_equal(compute_distance_matrix(Xv, Av, "manhattan", 2, Z),
               unname(between("manhattan")))
  expect_equal(compute_distance_matrix(Xv, Av, "minkowski", 3, Z),
               unname(between("minkowski", p = 3)))
  # p = 2 takes the BLAS Euclidean path
  expect_equal(compute_distance_matrix(Xv, Av, "minkowski", 2, Z),
               compute_distance_matrix(Xv, Av, "euclidean", 2, Z))
  expect_equal(compute_distance_matrix(Xv, Av, "cosine", 2, Z),
               1 - ref_pairs(Xv, Av, cos_sim))
  expect_equal(compute_distance_matrix(Xv, Av, "mahalanobis", 2, solve(S)),
               sqrt(vapply(seq_len(nrow(Av)), function(j)
                 stats::mahalanobis(Xv, Av[j, ], S), numeric(nrow(Xv)))))
  expect_equal(compute_distance_matrix(Hx, Ha, "hamming", 2, Z),
               ref_pairs(Hx, Ha, function(x, a) sum(x != a)))
})

test_that("pairwise distance equals the full matrix, with an exact zero diagonal", {
  for (type in c("euclidean", "manhattan", "minkowski", "cosine",
                 "mahalanobis")) {
    D <- compute_pairwise_distance(Xv, type, 3, solve(S))
    expect_equal(D, compute_distance_matrix(Xv, Xv, type, 3, solve(S)),
                 label = type)
    # Passing one matrix as both arguments takes BLAS's symmetric path,
    # which used to leave ~6e-8 on the Euclidean diagonal.
    expect_identical(diag(D), rep(0, nrow(Xv)), label = type)
  }
  expect_equal(compute_pairwise_distance(Hx, "hamming", 2, Z),
               compute_distance_matrix(Hx, Hx, "hamming", 2, Z))
})

test_that("unknown metric names are rejected", {
  expect_error(compute_affinity_matrix(Xv, Av, "rbf"),
               "Invalid affinity type")
  expect_error(compute_distance_matrix(Xv, Av, "chebyshev", 2, Z),
               "Invalid distance type")
  expect_error(network_suppression_cpp(Av, "chebyshev", 0.1, 2, Z),
               "Invalid distance type")
})

# ==== Network suppression ====

test_that("suppression removes a near-duplicate under every metric", {
  # Row 2 is row 1 scaled by 1.0001: a near-duplicate in every geometry
  # (its integer-coded value is identical, for Hamming). Row 3 is far from
  # both, and orthogonal to row 1 for cosine.
  a <- c(1, 2, 3, 4)
  A <- rbind(a, a * 1.0001, c(-4, 3, -2, 1))
  for (type in c("euclidean", "manhattan", "minkowski", "cosine",
                 "mahalanobis", "hamming")) {
    expect_identical(network_suppression_cpp(A, type, 0.01, 3, diag(4)),
                     c(TRUE, FALSE, TRUE), label = type)
  }
})

test_that("suppression uses the strict epsilon boundary", {
  A <- rbind(c(0, 0), c(3, 4))  # Euclidean distance exactly 5
  expect_identical(network_suppression_cpp(A, "euclidean", 5, 2, Z),
                   c(TRUE, TRUE))
  expect_identical(network_suppression_cpp(A, "euclidean", 5.0001, 2, Z),
                   c(TRUE, FALSE))
})

test_that("cosine suppression treats zero vectors as maximally distant", {
  A <- rbind(c(0, 0, 0), c(0, 0, 0), c(1, 1, 1))
  expect_identical(network_suppression_cpp(A, "cosine", 0.5, 2, Z),
                   c(TRUE, TRUE, TRUE))
})

# ==== Scalar affinity on the clonal-selection hot path ====

test_that("clonal selection evaluates clones with every affinity kernel", {
  set.seed(8)
  A <- X[sample(150, 6), ]
  Xi <- round(X)  # integer-coded for Hamming
  for (af in c("gaussian", "laplace", "polynomial", "cosine", "hamming")) {
    data <- if (af == "hamming") Xi else X
    res <- clonal_selection_iteration_cpp(
      if (af == "hamming") round(A) else A, data, rep(0, 150), 0L, 3L, 5.0,
      Inf, 1.0, 0.01, 1L, af, 0.5, 1.0, 2.0, 0L,
      "uniform", 1.0, 0.5, 1.0, 0.1, 0.9, 0.999, 1e-8, Z, Z)
    expect_true(all(is.finite(res$A)), label = af)
    expect_gte(nrow(res$A), 6L)
  }
})

test_that("polynomial affinity is bounded and peaks at the data point", {
  # The raw kernel (x.a + c)^p grows without limit in ||a||, so lengthening
  # an antibody along x raised its affinity to x. The bounded kernel must
  # score a == x highest and penalise a longer copy of it.
  set.seed(11)
  x <- matrix(runif(20), 1)
  A <- rbind(x, 3 * x, 0.5 * x, matrix(runif(60), 3))
  for (cc in c(0, 1)) {
    aff <- compute_affinity_matrix(x, A, "polynomial", c = cc, p = 3)
    expect_true(all(aff >= 0 & aff <= 1), label = paste("c =", cc))
    expect_equal(aff[1, 1], 1, label = paste("c =", cc))
    expect_equal(which.max(aff[1, ]), 1L, label = paste("c =", cc))
    expect_lt(aff[1, 2], 0.5)  # three times too long scores well below a == x
  }
})
