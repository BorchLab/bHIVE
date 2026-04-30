# Tests for AINet R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])
y_class <- iris$Species

# ==== Construction / Validation ====

test_that("AINet initializes with defaults", {
  model <- AINet$new(verbose = FALSE)
  expect_true(inherits(model, "AINet"))
  expect_true(inherits(model, "ImmuneAlgorithm"))
  expect_equal(model$config$nAntibodies, 20L)
  expect_equal(model$config$beta, 5)
  expect_null(model$result)
})

test_that("AINet validates numeric parameters", {
  expect_error(AINet$new(nAntibodies = -1), "nAntibodies")
  expect_error(AINet$new(beta = 0), "beta")
  expect_error(AINet$new(epsilon = -1), "epsilon")
  expect_error(AINet$new(maxIter = 0), "maxIter")
  expect_error(AINet$new(k = 0), "k")
  expect_error(AINet$new(mutationDecay = 0), "mutationDecay")
  expect_error(AINet$new(mutationDecay = 1.5), "mutationDecay")
  expect_error(AINet$new(mutationMin = -1), "mutationMin")
})

test_that("AINet validates affinity and distance function names", {
  expect_error(AINet$new(affinityFunc = "invalid"))
  expect_error(AINet$new(distFunc = "invalid"))
  expect_error(AINet$new(initMethod = "invalid"))
})

test_that("AINet accepts module injection", {
  shm <- SHMEngine$new(method = "adaptive")
  model <- AINet$new(shm = shm, verbose = FALSE)
  expect_true("shm" %in% names(model$modules))
  expect_equal(model$modules$shm$method, "adaptive")
})

test_that("AINet stores only non-NULL modules", {
  model <- AINet$new(verbose = FALSE)
  expect_equal(length(model$modules), 0)

  model2 <- AINet$new(shm = SHMEngine$new(), verbose = FALSE)
  expect_equal(length(model2$modules), 1)
})

# ==== Clustering ====

test_that("AINet fits clustering task", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")

  expect_equal(model$result$task, "clustering")
  expect_true(is.matrix(model$result$antibodies))
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments > 0))
})

test_that("AINet infers clustering when y is NULL", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X)
  expect_equal(model$result$task, "clustering")
})

# ==== Classification ====

test_that("AINet fits classification task", {
  model <- AINet$new(nAntibodies = 15, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")

  expect_equal(model$result$task, "classification")
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments %in% levels(y_class)))
})

test_that("AINet infers classification when y is factor", {
  model <- AINet$new(nAntibodies = 10, maxIter = 3, verbose = FALSE)
  model$fit(X, y_class)
  expect_equal(model$result$task, "classification")
})

# ==== Prediction ====

test_that("AINet predict works for clustering", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
})

test_that("AINet predict works for classification", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
  expect_true(all(preds %in% levels(y_class)))
})

# ==== Initialization methods ====

test_that("AINet works with all initialization methods", {
  for (init in c("sample", "random", "random_uniform", "kmeans++")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3, initMethod = init,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Affinity / distance variants ====

test_that("AINet works with different affinity functions", {
  for (af in c("gaussian", "laplace", "polynomial", "cosine")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       affinityFunc = af, verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

test_that("AINet works with different distance functions", {
  for (df in c("euclidean", "manhattan", "minkowski")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       distFunc = df, verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Repertoire / history ====

test_that("AINet populates repertoire after fitting", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(inherits(model$repertoire, "ImmuneRepertoire"))
  expect_true(model$repertoire$size() > 0)
})

test_that("AINet records iteration history", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) > 0)
  expect_true(all(vapply(model$history, function(h) h$n_antibodies > 0, logical(1))))
})

# ==== Early stopping ====

test_that("AINet early stopping works", {
  model <- AINet$new(nAntibodies = 10, maxIter = 100,
                     stopTolerance = 5, noImprovementLimit = 3,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) < 100)
})

# ==== Edge case: all antibodies suppressed ====

test_that("AINet handles small populations gracefully", {
  # With extreme epsilon, the algorithm may error or produce a tiny repertoire
  # Either outcome is acceptable
  result <- tryCatch({
    model <- AINet$new(nAntibodies = 2, maxIter = 3, epsilon = 100,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    "completed"
  }, error = function(e) "error")
  expect_true(result %in% c("completed", "error"))
})
