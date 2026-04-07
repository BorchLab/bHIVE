# Tests for ActivationGate R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("ActivationGate initializes with defaults", {
  gate <- ActivationGate$new()
  expect_equal(gate$signal2_type, "density")
  expect_equal(gate$threshold1, 0.1)
  expect_equal(gate$threshold2, 0.3)
  expect_null(gate$danger_signals)
})

test_that("ActivationGate accepts all signal types", {
  for (s2 in c("density", "danger", "entropy")) {
    gate <- ActivationGate$new(signal2_type = s2)
    expect_equal(gate$signal2_type, s2)
  }
})

test_that("ActivationGate rejects invalid signal type", {
  expect_error(ActivationGate$new(signal2_type = "invalid"))
})

# ==== Evaluate ====

test_that("ActivationGate evaluate returns correct dimensions", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  gate <- ActivationGate$new()
  aff <- rep$affinity_matrix(X, "gaussian")
  result <- gate$evaluate(aff, X, A)

  expect_true(is.logical(result))
  expect_equal(dim(result), c(150, 10))
})

test_that("ActivationGate density signal produces activations", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  gate <- ActivationGate$new(signal2_type = "density", threshold1 = 0.01,
                              threshold2 = 0.1)
  aff <- rep$affinity_matrix(X, "gaussian")
  result <- gate$evaluate(aff, X, A)
  expect_true(sum(result) > 0)
})

test_that("ActivationGate high thresholds reduce activations", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  aff <- rep$affinity_matrix(X, "gaussian")

  gate_low <- ActivationGate$new(threshold1 = 0.01, threshold2 = 0.1)
  gate_high <- ActivationGate$new(threshold1 = 0.9, threshold2 = 0.9)

  act_low <- sum(gate_low$evaluate(aff, X, A))
  act_high <- sum(gate_high$evaluate(aff, X, A))

  expect_true(act_low >= act_high)
})

test_that("ActivationGate danger signal works", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  danger <- runif(150)
  gate <- ActivationGate$new(signal2_type = "danger", threshold2 = 0.5,
                              danger_signals = danger)
  aff <- rep$affinity_matrix(X, "gaussian")
  result <- gate$evaluate(aff, X, A)
  expect_true(is.logical(result))
})

test_that("ActivationGate danger signal errors without danger_signals", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  gate <- ActivationGate$new(signal2_type = "danger")
  aff <- rep$affinity_matrix(X, "gaussian")
  expect_error(gate$evaluate(aff, X, A), "danger_signals")
})

test_that("ActivationGate entropy signal works for classification", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  gate <- ActivationGate$new(signal2_type = "entropy", threshold2 = 1.0)
  aff <- rep$affinity_matrix(X, "gaussian")
  result <- gate$evaluate(aff, X, A, y = iris$Species, task = "classification")
  expect_true(is.logical(result))
  expect_equal(dim(result), c(150, 10))
})

test_that("ActivationGate entropy falls back to density for non-classification", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  gate <- ActivationGate$new(signal2_type = "entropy", threshold1 = 0.01,
                              threshold2 = 0.1)
  aff <- rep$affinity_matrix(X, "gaussian")
  result <- gate$evaluate(aff, X, A, task = "clustering")
  expect_true(is.logical(result))
  expect_true(sum(result) > 0)
})

# ==== Print ====

test_that("ActivationGate print works", {
  gate <- ActivationGate$new(signal2_type = "density", threshold1 = 0.2)
  expect_output(print(gate), "ActivationGate")
  expect_output(print(gate), "density")
})
