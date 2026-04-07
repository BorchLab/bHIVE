# Tests for IdiotypicNetwork R6 class

set.seed(42)

# ==== Construction ====

test_that("IdiotypicNetwork initializes with valid defaults", {
  idi <- IdiotypicNetwork$new()
  expect_equal(idi$theta_low, 0.01)
  expect_equal(idi$theta_high, 0.5)
  expect_equal(idi$source_rate, 0.5)
  expect_equal(idi$decay_rate, 0.1)
  expect_equal(idi$dt, 0.1)
  expect_equal(idi$timeSteps, 20L)
  expect_equal(idi$survival_threshold, 0.5)
  expect_null(idi$last_dynamics)
})

test_that("IdiotypicNetwork validates parameters", {
  expect_error(IdiotypicNetwork$new(theta_low = -1), "theta_low")
  expect_error(IdiotypicNetwork$new(theta_low = 0.5, theta_high = 0.3), "theta_high")
  expect_error(IdiotypicNetwork$new(source_rate = -1), "source_rate")
  expect_error(IdiotypicNetwork$new(decay_rate = -1), "decay_rate")
  expect_error(IdiotypicNetwork$new(dt = 0), "dt")
  expect_error(IdiotypicNetwork$new(timeSteps = 0), "timeSteps")
  expect_error(IdiotypicNetwork$new(survival_threshold = -1), "survival_threshold")
})

# ==== Regulation ====

test_that("IdiotypicNetwork regulate runs on a repertoire", {
  A <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new(theta_low = 0.01, theta_high = 0.5)
  idi$regulate(rep, "gaussian", list(alpha = 0.5))

  expect_true(!is.null(idi$last_dynamics))
  expect_true(rep$size() <= 10)  # some may be removed
  expect_true(rep$size() > 0 || !is.null(idi$last_dynamics))
})

test_that("IdiotypicNetwork regulate modifies repertoire in place", {
  A <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  initial_size <- rep$size()
  idi <- IdiotypicNetwork$new(theta_low = 0.001, theta_high = 0.99,
                               survival_threshold = 0.1)
  idi$regulate(rep, "gaussian", list(alpha = 1))
  # Regulation happened; metadata should be in sync

  expect_equal(rep$size(), nrow(rep$metadata))
})

test_that("IdiotypicNetwork get_network returns matrix after regulation", {
  A <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new()
  expect_null(idi$get_network())

  idi$regulate(rep, "gaussian")
  net <- idi$get_network()
  expect_true(is.matrix(net))
  expect_equal(nrow(net), 6)
  expect_equal(ncol(net), 6)
})

test_that("IdiotypicNetwork get_population returns vector after regulation", {
  A <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new()
  expect_null(idi$get_population())

  idi$regulate(rep, "gaussian")
  pop <- idi$get_population()
  expect_true(is.numeric(pop))
  expect_equal(length(pop), 6)
})

test_that("IdiotypicNetwork warns when all antibodies removed", {
  # Create very similar antibodies with strict thresholds
  A <- matrix(0.5, nrow = 5, ncol = 3) + matrix(rnorm(15, sd = 0.01), 5, 3)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new(theta_low = 0.01, theta_high = 0.5,
                               survival_threshold = 10)
  expect_warning(idi$regulate(rep, "gaussian", list(alpha = 1)),
                 "removed all")
})

# ==== Print ====

test_that("IdiotypicNetwork print works before regulation", {
  idi <- IdiotypicNetwork$new()
  expect_output(print(idi), "IdiotypicNetwork")
})

test_that("IdiotypicNetwork print works after regulation", {
  A <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  idi <- IdiotypicNetwork$new()
  idi$regulate(rep, "gaussian")
  expect_output(print(idi), "survived")
})
