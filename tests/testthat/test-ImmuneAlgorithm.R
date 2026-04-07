# Tests for ImmuneAlgorithm R6 class (abstract base)

set.seed(42)

test_that("ImmuneAlgorithm initializes with defaults", {
  algo <- ImmuneAlgorithm$new()
  expect_true(is.list(algo$config))
  expect_true(is.list(algo$modules))
  expect_true(is.list(algo$history))
  expect_null(algo$repertoire)
  expect_null(algo$result)
})

test_that("ImmuneAlgorithm initializes with custom config and modules", {
  algo <- ImmuneAlgorithm$new(
    config = list(nAntibodies = 10, beta = 5),
    modules = list(shm = SHMEngine$new())
  )
  expect_equal(algo$config$nAntibodies, 10)
  expect_true(inherits(algo$modules$shm, "SHMEngine"))
})

test_that("ImmuneAlgorithm$fit is abstract and throws error", {
  algo <- ImmuneAlgorithm$new()
  expect_error(algo$fit(matrix(1)), "abstract")
})

test_that("ImmuneAlgorithm$predict fails when not fitted", {
  algo <- ImmuneAlgorithm$new()
  expect_error(algo$predict(matrix(1)), "not been fitted")
})

test_that("ImmuneAlgorithm$summary returns NULL when no history", {
  algo <- ImmuneAlgorithm$new()
  expect_null(algo$summary())
})

test_that("ImmuneAlgorithm$print works before fitting", {
  algo <- ImmuneAlgorithm$new()
  expect_output(print(algo), "not yet fitted")
})

test_that("ImmuneAlgorithm$print works after fitting (via AINet)", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  model <- AINet$new(nAntibodies = 5, maxIter = 3, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_output(print(model), "AINet")
  expect_output(print(model), "clustering")
})

test_that("ImmuneAlgorithm$summary returns data frame after fitting", {
  data(iris)
  X <- as.matrix(iris[, 1:4])
  model <- AINet$new(nAntibodies = 5, maxIter = 3, verbose = FALSE)
  model$fit(X, task = "clustering")
  s <- model$summary()
  expect_true(is.data.frame(s))
  expect_true("iteration" %in% names(s))
  expect_true("n_antibodies" %in% names(s))
  expect_equal(nrow(s), length(model$history))
})
