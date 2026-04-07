# Tests for ConvergentSelector R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("ConvergentSelector initializes with defaults", {
  conv <- ConvergentSelector$new()
  expect_equal(conv$tolerance, 0.5)
  expect_equal(conv$min_appearances, 2)
  expect_null(conv$public_antibodies)
})

test_that("ConvergentSelector initializes with custom params", {
  conv <- ConvergentSelector$new(tolerance = 1.0, min_appearances = 3)
  expect_equal(conv$tolerance, 1.0)
  expect_equal(conv$min_appearances, 3)
})

# ==== find_public ====

test_that("ConvergentSelector find_public with identical repertoires", {
  A1 <- X[1:10, ]
  A2 <- X[1:10, ] + rnorm(40, sd = 0.01)  # near-identical
  conv <- ConvergentSelector$new(tolerance = 0.5, min_appearances = 2)
  public <- conv$find_public(list(A1, A2))
  expect_true(is.matrix(public))
  expect_true(nrow(public) > 0)
  expect_true(!is.null(conv$public_antibodies))
})

test_that("ConvergentSelector find_public with disjoint repertoires", {
  A1 <- matrix(0, nrow = 5, ncol = 4)
  A2 <- matrix(100, nrow = 5, ncol = 4)
  conv <- ConvergentSelector$new(tolerance = 0.1, min_appearances = 2)
  expect_warning(
    public <- conv$find_public(list(A1, A2)),
    "No public"
  )
  # Falls back to reference repertoire
  expect_true(nrow(public) > 0)
})

test_that("ConvergentSelector find_public warns with single repertoire", {
  A <- X[1:5, ]
  conv <- ConvergentSelector$new()
  expect_warning(
    public <- conv$find_public(list(A)),
    "at least 2"
  )
  expect_equal(nrow(public), 5)
})

test_that("ConvergentSelector find_public with ImmuneRepertoire objects", {
  rep1 <- ImmuneRepertoire$new(X[1:10, ])
  rep2 <- ImmuneRepertoire$new(X[1:10, ] + rnorm(40, sd = 0.01))
  conv <- ConvergentSelector$new(tolerance = 0.5, min_appearances = 2)
  public <- conv$find_public(list(rep1, rep2))
  expect_true(is.matrix(public))
  expect_true(nrow(public) > 0)
})

test_that("ConvergentSelector find_public with 3+ repertoires", {
  A1 <- X[1:10, ]
  A2 <- X[1:10, ] + rnorm(40, sd = 0.01)
  A3 <- X[1:10, ] + rnorm(40, sd = 0.01)
  conv <- ConvergentSelector$new(tolerance = 0.5, min_appearances = 3)
  public <- conv$find_public(list(A1, A2, A3))
  expect_true(is.matrix(public))
})

# ==== from_results ====

test_that("ConvergentSelector from_results works with bHIVE results", {
  results <- lapply(seq_len(3), function(i) {
    m <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
    m$fit(X, task = "clustering")
    m$result
  })
  conv <- ConvergentSelector$new(tolerance = 2.0, min_appearances = 2)
  public <- conv$from_results(results)
  expect_true(is.matrix(public))
  expect_true(nrow(public) > 0)
})

# ==== Print ====

test_that("ConvergentSelector print works before and after selection", {
  conv <- ConvergentSelector$new()
  expect_output(print(conv), "ConvergentSelector")

  A1 <- X[1:10, ]
  A2 <- X[1:10, ] + rnorm(40, sd = 0.01)
  conv$find_public(list(A1, A2))
  expect_output(print(conv), "Public antibodies")
})
