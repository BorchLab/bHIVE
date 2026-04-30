# Tests for GerminalCenter R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("GerminalCenter initializes with defaults", {
  gc <- GerminalCenter$new()
  expect_equal(gc$nTfh, 10L)
  expect_equal(gc$selectionPressure, 0.5)
  expect_equal(gc$rounds, 1L)
})

test_that("GerminalCenter validates parameters", {
  expect_error(GerminalCenter$new(nTfh = 0), "nTfh")
  expect_error(GerminalCenter$new(selectionPressure = -0.1), "selectionPressure")
  expect_error(GerminalCenter$new(selectionPressure = 1.5), "selectionPressure")
  expect_error(GerminalCenter$new(rounds = 0), "rounds")
})

# ==== Selection ====

test_that("GerminalCenter select reduces repertoire size (clustering)", {
  A <- X[sample(150, 30), ]
  rep <- ImmuneRepertoire$new(A)
  gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.8)
  gc$select(rep, X, task = "clustering")
  expect_true(rep$size() < 30)
  expect_true(rep$size() > 0)
  expect_equal(rep$size(), nrow(rep$metadata))
})

test_that("GerminalCenter select works for classification", {
  A <- X[sample(150, 20), ]
  rep <- ImmuneRepertoire$new(A)
  gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.5)
  gc$select(rep, X, iris$Species, "classification")
  expect_true(rep$size() <= 20)
  expect_true(rep$size() > 0)
})

test_that("GerminalCenter does nothing when repertoire <= nTfh", {
  A <- X[sample(150, 5), ]
  rep <- ImmuneRepertoire$new(A)
  gc <- GerminalCenter$new(nTfh = 10)
  gc$select(rep, X, task = "clustering")
  expect_equal(rep$size(), 5)  # unchanged
})

test_that("GerminalCenter multiple rounds reduces more", {
  set.seed(42)
  A <- X[sample(150, 30), ]
  rep1 <- ImmuneRepertoire$new(A)
  gc1 <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.5, rounds = 1)
  gc1$select(rep1, X, task = "clustering")
  size1 <- rep1$size()

  set.seed(42)
  rep3 <- ImmuneRepertoire$new(A)
  gc3 <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.5, rounds = 3)
  gc3$select(rep3, X, task = "clustering")
  size3 <- rep3$size()

  expect_true(size3 <= size1)
})

test_that("GerminalCenter selectionPressure=0 keeps most antibodies", {
  A <- X[sample(150, 20), ]
  rep <- ImmuneRepertoire$new(A)
  gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0)
  gc$select(rep, X, task = "clustering")
  expect_equal(rep$size(), 20)  # all survive
})

# ==== Print ====

test_that("GerminalCenter print works", {
  gc <- GerminalCenter$new(nTfh = 8, selectionPressure = 0.3)
  expect_output(print(gc), "GerminalCenter")
  expect_output(print(gc), "0.30")
})
