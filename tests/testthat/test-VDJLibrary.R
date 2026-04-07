# Tests for VDJLibrary R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("VDJLibrary initializes with defaults", {
  vdj <- VDJLibrary$new()
  expect_equal(vdj$nV, 10)
  expect_equal(vdj$nD, 5)
  expect_equal(vdj$nJ, 4)
  expect_equal(vdj$method, "pca")
  expect_null(vdj$library)
})

test_that("VDJLibrary accepts all valid methods", {
  for (m in c("pca", "cluster", "random_partition")) {
    vdj <- VDJLibrary$new(method = m)
    expect_equal(vdj$method, m)
  }
})

test_that("VDJLibrary rejects invalid method", {
  expect_error(VDJLibrary$new(method = "invalid"))
})

# ==== Build ====

test_that("VDJLibrary build creates library structure", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  vdj$build(X)
  expect_true(!is.null(vdj$library))
  expect_true(all(c("V", "D", "J", "V_dims", "D_dims", "J_dims") %in%
                    names(vdj$library)))
  expect_true(is.matrix(vdj$library$V))
  expect_true(is.matrix(vdj$library$D))
  expect_true(is.matrix(vdj$library$J))
})

test_that("VDJLibrary build with PCA method", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3, method = "pca")
  vdj$build(X)
  # Dimensions should cover all columns
  all_dims <- c(vdj$library$V_dims, vdj$library$D_dims, vdj$library$J_dims)
  expect_true(all(seq_len(4) %in% all_dims))
})

test_that("VDJLibrary build with cluster method", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3, method = "cluster")
  vdj$build(X)
  expect_true(!is.null(vdj$library))
})

test_that("VDJLibrary build with random_partition method", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3, method = "random_partition")
  vdj$build(X)
  expect_true(!is.null(vdj$library))
})

# ==== Generate ====

test_that("VDJLibrary generate returns correct dimensions", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  A <- vdj$generate(20, X)
  expect_equal(dim(A), c(20, 4))
  expect_true(all(is.finite(A)))
})

test_that("VDJLibrary generate auto-builds library", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  expect_null(vdj$library)
  A <- vdj$generate(10, X)
  expect_true(!is.null(vdj$library))
  expect_equal(nrow(A), 10)
})

test_that("VDJLibrary generate with all methods", {
  for (m in c("pca", "cluster", "random_partition")) {
    vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3, method = m)
    A <- vdj$generate(15, X)
    expect_equal(dim(A), c(15, 4))
    expect_true(all(is.finite(A)))
  }
})

test_that("VDJLibrary generate produces diverse antibodies", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  A <- vdj$generate(20, X)
  # Not all rows identical
  expect_true(nrow(unique(A)) > 1)
})

# ==== Edge cases ====

test_that("VDJLibrary works with 2-dimensional data", {
  X2 <- X[, 1:2]
  vdj <- VDJLibrary$new(nV = 3, nD = 2, nJ = 2)
  A <- vdj$generate(10, X2)
  expect_equal(dim(A), c(10, 2))
})

test_that("VDJLibrary works with 1-dimensional data", {
  X1 <- X[, 1, drop = FALSE]
  vdj <- VDJLibrary$new(nV = 3, nD = 2, nJ = 2)
  A <- vdj$generate(10, X1)
  expect_equal(dim(A), c(10, 1))
})

# ==== Print ====

test_that("VDJLibrary print works before build", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  expect_output(print(vdj), "VDJLibrary")
})

test_that("VDJLibrary print works after build", {
  vdj <- VDJLibrary$new(nV = 5, nD = 3, nJ = 3)
  vdj$build(X)
  expect_output(print(vdj), "Combinatorial space")
})
