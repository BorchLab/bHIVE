# Tests for MemoryPool R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("MemoryPool initializes with defaults", {
  mp <- MemoryPool$new()
  expect_equal(mp$archive_threshold, 0.5)
  expect_equal(mp$max_memory, 100)
  expect_equal(mp$recall_threshold, 0.3)
  expect_equal(mp$size(), 0)
  expect_null(mp$memory_cells)
  expect_null(mp$memory_metadata)
})

test_that("MemoryPool initializes with custom params", {
  mp <- MemoryPool$new(archive_threshold = 0.1, max_memory = 50,
                        recall_threshold = 0.05)
  expect_equal(mp$archive_threshold, 0.1)
  expect_equal(mp$max_memory, 50)
})

# ==== Archive ====

test_that("MemoryPool archive stores high-affinity antibodies", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  mp <- MemoryPool$new(archive_threshold = 0.01)  # low threshold to archive
  n_archived <- mp$archive(rep, X)
  expect_true(n_archived > 0)
  expect_true(mp$size() > 0)
  expect_true(is.matrix(mp$memory_cells))
  expect_true(is.data.frame(mp$memory_metadata))
})

test_that("MemoryPool archive returns 0 when no antibodies meet threshold", {
  A <- matrix(100, nrow = 3, ncol = 4)  # far from data
  rep <- ImmuneRepertoire$new(A)
  mp <- MemoryPool$new(archive_threshold = 0.99)
  n_archived <- mp$archive(rep, X)
  expect_equal(n_archived, 0)
  expect_equal(mp$size(), 0)
})

test_that("MemoryPool archive accumulates across calls", {
  mp <- MemoryPool$new(archive_threshold = 0.01)

  A1 <- X[sample(150, 5), ]
  rep1 <- ImmuneRepertoire$new(A1)
  mp$archive(rep1, X)
  size1 <- mp$size()

  A2 <- X[sample(150, 5), ]
  rep2 <- ImmuneRepertoire$new(A2)
  mp$archive(rep2, X)
  size2 <- mp$size()

  expect_true(size2 >= size1)
})

test_that("MemoryPool archive respects max_memory", {
  mp <- MemoryPool$new(archive_threshold = 0.001, max_memory = 5)
  A <- X[sample(150, 20), ]
  rep <- ImmuneRepertoire$new(A)
  mp$archive(rep, X)
  expect_true(mp$size() <= 5)
})

test_that("MemoryPool archive sets state to memory", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  mp <- MemoryPool$new(archive_threshold = 0.01)
  mp$archive(rep, X)
  expect_true(all(mp$memory_metadata$state == "memory"))
})

# ==== Recall ====

test_that("MemoryPool recall returns relevant memories", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  mp <- MemoryPool$new(archive_threshold = 0.01, recall_threshold = 0.01)
  mp$archive(rep, X)
  recalled <- mp$recall(X[1:10, ])
  expect_true(is.matrix(recalled))
  expect_true(nrow(recalled) > 0)
  expect_equal(ncol(recalled), 4)
})

test_that("MemoryPool recall returns empty matrix when no memories", {
  mp <- MemoryPool$new()
  recalled <- mp$recall(X[1:10, ])
  expect_true(is.matrix(recalled))
  expect_equal(nrow(recalled), 0)
})

test_that("MemoryPool recall returns empty when no memories meet threshold", {
  A <- X[sample(150, 5), ]
  rep <- ImmuneRepertoire$new(A)
  mp <- MemoryPool$new(archive_threshold = 0.01, recall_threshold = 0.99)
  mp$archive(rep, X)
  # Very high recall threshold
  recalled <- mp$recall(matrix(100, nrow = 5, ncol = 4))
  expect_equal(nrow(recalled), 0)
})

# ==== Size ====

test_that("MemoryPool size tracks correctly", {
  mp <- MemoryPool$new(archive_threshold = 0.01)
  expect_equal(mp$size(), 0)

  A <- X[sample(150, 5), ]
  rep <- ImmuneRepertoire$new(A)
  mp$archive(rep, X)
  expect_true(mp$size() > 0)
})

# ==== Print ====

test_that("MemoryPool print works", {
  mp <- MemoryPool$new()
  expect_output(print(mp), "MemoryPool")
  expect_output(print(mp), "0 cells")
})
