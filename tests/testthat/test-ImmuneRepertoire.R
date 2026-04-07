# Tests for ImmuneRepertoire R6 class

set.seed(42)

# ---- Setup ----
X <- matrix(rnorm(100), nrow = 20, ncol = 5)
A <- matrix(rnorm(50), nrow = 10, ncol = 5)

# ==== Construction ====

test_that("ImmuneRepertoire initializes with correct dimensions", {
  rep <- ImmuneRepertoire$new(A)
  expect_equal(rep$size(), 10)
  expect_equal(rep$n_features(), 5)
  expect_equal(dim(rep$cells), c(10, 5))
})

test_that("ImmuneRepertoire creates default metadata", {
  rep <- ImmuneRepertoire$new(A)
  expect_true(is.data.frame(rep$metadata))
  expect_equal(nrow(rep$metadata), 10)
  expect_true(all(rep$metadata$isotype == "IgM"))
  expect_true(all(rep$metadata$state == "naive"))
  expect_true(all(rep$metadata$age == 0L))
  expect_true(all(is.na(rep$metadata$lineage)))
})

test_that("ImmuneRepertoire accepts custom metadata", {
  meta <- data.frame(
    isotype = rep("IgG", 10),
    state = rep("memory", 10),
    age = seq_len(10),
    lineage = paste0("L", seq_len(10)),
    stringsAsFactors = FALSE
  )
  rep <- ImmuneRepertoire$new(A, metadata = meta)
  expect_equal(rep$metadata$isotype[1], "IgG")
  expect_equal(rep$metadata$age[5], 5L)
})

test_that("ImmuneRepertoire rejects non-numeric input", {
  bad <- matrix("a", nrow = 5, ncol = 3)
  expect_error(ImmuneRepertoire$new(bad))
})

test_that("ImmuneRepertoire rejects metadata with wrong row count", {
  bad_meta <- data.frame(isotype = "IgM", state = "naive",
                         age = 0L, lineage = NA_character_,
                         stringsAsFactors = FALSE)
  expect_error(ImmuneRepertoire$new(A, metadata = bad_meta))
})

# ==== Affinity / Distance ====

test_that("affinity_matrix returns correct dimensions", {
  rep <- ImmuneRepertoire$new(A)
  aff <- rep$affinity_matrix(X, "gaussian")
  expect_equal(dim(aff), c(20, 10))
  expect_true(all(aff >= 0))
  expect_true(all(aff <= 1))
})

test_that("affinity_matrix works with different methods", {
  rep <- ImmuneRepertoire$new(A)
  for (method in c("gaussian", "laplace", "polynomial", "cosine")) {
    aff <- rep$affinity_matrix(X, method)
    expect_equal(dim(aff), c(20, 10))
    expect_true(all(is.finite(aff)))
  }
})

test_that("distance_matrix returns correct dimensions", {
  rep <- ImmuneRepertoire$new(A)
  dist <- rep$distance_matrix(X, "euclidean")
  expect_equal(dim(dist), c(20, 10))
  expect_true(all(dist >= 0))
})

test_that("distance_matrix works with different methods", {
  rep <- ImmuneRepertoire$new(A)
  for (method in c("euclidean", "manhattan")) {
    dist <- rep$distance_matrix(X, method)
    expect_equal(dim(dist), c(20, 10))
    expect_true(all(is.finite(dist)))
  }
})

# ==== Suppression ====

test_that("suppress reduces antibody count", {
  A_close <- rbind(A, A + 0.001)  # near-duplicates
  rep <- ImmuneRepertoire$new(A_close)
  initial_size <- rep$size()
  rep$suppress(epsilon = 0.1, method = "euclidean")
  expect_true(rep$size() <= initial_size)
  expect_equal(rep$size(), nrow(rep$metadata))  # metadata synced
})

test_that("suppress preserves distinct antibodies", {
  rep <- ImmuneRepertoire$new(A)
  rep$suppress(epsilon = 0.001, method = "euclidean")
  expect_equal(rep$size(), 10)  # very small epsilon => nothing suppressed
})

# ==== Subset ====

test_that("subset correctly reduces repertoire", {
  rep <- ImmuneRepertoire$new(A)
  rep$subset(c(1, 3, 5))
  expect_equal(rep$size(), 3)
  expect_equal(nrow(rep$metadata), 3)
})

# ==== Add ====

test_that("add increases repertoire size", {
  rep <- ImmuneRepertoire$new(A)
  new_cells <- matrix(rnorm(10), nrow = 2, ncol = 5)
  rep$add(new_cells)
  expect_equal(rep$size(), 12)
  expect_equal(nrow(rep$metadata), 12)
})

test_that("add with custom metadata", {
  rep <- ImmuneRepertoire$new(A)
  new_cells <- matrix(rnorm(5), nrow = 1, ncol = 5)
  new_meta <- data.frame(
    isotype = "IgG", state = "memory", age = 5L, lineage = "clone1",
    stringsAsFactors = FALSE
  )
  rep$add(new_cells, new_meta)
  expect_equal(rep$metadata$isotype[11], "IgG")
})

test_that("add rejects wrong number of columns", {
  rep <- ImmuneRepertoire$new(A)
  bad <- matrix(1, nrow = 1, ncol = 3)
  expect_error(rep$add(bad))
})

# ==== Age ====

test_that("age_all increments all antibody ages", {
  rep <- ImmuneRepertoire$new(A)
  rep$age_all()
  expect_true(all(rep$metadata$age == 1))
  rep$age_all()
  expect_true(all(rep$metadata$age == 2))
})

# ==== as_matrix ====

test_that("as_matrix returns the cells matrix", {
  rep <- ImmuneRepertoire$new(A)
  expect_identical(rep$as_matrix(), rep$cells)
})

# ==== print ====

test_that("print runs without error", {
  rep <- ImmuneRepertoire$new(A)
  expect_output(print(rep), "ImmuneRepertoire")
})
