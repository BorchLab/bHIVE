# Tests for Microenvironment R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])

# ==== Construction ====

test_that("Microenvironment initializes with defaults", {
  me <- Microenvironment$new()
  expect_null(me$density_bandwidth)
  expect_equal(me$high_density_threshold, 0.75)
  expect_equal(me$low_density_threshold, 0.25)
  expect_equal(me$stabilization_factor, 0.3)
  expect_equal(me$exploration_factor, 2.0)
  expect_null(me$last_densities)
  expect_null(me$last_zones)
})

test_that("Microenvironment initializes with custom params", {
  me <- Microenvironment$new(
    high_density_threshold = 0.9,
    low_density_threshold = 0.1,
    stabilization_factor = 0.5,
    exploration_factor = 3.0
  )
  expect_equal(me$high_density_threshold, 0.9)
  expect_equal(me$exploration_factor, 3.0)
})

# ==== Assessment ====

test_that("Microenvironment assess returns correct structure", {
  A <- X[sample(150, 15), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new()
  env <- me$assess(rep, X)

  expect_true(is.list(env))
  expect_named(env, c("densities", "zones", "mutation_modifiers", "gradients"))
  expect_equal(length(env$densities), 15)
  expect_equal(length(env$zones), 15)
  expect_equal(length(env$mutation_modifiers), 15)
  expect_equal(dim(env$gradients), c(15, 4))
})

test_that("Microenvironment classifies zones correctly", {
  A <- X[sample(150, 20), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new()
  env <- me$assess(rep, X)

  expect_true(all(env$zones %in% c("stable", "explore", "boundary")))
  # With 20 antibodies at 25/75 thresholds, we expect some of each
  expect_true(any(env$zones == "stable"))
  expect_true(any(env$zones == "explore"))
})

test_that("Microenvironment mutation_modifiers match zones", {
  A <- X[sample(150, 15), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new(stabilization_factor = 0.3, exploration_factor = 2.0)
  env <- me$assess(rep, X)

  stable_idx <- which(env$zones == "stable")
  if (length(stable_idx) > 0) {
    expect_true(all(env$mutation_modifiers[stable_idx] == 0.3))
  }
  explore_idx <- which(env$zones == "explore")
  if (length(explore_idx) > 0) {
    expect_true(all(env$mutation_modifiers[explore_idx] == 2.0))
  }
  boundary_idx <- which(env$zones == "boundary")
  if (length(boundary_idx) > 0) {
    expect_true(all(env$mutation_modifiers[boundary_idx] == 1.0))
  }
})

test_that("Microenvironment updates repertoire metadata states", {
  A <- X[sample(150, 15), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new()
  env <- me$assess(rep, X)

  stable_idx <- which(env$zones == "stable")
  if (length(stable_idx) > 0) {
    expect_true(all(rep$metadata$state[stable_idx] == "memory"))
  }
  explore_idx <- which(env$zones == "explore")
  if (length(explore_idx) > 0) {
    expect_true(all(rep$metadata$state[explore_idx] == "activated"))
  }
})

test_that("Microenvironment stores last_densities and last_zones", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new()
  me$assess(rep, X)

  expect_true(!is.null(me$last_densities))
  expect_true(!is.null(me$last_zones))
  expect_equal(length(me$last_densities), 10)
  expect_equal(length(me$last_zones), 10)
})

test_that("Microenvironment gradients are finite", {
  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  me <- Microenvironment$new()
  env <- me$assess(rep, X)
  expect_true(all(is.finite(env$gradients)))
})

# ==== Print ====

test_that("Microenvironment print works before and after assessment", {
  me <- Microenvironment$new()
  expect_output(print(me), "Microenvironment")

  A <- X[sample(150, 10), ]
  rep <- ImmuneRepertoire$new(A)
  me$assess(rep, X)
  expect_output(print(me), "assessment")
})
