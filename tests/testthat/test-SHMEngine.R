# Tests for SHMEngine R6 class

set.seed(42)

# ==== Construction ====

test_that("SHMEngine initializes with valid defaults", {
  shm <- SHMEngine$new()
  expect_equal(shm$method, "uniform")
  expect_true(is.list(shm$params))
  expect_equal(shm$params$decay, 1.0)
  expect_equal(shm$params$mutationMin, 0.01)
})

test_that("SHMEngine initializes with all valid methods", {
  methods <- c("uniform", "airs", "hotspot", "energy", "adaptive")
  for (m in methods) {
    shm <- SHMEngine$new(method = m)
    expect_equal(shm$method, m)
  }
})

test_that("SHMEngine rejects invalid method", {
  expect_error(SHMEngine$new(method = "invalid"))
})

test_that("SHMEngine validates parameter ranges", {
  expect_error(SHMEngine$new(decay = 0), "decay")
  expect_error(SHMEngine$new(decay = 1.5), "decay")
  expect_error(SHMEngine$new(mutationMin = -1), "mutationMin")
  expect_error(SHMEngine$new(temperature = 0), "temperature")
  expect_error(SHMEngine$new(temperature = -1), "temperature")
  expect_error(SHMEngine$new(E_0 = 0), "E_0")
  expect_error(SHMEngine$new(base_rate = 0), "base_rate")
  expect_error(SHMEngine$new(beta1 = 1), "beta1")
  expect_error(SHMEngine$new(beta2 = -0.1), "beta2")
  expect_error(SHMEngine$new(adam_epsilon = 0), "adam_epsilon")
})

test_that("SHMEngine stores custom parameters", {
  shm <- SHMEngine$new(method = "airs", c_rate = 2.0, temperature = 0.3)
  expect_equal(shm$params$c_rate, 2.0)
  expect_equal(shm$params$temperature, 0.3)
})

# ==== State management ====

test_that("SHMEngine init_state creates matrices for adaptive method", {
  shm <- SHMEngine$new(method = "adaptive")
  shm$init_state(5, 3)
  expect_equal(dim(shm$m1_state), c(5, 3))
  expect_equal(dim(shm$m2_state), c(5, 3))
  expect_true(all(shm$m1_state == 0))
  expect_true(all(shm$m2_state == 0))
})

test_that("SHMEngine init_state does nothing for non-adaptive", {
  shm <- SHMEngine$new(method = "uniform")
  shm$init_state(5, 3)
  expect_null(shm$m1_state)
  expect_null(shm$m2_state)
})

test_that("SHMEngine reset_state subsets existing state", {
  shm <- SHMEngine$new(method = "adaptive")
  shm$init_state(5, 3)
  shm$m1_state[2, ] <- c(1, 2, 3)
  shm$reset_state(2, 3, kept_idx = c(2, 4))
  expect_equal(nrow(shm$m1_state), 2)
  expect_equal(shm$m1_state[1, ], c(1, 2, 3))
})

test_that("SHMEngine reset_state creates fresh state when no kept_idx", {
  shm <- SHMEngine$new(method = "adaptive")
  shm$init_state(5, 3)
  shm$reset_state(10, 3)
  expect_equal(dim(shm$m1_state), c(10, 3))
  expect_true(all(shm$m1_state == 0))
})

test_that("SHMEngine reset_state does nothing for non-adaptive", {
  shm <- SHMEngine$new(method = "uniform")
  shm$reset_state(5, 3)
  expect_null(shm$m1_state)
})

# ==== Print ====

test_that("SHMEngine print shows method-specific params", {
  expect_output(print(SHMEngine$new(method = "uniform")), "decay")
  expect_output(print(SHMEngine$new(method = "airs")), "temperature")
  expect_output(print(SHMEngine$new(method = "hotspot")), "base_rate")
  expect_output(print(SHMEngine$new(method = "energy")), "E_0")
  expect_output(print(SHMEngine$new(method = "adaptive")), "beta1")
})

# ==== Integration: SHM with AINet ====

test_that("SHMEngine integrates with AINet for all methods", {
  data(iris)
  X <- as.matrix(iris[, 1:4])

  for (m in c("uniform", "airs", "hotspot", "energy", "adaptive")) {
    shm <- SHMEngine$new(method = m)
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       shm = shm, verbose = FALSE)
    # Even though AINet doesn't yet use the SHM module in its main loop,
    # it should store it correctly
    expect_equal(model$modules$shm$method, m)
  }
})
