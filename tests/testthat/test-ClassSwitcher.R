# Tests for ClassSwitcher R6 class

set.seed(42)

# ==== Construction ====

test_that("ClassSwitcher initializes with defaults", {
  cs <- ClassSwitcher$new()
  expect_equal(cs$alpha_IgM, 0.1)
  expect_equal(cs$alpha_IgG, 5.0)
  expect_equal(cs$alpha_IgA, 1.0)
})

test_that("ClassSwitcher initializes with custom params", {
  cs <- ClassSwitcher$new(alpha_IgM = 0.5, alpha_IgG = 10, alpha_IgA = 2.0)
  expect_equal(cs$alpha_IgM, 0.5)
  expect_equal(cs$alpha_IgG, 10)
  expect_equal(cs$alpha_IgA, 2.0)
})

# ==== Switch isotypes ====

test_that("ClassSwitcher switch_isotypes returns correct alphas", {
  A <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  cs <- ClassSwitcher$new(alpha_IgM = 0.1, alpha_IgG = 5.0, alpha_IgA = 1.0)
  zones <- c("stable", "explore", "boundary", "stable", "explore",
             "boundary", "stable", "explore", "boundary", "stable")
  alphas <- cs$switch_isotypes(rep, zones)

  expect_equal(length(alphas), 10)
  # stable -> IgG (5.0)
  expect_equal(alphas[1], 5.0)
  # explore -> IgM (0.1)
  expect_equal(alphas[2], 0.1)
  # boundary -> IgA (1.0)
  expect_equal(alphas[3], 1.0)
})

test_that("ClassSwitcher switch_isotypes updates repertoire metadata", {
  A <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  cs <- ClassSwitcher$new()
  zones <- c("stable", "explore", "boundary", "stable", "explore", "boundary")
  cs$switch_isotypes(rep, zones)

  expect_equal(rep$metadata$isotype[1], "IgG")
  expect_equal(rep$metadata$isotype[2], "IgM")
  expect_equal(rep$metadata$isotype[3], "IgA")
})

test_that("ClassSwitcher all stable -> all IgG", {
  A <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  cs <- ClassSwitcher$new()
  alphas <- cs$switch_isotypes(rep, rep("stable", 5))
  expect_true(all(alphas == 5.0))
  expect_true(all(rep$metadata$isotype == "IgG"))
})

test_that("ClassSwitcher all explore -> all IgM", {
  A <- matrix(rnorm(25), nrow = 5, ncol = 5)
  rep <- ImmuneRepertoire$new(A)
  cs <- ClassSwitcher$new()
  alphas <- cs$switch_isotypes(rep, rep("explore", 5))
  expect_true(all(alphas == 0.1))
  expect_true(all(rep$metadata$isotype == "IgM"))
})

# ==== get_alpha ====

test_that("ClassSwitcher get_alpha returns correct values", {
  cs <- ClassSwitcher$new(alpha_IgM = 0.1, alpha_IgG = 5.0, alpha_IgA = 1.0)
  expect_equal(cs$get_alpha("IgM"), 0.1)
  expect_equal(cs$get_alpha("IgG"), 5.0)
  expect_equal(cs$get_alpha("IgA"), 1.0)
})

test_that("ClassSwitcher get_alpha returns IgA for unknown isotype", {
  cs <- ClassSwitcher$new()
  expect_equal(cs$get_alpha("IgE"), 1.0)
})

# ==== Print ====

test_that("ClassSwitcher print works", {
  cs <- ClassSwitcher$new()
  expect_output(print(cs), "ClassSwitcher")
  expect_output(print(cs), "IgM")
  expect_output(print(cs), "IgG")
  expect_output(print(cs), "IgA")
})
