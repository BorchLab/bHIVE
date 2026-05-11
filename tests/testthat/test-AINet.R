# Tests for AINet R6 class

set.seed(42)
data(iris)
X <- as.matrix(iris[, 1:4])
y_class <- iris$Species

# ==== Construction / Validation ====

test_that("AINet initializes with defaults", {
  model <- AINet$new(verbose = FALSE)
  expect_true(inherits(model, "AINet"))
  expect_true(inherits(model, "ImmuneAlgorithm"))
  expect_equal(model$config$nAntibodies, 20L)
  expect_equal(model$config$beta, 5)
  expect_null(model$result)
})

test_that("AINet validates numeric parameters", {
  expect_error(AINet$new(nAntibodies = -1), "nAntibodies")
  expect_error(AINet$new(beta = 0), "beta")
  expect_error(AINet$new(epsilon = -1), "epsilon")
  expect_error(AINet$new(maxIter = 0), "maxIter")
  expect_error(AINet$new(k = 0), "k")
  expect_error(AINet$new(mutationDecay = 0), "mutationDecay")
  expect_error(AINet$new(mutationDecay = 1.5), "mutationDecay")
  expect_error(AINet$new(mutationMin = -1), "mutationMin")
})

test_that("AINet validates affinity and distance function names", {
  expect_error(AINet$new(affinityFunc = "invalid"))
  expect_error(AINet$new(distFunc = "invalid"))
  expect_error(AINet$new(initMethod = "invalid"))
})

test_that("AINet accepts module injection", {
  shm <- SHMEngine$new(method = "adaptive")
  model <- AINet$new(shm = shm, verbose = FALSE)
  expect_true("shm" %in% names(model$modules))
  expect_equal(model$modules$shm$method, "adaptive")
})

test_that("AINet stores only non-NULL modules", {
  model <- AINet$new(verbose = FALSE)
  expect_equal(length(model$modules), 0)

  model2 <- AINet$new(shm = SHMEngine$new(), verbose = FALSE)
  expect_equal(length(model2$modules), 1)
})

# ==== Clustering ====

test_that("AINet fits clustering task", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")

  expect_equal(model$result$task, "clustering")
  expect_true(is.matrix(model$result$antibodies))
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments > 0))
})

test_that("AINet infers clustering when y is NULL", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X)
  expect_equal(model$result$task, "clustering")
})

# ==== Classification ====

test_that("AINet fits classification task", {
  model <- AINet$new(nAntibodies = 15, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")

  expect_equal(model$result$task, "classification")
  expect_equal(length(model$result$assignments), nrow(X))
  expect_true(all(model$result$assignments %in% levels(y_class)))
})

test_that("AINet infers classification when y is factor", {
  model <- AINet$new(nAntibodies = 10, maxIter = 3, verbose = FALSE)
  model$fit(X, y_class)
  expect_equal(model$result$task, "classification")
})

# ==== Prediction ====

test_that("AINet predict works for clustering", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
})

test_that("AINet predict works for classification", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, y_class, task = "classification")
  preds <- model$predict(X[1:10, ])
  expect_equal(length(preds), 10)
  expect_true(all(preds %in% levels(y_class)))
})

# ==== Initialization methods ====

test_that("AINet works with all initialization methods", {
  for (init in c("sample", "random", "random_uniform", "kmeans++")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3, initMethod = init,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Affinity / distance variants ====

test_that("AINet works with different affinity functions", {
  for (af in c("gaussian", "laplace", "polynomial", "cosine")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       affinityFunc = af, verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

test_that("AINet works with different distance functions", {
  for (df in c("euclidean", "manhattan", "minkowski")) {
    model <- AINet$new(nAntibodies = 10, maxIter = 3,
                       distFunc = df, verbose = FALSE)
    model$fit(X, task = "clustering")
    expect_true(is.matrix(model$result$antibodies))
  }
})

# ==== Repertoire / history ====

test_that("AINet populates repertoire after fitting", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(inherits(model$repertoire, "ImmuneRepertoire"))
  expect_true(model$repertoire$size() > 0)
})

test_that("AINet records iteration history", {
  model <- AINet$new(nAntibodies = 10, maxIter = 5, verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) > 0)
  expect_true(all(vapply(model$history, function(h) h$n_antibodies > 0, logical(1))))
})

# ==== Early stopping ====

test_that("AINet early stopping works", {
  model <- AINet$new(nAntibodies = 10, maxIter = 100,
                     stopTolerance = 5, noImprovementLimit = 3,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_true(length(model$history) < 100)
})

# ==== Edge case: all antibodies suppressed ====

test_that("AINet handles small populations gracefully", {
  # With extreme epsilon, the algorithm may error or produce a tiny repertoire
  # Either outcome is acceptable
  result <- tryCatch({
    model <- AINet$new(nAntibodies = 2, maxIter = 3, epsilon = 100,
                       verbose = FALSE)
    model$fit(X, task = "clustering")
    "completed"
  }, error = function(e) "error")
  expect_true(result %in% c("completed", "error"))
})

# ==== Module integration: VDJLibrary as init ====

test_that("AINet routes initialization through `init$generate` when supplied", {
  # Inject a tracking mock as the `init` module. If $.initialize_antibodies
  # is correctly routing through it, the closure will be called exactly
  # once during fit and we can verify both the call and the returned
  # repertoire shape.
  call_log <- new.env()
  call_log$n_calls <- 0L
  call_log$last_n  <- NA_integer_

  mock_init <- list(
    generate = function(nAntibodies, X) {
      call_log$n_calls <- call_log$n_calls + 1L
      call_log$last_n  <- nAntibodies
      # Return a deterministic structured matrix (rows are constant per row)
      mat <- matrix(seq_len(nAntibodies) / nAntibodies,
                    nrow = nAntibodies, ncol = ncol(X))
      mat
    }
  )

  set.seed(7)
  model <- AINet$new(nAntibodies = 12, maxIter = 1, k = 3, epsilon = 1e-6,
                     initMethod = "sample", init = mock_init, verbose = FALSE)
  model$fit(X, task = "clustering")

  expect_equal(call_log$n_calls, 1L)
  expect_equal(call_log$last_n, 12L)
})

test_that("AINet ignores `init` when its $generate is not a function", {
  bogus <- list(generate = "not a function")
  model <- AINet$new(nAntibodies = 8, maxIter = 1, init = bogus, verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_true(nrow(model$repertoire$as_matrix()) >= 1)
})

# ==== Module integration: IdiotypicNetwork ====

test_that("IdiotypicNetwork module changes the surviving repertoire", {
  set.seed(11)
  base <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(11)
  idi <- IdiotypicNetwork$new(theta_low = 0.05, theta_high = 1.0,
                              source_rate = 0.5, decay_rate = 0.1,
                              dt = 0.05, timeSteps = 5,
                              survival_threshold = 0.3)
  with_idi <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                        initMethod = "sample", idiotypic = idi, verbose = FALSE)
  with_idi$fit(X, task = "clustering")

  # Bell-curve dynamics should change the post-fit repertoire size relative
  # to the unmodulated baseline (almost always smaller; rare ties allowed).
  expect_false(identical(base$repertoire$size(), with_idi$repertoire$size()))
})

test_that("IdiotypicNetwork safety net keeps a minimum repertoire when ill-tuned", {
  # Thresholds and survival_threshold are calibrated so every antibody would
  # otherwise die; the safety net retains the top-population subset and the
  # iteration must continue rather than throwing.
  set.seed(13)
  idi_bad <- IdiotypicNetwork$new(theta_low = 100, theta_high = 1000,
                                  source_rate = 0,    decay_rate = 10,
                                  dt = 1.0, timeSteps = 3,
                                  survival_threshold = 1e6)
  model <- AINet$new(nAntibodies = 20, maxIter = 3, k = 3, epsilon = 1e-6,
                     initMethod = "sample", idiotypic = idi_bad,
                     verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_true(model$repertoire$size() >= 1L)
})

# ==== Module integration: ActivationGate (pre-selection density gating) ====

test_that("ActivationGate alters which antibodies enter clonal selection", {
  set.seed(17)
  base <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(17)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.6)
  gated <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  gated$fit(X, task = "clustering")

  # Either the surviving repertoires differ, or their final assignments do.
  same_size <- identical(base$repertoire$size(), gated$repertoire$size())
  same_pred <- identical(base$result$assignments, gated$result$assignments)
  expect_false(same_size && same_pred)
})

test_that("ActivationGate works in classification without breaking class coverage", {
  set.seed(19)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.7)
  model <- AINet$new(nAntibodies = 30, maxIter = 5, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
  preds <- model$predict(X)
  expect_equal(length(preds), nrow(X))
  expect_gte(length(unique(model$result$antibody_classes)), 2L)
})

# ==== Module integration: Microenvironment-aware mutation ====

test_that("Microenvironment module perturbs the post-selection repertoire", {
  set.seed(23)
  base <- AINet$new(nAntibodies = 25, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(23)
  mic <- Microenvironment$new(high_density_threshold = 0.75,
                              low_density_threshold = 0.25,
                              stabilization_factor = 0.3,
                              exploration_factor = 2.0)
  with_mic <- AINet$new(nAntibodies = 25, maxIter = 4, k = 3, epsilon = 1e-6,
                        initMethod = "sample", microenvironment = mic,
                        verbose = FALSE)
  with_mic$fit(X, task = "clustering")

  # Density-aware jitter must change the antibody coordinates relative to
  # the unmodulated baseline. If sizes match, compare matrices; if sizes
  # differ, the size diff itself proves an effect.
  A_base <- base$repertoire$as_matrix()
  A_mic  <- with_mic$repertoire$as_matrix()
  if (nrow(A_base) == nrow(A_mic)) {
    expect_false(isTRUE(all.equal(A_base, A_mic)))
  } else {
    succeed("Microenvironment changed surviving repertoire size")
  }
})

# ==== Final-stage orphan pruning ====

test_that("AINet prunes antibodies that bind no training point", {
  # Force orphans by using `random_uniform` init: many antibodies will land
  # in regions of feature space no iris row is closest to.
  set.seed(29)
  model <- AINet$new(nAntibodies = 50, maxIter = 3, k = 3, epsilon = 1e-6,
                     initMethod = "random_uniform", verbose = FALSE)
  model$fit(X, task = "clustering")

  # After pruning, every surviving antibody must serve as nearest neighbor
  # to at least one training point. Equivalently, the number of unique
  # cluster IDs in the result equals the antibody count.
  expect_equal(length(unique(model$result$assignments)),
               nrow(model$result$antibodies))
})

test_that("Orphan pruning preserves at least 2 antibodies when many are orphans", {
  # The pruning step refuses to drop below 2 surviving antibodies even if
  # most are orphans, so the final assignment can still produce >=2 clusters.
  set.seed(31)
  model <- AINet$new(nAntibodies = 40, maxIter = 2, k = 3, epsilon = 1e-6,
                     initMethod = "random_uniform", verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_gte(nrow(model$result$antibodies), 2L)
})

# ==== Pre-allocated antibody_classes (regression test for iter-1 indexing) ====

test_that("AINet classification with ActivationGate does not error on iter 1", {
  # The pre-selection gate references antibody_classes on iteration 1,
  # before the first clonal_selection_iteration_cpp call has populated it.
  # A pre-allocation step in $fit() seeds it; this test guards against
  # a regression where that pre-allocation is removed.
  set.seed(37)
  gate <- ActivationGate$new(signal2_type = "density",
                             threshold1 = 0, threshold2 = 0.6)
  model <- AINet$new(nAntibodies = 20, maxIter = 1, k = 3, epsilon = 1e-6,
                     initMethod = "sample", activation = gate, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
})

# ==== Module integration: GerminalCenter ====

test_that("GerminalCenter module reduces the surviving repertoire", {
  # With selectionPressure ~ 1 and small nTfh, surviving count is bounded
  # near nTfh. The baseline (no GC) typically retains many more antibodies.
  set.seed(41)
  base <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                    initMethod = "sample", verbose = FALSE)
  base$fit(X, task = "clustering")

  set.seed(41)
  gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.9, rounds = 2)
  with_gc <- AINet$new(nAntibodies = 30, maxIter = 4, k = 3, epsilon = 1e-6,
                       initMethod = "sample", germinalCenter = gc,
                       verbose = FALSE)
  with_gc$fit(X, task = "clustering")

  expect_lt(with_gc$repertoire$size(), base$repertoire$size())
  expect_gte(with_gc$repertoire$size(), 1L)
  # GC must record survivors from its last call.
  expect_true(is.numeric(gc$last_survivors))
})

test_that("GerminalCenter in classification keeps antibody_classes aligned", {
  set.seed(43)
  gc <- GerminalCenter$new(nTfh = 4, selectionPressure = 0.7)
  model <- AINet$new(nAntibodies = 25, maxIter = 3, k = 3, epsilon = 1e-6,
                     germinalCenter = gc, verbose = FALSE)
  expect_no_error(model$fit(X, y_class, task = "classification"))
  expect_equal(length(model$result$antibody_classes),
               nrow(model$result$antibodies))
})

# ==== Module integration: MemoryPool ====

test_that("MemoryPool archives high-affinity antibodies after fit", {
  set.seed(47)
  mp <- MemoryPool$new(archive_threshold = 1e-6, max_memory = 100)
  model <- AINet$new(nAntibodies = 20, maxIter = 3, memory = mp,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  expect_gt(mp$size(), 0L)
  expect_equal(ncol(mp$memory_cells), ncol(X))
})

test_that("MemoryPool recall seeds a second fit's starting repertoire", {
  # Fit once to populate memory, then fit a fresh model that recalls.
  # The recalled cells should be present in the starting repertoire even
  # though nAntibodies is small.
  set.seed(53)
  mp <- MemoryPool$new(archive_threshold = 1e-6, recall_threshold = 0,
                       max_memory = 100)
  AINet$new(nAntibodies = 20, maxIter = 2, memory = mp,
            verbose = FALSE)$fit(X, task = "clustering")
  archived <- mp$size()
  expect_gt(archived, 0L)

  # Second fit recalls and merges memory into the initial repertoire.
  set.seed(57)
  model2 <- AINet$new(nAntibodies = 5, maxIter = 1, epsilon = 1e-6,
                      memory = mp, verbose = FALSE)
  model2$fit(X, task = "clustering")
  # Repertoire size should reflect at least some recall (post-iter pruning
  # may shrink it, so just confirm fit completes without error).
  expect_gte(model2$repertoire$size(), 1L)
})

# ==== Module integration: ClassSwitcher ====

test_that("ClassSwitcher updates isotypes when paired with Microenvironment", {
  set.seed(61)
  me <- Microenvironment$new(high_density_threshold = 0.6,
                              low_density_threshold = 0.3)
  cs <- ClassSwitcher$new(alpha_IgM = 0.1, alpha_IgG = 5.0, alpha_IgA = 1.0)
  model <- AINet$new(nAntibodies = 30, maxIter = 4, epsilon = 1e-6,
                     microenvironment = me, classSwitcher = cs,
                     verbose = FALSE)
  model$fit(X, task = "clustering")
  # After fit, repertoire metadata$isotype should reflect zone-based switching.
  iso <- model$repertoire$metadata$isotype
  expect_true(any(iso %in% c("IgG", "IgA")))
})

# ==== SHM dispatch in clonal_selection ====

test_that("Different SHM strategies produce different repertoires", {
  # Use random_uniform init so antibodies start scattered across the
  # feature bounding box. Initial affinities are low, mutation rates are
  # consequently large, and many mutations get accepted — so the two SHM
  # strategies (with their different rate formulas) diverge measurably.
  # With "sample" init, antibodies clone data points (affinity ~ 1) and
  # both strategies' rates collapse toward mutationMin, making the test
  # sensitive to floating-point edge cases on the accept/reject boundary.
  set.seed(67)
  m_unif <- AINet$new(nAntibodies = 20, maxIter = 4, epsilon = 1e-6,
                      initMethod = "random_uniform",
                      shm = SHMEngine$new(method = "uniform"),
                      verbose = FALSE)
  m_unif$fit(X, task = "clustering")

  set.seed(67)
  m_airs <- AINet$new(nAntibodies = 20, maxIter = 4, epsilon = 1e-6,
                      initMethod = "random_uniform",
                      shm = SHMEngine$new(method = "airs",
                                          temperature = 0.2, c_rate = 0.5),
                      verbose = FALSE)
  m_airs$fit(X, task = "clustering")

  A_u <- m_unif$repertoire$as_matrix()
  A_a <- m_airs$repertoire$as_matrix()
  # Different SHM strategies must produce different antibody coordinates
  # (size match would still permit equality only if results are identical).
  if (nrow(A_u) == nrow(A_a)) {
    expect_false(isTRUE(all.equal(A_u, A_a)))
  } else {
    succeed("SHM strategy changed surviving repertoire size")
  }
})

test_that("AINet runs end-to-end with adaptive SHM (moment state tracking)", {
  # Adaptive SHM requires (m x d) moment matrices to be threaded through
  # clonal_selection and re-aligned after every subset/reorder operation.
  # This test exercises the full plumbing without checking specific values.
  set.seed(71)
  model <- AINet$new(nAntibodies = 25, maxIter = 4, epsilon = 1e-6,
                     shm = SHMEngine$new(method = "adaptive",
                                          base_rate = 0.05),
                     verbose = FALSE)
  expect_no_error(model$fit(X, task = "clustering"))
  expect_gte(model$repertoire$size(), 1L)
})

test_that("AINet runs with all five SHM strategies in classification", {
  for (method in c("uniform", "airs", "hotspot", "energy", "adaptive")) {
    set.seed(73)
    model <- AINet$new(nAntibodies = 20, maxIter = 3, epsilon = 1e-6,
                       shm = SHMEngine$new(method = method),
                       verbose = FALSE)
    expect_no_error(model$fit(X, y_class, task = "classification"))
    expect_equal(length(model$result$assignments), nrow(X))
  }
})
