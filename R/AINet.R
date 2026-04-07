#' @title AINet
#' @description R6 implementation of the Artificial Immune Network algorithm.
#' This is the core bHIVE algorithm using C++ backends for performance-critical
#' operations. Supports composable modules for somatic hypermutation, idiotypic
#' network regulation, germinal center selection, and more.
#'
#' @examples
#' # Clustering with Iris data
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' model <- AINet$new(nAntibodies = 15, maxIter = 10, verbose = FALSE)
#' model$fit(X, task = "clustering")
#' table(model$result$assignments)
#'
#' # Classification
#' model2 <- AINet$new(nAntibodies = 20, maxIter = 10, verbose = FALSE)
#' model2$fit(X, iris$Species, task = "classification")
#' mean(model2$result$assignments == as.character(iris$Species))
#'
#' # Predict on new data
#' preds <- model2$predict(X[1:10, ])
#'
#' @importFrom R6 R6Class
#' @importFrom stats rnorm runif sd
#' @export
AINet <- R6::R6Class(
  "AINet",
  inherit = ImmuneAlgorithm,
  public = list(

    #' @description Create a new AINet algorithm instance.
    #' @param nAntibodies Integer. Initial antibody population size.
    #' @param beta Numeric. Clone multiplier.
    #' @param epsilon Numeric. Suppression distance threshold.
    #' @param maxIter Integer. Maximum iterations.
    #' @param k Integer. Top-k antibodies to clone per data point.
    #' @param affinityFunc Character. Affinity function name.
    #' @param distFunc Character. Distance function name.
    #' @param affinityParams List. Parameters for affinity/distance functions.
    #' @param mutationDecay Numeric. Per-iteration mutation rate decay.
    #' @param mutationMin Numeric. Minimum mutation rate.
    #' @param maxClones Numeric. Maximum clones per antibody.
    #' @param stopTolerance Numeric. Early stopping tolerance.
    #' @param noImprovementLimit Integer. Early stopping patience.
    #' @param initMethod Character. Initialization method.
    #' @param shm An SHMEngine instance or NULL for default uniform mutation.
    #' @param init A VDJLibrary instance or NULL for default initialization.
    #' @param activation An ActivationGate instance or NULL.
    #' @param idiotypic An IdiotypicNetwork instance or NULL.
    #' @param germinalCenter A GerminalCenter instance or NULL.
    #' @param microenvironment A Microenvironment instance or NULL.
    #' @param memory A MemoryPool instance or NULL.
    #' @param verbose Logical. Print progress.
    initialize = function(nAntibodies = 20,
                          beta = 5,
                          epsilon = 0.01,
                          maxIter = 50,
                          k = 3,
                          affinityFunc = "gaussian",
                          distFunc = "euclidean",
                          affinityParams = list(alpha = 1, c = 1, p = 2,
                                                Sigma = NULL),
                          mutationDecay = 1.0,
                          mutationMin = 0.01,
                          maxClones = Inf,
                          stopTolerance = 0.0,
                          noImprovementLimit = Inf,
                          initMethod = "sample",
                          shm = NULL,
                          init = NULL,
                          activation = NULL,
                          idiotypic = NULL,
                          germinalCenter = NULL,
                          microenvironment = NULL,
                          memory = NULL,
                          verbose = TRUE) {

      # Validate numeric parameters
      stopifnot(
        "nAntibodies must be a positive integer" = is.numeric(nAntibodies) && nAntibodies >= 1,
        "beta must be positive" = is.numeric(beta) && beta > 0,
        "epsilon must be non-negative" = is.numeric(epsilon) && epsilon >= 0,
        "maxIter must be a positive integer" = is.numeric(maxIter) && maxIter >= 1,
        "k must be a positive integer" = is.numeric(k) && k >= 1,
        "mutationDecay must be in (0, 1]" = is.numeric(mutationDecay) && mutationDecay > 0 && mutationDecay <= 1,
        "mutationMin must be non-negative" = is.numeric(mutationMin) && mutationMin >= 0
      )

      config <- list(
        nAntibodies       = as.integer(nAntibodies),
        beta              = beta,
        epsilon           = epsilon,
        maxIter           = maxIter,
        k                 = k,
        affinityFunc      = match.arg(affinityFunc, c("gaussian", "laplace",
                                                       "polynomial", "cosine", "hamming")),
        distFunc          = match.arg(distFunc, c("euclidean", "manhattan",
                                                   "minkowski", "cosine",
                                                   "mahalanobis", "hamming")),
        affinityParams    = affinityParams,
        mutationDecay     = mutationDecay,
        mutationMin       = mutationMin,
        maxClones         = maxClones,
        stopTolerance     = stopTolerance,
        noImprovementLimit = noImprovementLimit,
        initMethod        = match.arg(initMethod, c("sample", "random",
                                                     "random_uniform", "kmeans++")),
        verbose           = verbose
      )

      modules <- list(
        shm              = shm,
        init             = init,
        activation       = activation,
        idiotypic        = idiotypic,
        germinalCenter   = germinalCenter,
        microenvironment = microenvironment,
        memory           = memory
      )
      # Remove NULL modules
      modules <- modules[!vapply(modules, is.null, logical(1))]

      super$initialize(config = config, modules = modules)
    },

    #' @description Fit the AINet algorithm to data.
    #' @param X Numeric matrix or data frame (n x d).
    #' @param y Optional target: factor (classification) or numeric (regression).
    #' @param task Character: "clustering", "classification", or "regression".
    #'   Inferred from y if NULL.
    #' @param ... Additional arguments (currently unused).
    #' @return Invisible self, with \code{result} populated.
    fit = function(X, y = NULL, task = NULL, ...) {
      # Validate inputs
      .validate_bHIVE_input(X, y)

      # Infer task
      if (is.null(task)) {
        task <- if (is.null(y)) "clustering"
                else if (is.factor(y)) "classification"
                else "regression"
      }
      task <- match.arg(task, c("clustering", "classification", "regression"))

      X <- as.matrix(X)
      n <- nrow(X)
      d <- ncol(X)
      cfg <- self$config

      # Standardize regression target
      y_orig <- y
      y_mean <- 0; y_sd <- 1
      if (task == "regression") {
        y_mean <- mean(y, na.rm = TRUE)
        y_sd   <- sd(y, na.rm = TRUE)
        if (y_sd == 0) y_sd <- 1
        y <- (y - y_mean) / y_sd
      }

      # ================================
      # 1. Initialize antibody population
      # ================================
      A <- private$.initialize_antibodies(X, cfg$nAntibodies, cfg$initMethod)
      self$repertoire <- ImmuneRepertoire$new(A)
      m <- self$repertoire$size()

      # ================================
      # 2. Task-specific setup
      # ================================
      task_int <- switch(task, clustering = 0L, classification = 1L, regression = 2L)
      nClasses <- 0L
      classes  <- NULL
      if (task == "classification") {
        classes  <- levels(y)
        nClasses <- length(classes)
        y_num    <- as.numeric(y) - 1  # 0-indexed class
      } else if (task == "regression") {
        y_num <- y
      } else {
        y_num <- rep(0, n)
      }

      # Affinity/distance params
      alpha <- cfg$affinityParams$alpha %||% 1
      c_p   <- cfg$affinityParams$c %||% 1
      p_p   <- cfg$affinityParams$p %||% 2
      Sigma_inv <- if (!is.null(cfg$affinityParams$Sigma)) {
        solve(cfg$affinityParams$Sigma)
      } else {
        matrix(0, 0, 0)
      }

      # Early stopping state
      noImproveCount <- 0
      prevCount <- m

      # ================================
      # 3. Main iteration loop
      # ================================
      for (iter in seq_len(cfg$maxIter)) {
        A_current <- self$repertoire$as_matrix()
        m <- nrow(A_current)

        # (a) Clonal selection + mutation [C++]
        cs_result <- clonal_selection_iteration_cpp(
          A_current, X, y_num, task_int, cfg$k, cfg$beta,
          cfg$maxClones, cfg$mutationDecay, cfg$mutationMin,
          iter, cfg$affinityFunc, alpha, c_p, p_p, nClasses
        )
        self$repertoire$cells <- cs_result$A

        # Update labels
        if (task == "classification") {
          antibody_classes <- apply(cs_result$class_counts, 1, function(row) {
            if (all(row == 0)) classes[sample(nClasses, 1)]
            else classes[which.max(row)]
          })
        } else if (task == "regression") {
          antibody_values <- ifelse(cs_result$sum_aff > 0,
                                    cs_result$sum_y / cs_result$sum_aff,
                                    mean(y, na.rm = TRUE))
        }

        # (e) Network suppression [C++]
        # TODO: Replace with idiotypic network dynamics when module is available
        keep <- network_suppression_cpp(
          self$repertoire$cells, cfg$distFunc, cfg$epsilon,
          p_p, Sigma_inv
        )
        kept_idx <- which(keep)
        self$repertoire$subset(kept_idx)
        m_new <- self$repertoire$size()

        if (task == "classification") {
          antibody_classes <- antibody_classes[kept_idx]
        } else if (task == "regression") {
          antibody_values <- antibody_values[kept_idx]
        }

        if (m_new == 0) {
          stop("All antibodies were suppressed. Increase nAntibodies or decrease epsilon.")
        }

        # Record iteration history
        self$history[[iter]] <- list(n_antibodies = m_new)

        # Early stopping
        changeCount <- abs(m_new - prevCount)
        if (changeCount <= cfg$stopTolerance) {
          noImproveCount <- noImproveCount + 1
        } else {
          noImproveCount <- 0
        }
        prevCount <- m_new

        if (noImproveCount >= cfg$noImprovementLimit) {
          if (cfg$verbose) {
            cat("Early stopping: no improvement for", noImproveCount, "iterations.\n")
          }
          break
        }

        if (cfg$verbose) {
          cat(sprintf("Iteration %d | #Antibodies: %d | noImproveCount: %d\n",
                      iter, m_new, noImproveCount))
        }
      }

      # ================================
      # 4. Final assignment [C++]
      # ================================
      A_final <- self$repertoire$as_matrix()
      m <- nrow(A_final)

      if (task == "clustering") {
        fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                   0L, alpha, c_p, p_p, Sigma_inv,
                                   numeric(0), 0.0)
        assignments <- as.numeric(factor(fa$assignments))
        self$result <- list(
          antibodies  = A_final,
          assignments = assignments,
          task        = task
        )
      } else if (task == "classification") {
        fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                   1L, alpha, c_p, p_p, Sigma_inv,
                                   numeric(0), 0.0)
        assignments <- antibody_classes[fa$best_antibody_idx]
        self$result <- list(
          antibodies       = A_final,
          assignments      = assignments,
          antibody_classes = antibody_classes,
          task             = task
        )
      } else {
        fa <- final_assignment_cpp(X, A_final, cfg$affinityFunc, cfg$distFunc,
                                   2L, alpha, c_p, p_p, Sigma_inv,
                                   antibody_values, mean(y, na.rm = TRUE))
        predictions <- fa$predictions * y_sd + y_mean
        self$result <- list(
          antibodies      = A_final,
          assignments     = fa$cluster_assign,
          predictions     = predictions,
          antibody_values = antibody_values,
          overall_mean    = mean(y, na.rm = TRUE),
          task            = task
        )
      }

      invisible(self)
    }
  ),

  private = list(

    .initialize_antibodies = function(X, nAntibodies, method) {
      n <- nrow(X)
      d <- ncol(X)
      switch(
        method,
        "sample" = X[sample.int(n, size = nAntibodies, replace = TRUE), , drop = FALSE],
        "random" = {
          xMean <- colMeans(X)
          xSd   <- apply(X, 2, sd) + 1e-8
          mat   <- matrix(rnorm(nAntibodies * d), nrow = nAntibodies)
          mat   <- sweep(mat, 2, xSd, `*`)
          sweep(mat, 2, xMean, `+`)
        },
        "random_uniform" = {
          xMin <- apply(X, 2, min)
          xMax <- apply(X, 2, max)
          mat  <- matrix(runif(nAntibodies * d), nrow = nAntibodies)
          sweep(sweep(mat, 2, xMax - xMin, `*`), 2, xMin, `+`)
        },
        "kmeans++" = init_kmeanspp_cpp(X, nAntibodies)
      )
    }
  )
)
