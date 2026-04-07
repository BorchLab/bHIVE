#' @title ImmuneRepertoire
#' @description R6 class representing a collection of antibodies (immune cells)
#' with associated metadata. Core data structure for bHIVE algorithms.
#'
#' @details
#' An ImmuneRepertoire holds a matrix of antibody vectors (each row is one
#' antibody in feature space) plus per-antibody metadata (isotype, state, age,
#' lineage). All heavy computation is dispatched to C++ via RcppArmadillo.
#'
#' @examples
#' # Create a repertoire from random antibodies
#' A <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' rep <- ImmuneRepertoire$new(A)
#' print(rep)
#'
#' # Compute affinity to data
#' X <- matrix(rnorm(100), nrow = 20, ncol = 5)
#' aff <- rep$affinity_matrix(X, "gaussian", list(alpha = 1))
#' dim(aff)  # 20 x 10
#'
#' # Network suppression
#' rep$suppress(epsilon = 1.5, method = "euclidean")
#' rep$size()  # fewer antibodies after suppression
#'
#' @param cells Numeric matrix (nAntibodies x nFeatures).
#' @param metadata Optional data frame with columns: isotype, state, age, lineage.
#'
#' @importFrom R6 R6Class
#' @export
ImmuneRepertoire <- R6::R6Class(
  "ImmuneRepertoire",
  public = list(

    #' @field cells Numeric matrix (nAntibodies x nFeatures).
    cells = NULL,

    #' @field metadata Data frame with per-antibody attributes.
    metadata = NULL,

    #' @description Create a new ImmuneRepertoire.
    #' @param cells Numeric matrix (nAntibodies x nFeatures).
    #' @param metadata Optional data frame with columns: isotype, state, age, lineage.
    initialize = function(cells, metadata = NULL) {
      cells <- as.matrix(cells)
      stopifnot(is.numeric(cells), length(dim(cells)) == 2)
      self$cells <- cells

      m <- nrow(cells)
      if (is.null(metadata)) {
        self$metadata <- data.frame(
          isotype = rep("IgM", m),
          state   = rep("naive", m),
          age     = rep(0L, m),
          lineage = rep(NA_character_, m),
          stringsAsFactors = FALSE
        )
      } else {
        stopifnot(is.data.frame(metadata), nrow(metadata) == m)
        self$metadata <- metadata
      }
    },

    #' @description Compute affinity matrix between data X and antibodies.
    #' @param X Numeric matrix (n x d) of data points.
    #' @param method Affinity function: "gaussian", "laplace", "polynomial",
    #'   "cosine", "hamming".
    #' @param params List with alpha, c, p parameters.
    #' @return Numeric matrix (n x m) of affinity values.
    affinity_matrix = function(X, method = "gaussian",
                               params = list(alpha = 1, c = 1, p = 2)) {
      X <- as.matrix(X)
      compute_affinity_matrix(
        X, self$cells, method,
        alpha = params$alpha %||% 1,
        c     = params$c %||% 1,
        p     = params$p %||% 2
      )
    },

    #' @description Compute distance matrix between data X and antibodies.
    #' @param X Numeric matrix (n x d).
    #' @param method Distance function: "euclidean", "manhattan", "minkowski",
    #'   "cosine", "mahalanobis", "hamming".
    #' @param params List with p, Sigma parameters.
    #' @return Numeric matrix (n x m) of distances.
    distance_matrix = function(X, method = "euclidean",
                               params = list(p = 2, Sigma = NULL)) {
      X <- as.matrix(X)
      Sigma_inv <- if (!is.null(params$Sigma)) solve(params$Sigma) else matrix(0, 0, 0)
      compute_distance_matrix(
        X, self$cells, method,
        p = params$p %||% 2,
        Sigma_inv = Sigma_inv
      )
    },

    #' @description Network suppression: remove redundant antibodies.
    #' @param epsilon Distance threshold for suppression.
    #' @param method Distance function for suppression.
    #' @param params List with p, Sigma parameters.
    #' @return Invisible self (modified in place).
    suppress = function(epsilon, method = "euclidean",
                        params = list(p = 2, Sigma = NULL)) {
      Sigma_inv <- if (!is.null(params$Sigma)) solve(params$Sigma) else matrix(0, 0, 0)
      keep <- network_suppression_cpp(
        self$cells, method, epsilon,
        p = params$p %||% 2,
        Sigma_inv = Sigma_inv
      )
      kept_idx <- which(keep)
      self$cells <- self$cells[kept_idx, , drop = FALSE]
      self$metadata <- self$metadata[kept_idx, , drop = FALSE]
      invisible(self)
    },

    #' @description Get number of antibodies.
    #' @return Integer.
    size = function() nrow(self$cells),

    #' @description Get number of features.
    #' @return Integer.
    n_features = function() ncol(self$cells),

    #' @description Subset the repertoire.
    #' @param idx Integer vector of row indices to keep.
    #' @return Invisible self (modified in place).
    subset = function(idx) {
      self$cells <- self$cells[idx, , drop = FALSE]
      self$metadata <- self$metadata[idx, , drop = FALSE]
      invisible(self)
    },

    #' @description Add antibodies to the repertoire.
    #' @param new_cells Numeric matrix (k x d) of new antibodies.
    #' @param new_metadata Optional data frame of metadata for new antibodies.
    #' @return Invisible self (modified in place).
    add = function(new_cells, new_metadata = NULL) {
      new_cells <- as.matrix(new_cells)
      stopifnot(ncol(new_cells) == self$n_features())
      k <- nrow(new_cells)
      if (is.null(new_metadata)) {
        new_metadata <- data.frame(
          isotype = rep("IgM", k),
          state   = rep("naive", k),
          age     = rep(0L, k),
          lineage = rep(NA_character_, k),
          stringsAsFactors = FALSE
        )
      }
      self$cells <- rbind(self$cells, new_cells)
      self$metadata <- rbind(self$metadata, new_metadata)
      invisible(self)
    },

    #' @description Increment age of all antibodies.
    #' @return Invisible self (modified in place).
    age_all = function() {
      self$metadata$age <- self$metadata$age + 1L
      invisible(self)
    },

    #' @description Convert to plain matrix.
    #' @return Numeric matrix (nAntibodies x nFeatures).
    as_matrix = function() self$cells,

    #' @description Print summary.
    #' @param ... Not used.
    print = function(...) {
      cat(sprintf("<ImmuneRepertoire> %d antibodies x %d features\n",
                  self$size(), self$n_features()))
      if (self$size() > 0) {
        cat("  Isotypes:", paste(table(self$metadata$isotype), collapse=", "), "\n")
        cat("  States:  ", paste(table(self$metadata$state), collapse=", "), "\n")
      }
      invisible(self)
    }
  ),

  private = list(
  )
)
