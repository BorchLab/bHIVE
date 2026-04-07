#' @title Microenvironment
#' @description Models local microenvironment cues that influence antibody
#' behavior based on the density and structure of nearby data points.
#'
#' @details
#' In real immunity, B cell fate is strongly influenced by local signals:
#' chemokines, cytokines, and interactions with stromal cells in specific
#' tissue microenvironments. This module translates that concept into
#' density-dependent adaptation:
#'
#' \itemize{
#'   \item \strong{High density zones}: Promote memory formation (stabilize
#'     antibodies, reduce mutation rate)
#'   \item \strong{Low density zones}: Promote exploration (increase mutation
#'     rate for broader search)
#'   \item \strong{Boundary zones}: Trigger class switching (change matching
#'     breadth between IgM-like broad and IgG-like specific modes)
#'   \item \strong{Chemokine-like gradients}: Bias mutation direction toward
#'     higher-density regions
#' }
#'
#' @examples
#' # Assess microenvironment around antibodies
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' me <- Microenvironment$new()
#' rep <- ImmuneRepertoire$new(X[sample(150, 15), ])
#' env <- me$assess(rep, X)
#' table(env$zones)  # stable, explore, boundary
#' env$mutation_modifiers  # per-antibody rate scaling
#'
#' @importFrom R6 R6Class
#' @export
Microenvironment <- R6::R6Class(
  "Microenvironment",
  public = list(

    #' @field density_bandwidth Bandwidth for kernel density estimation.
    density_bandwidth = NULL,

    #' @field high_density_threshold Density percentile above which antibodies stabilize.
    high_density_threshold = NULL,

    #' @field low_density_threshold Density percentile below which antibodies explore.
    low_density_threshold = NULL,

    #' @field stabilization_factor Mutation rate multiplier for high-density zones.
    stabilization_factor = NULL,

    #' @field exploration_factor Mutation rate multiplier for low-density zones.
    exploration_factor = NULL,

    #' @field last_densities Per-antibody local density from last assessment.
    last_densities = NULL,

    #' @field last_zones Per-antibody zone classification from last assessment.
    last_zones = NULL,

    #' @description Create a new Microenvironment module.
    #' @param density_bandwidth Numeric. KDE bandwidth (NULL for auto).
    #' @param high_density_threshold Numeric [0,1]. Percentile threshold for stabilization.
    #' @param low_density_threshold Numeric [0,1]. Percentile threshold for exploration.
    #' @param stabilization_factor Numeric. Mutation rate multiplier for stable zones.
    #' @param exploration_factor Numeric. Mutation rate multiplier for exploration zones.
    initialize = function(density_bandwidth = NULL,
                          high_density_threshold = 0.75,
                          low_density_threshold = 0.25,
                          stabilization_factor = 0.3,
                          exploration_factor = 2.0) {
      self$density_bandwidth       <- density_bandwidth
      self$high_density_threshold  <- high_density_threshold
      self$low_density_threshold   <- low_density_threshold
      self$stabilization_factor    <- stabilization_factor
      self$exploration_factor      <- exploration_factor
    },

    #' @description Assess the microenvironment around each antibody.
    #' @param repertoire An \code{\link{ImmuneRepertoire}} object.
    #' @param X Numeric matrix of training data (n x d).
    #' @param affinityFunc Character. Affinity function.
    #' @param affinityParams List. Affinity parameters.
    #' @return Named list with densities, zones, and mutation_modifiers per antibody.
    assess = function(repertoire, X, affinityFunc = "gaussian",
                      affinityParams = list(alpha = 1, c = 1, p = 2)) {
      A <- repertoire$as_matrix()
      m <- nrow(A)
      n <- nrow(X)

      # Compute affinity to data points
      aff <- compute_affinity_matrix(
        X, A, affinityFunc,
        alpha = affinityParams$alpha %||% 1,
        c = affinityParams$c %||% 1,
        p = affinityParams$p %||% 2
      )

      # Local density: sum of affinities to all data points (KDE-like)
      densities <- colSums(aff)  # m-length vector

      # Compute percentile thresholds
      low_thresh  <- quantile(densities, self$low_density_threshold)
      high_thresh <- quantile(densities, self$high_density_threshold)

      # Classify zones
      zones <- ifelse(densities >= high_thresh, "stable",
                       ifelse(densities <= low_thresh, "explore", "boundary"))

      # Compute mutation rate modifiers
      mutation_modifiers <- rep(1.0, m)
      mutation_modifiers[zones == "stable"]  <- self$stabilization_factor
      mutation_modifiers[zones == "explore"] <- self$exploration_factor
      mutation_modifiers[zones == "boundary"] <- 1.0  # neutral

      # Compute gradient direction (mean direction toward data from each antibody)
      # Weighted by affinity -- chemokine-like attraction
      gradients <- matrix(0, nrow = m, ncol = ncol(A))
      for (j in seq_len(m)) {
        weights <- aff[, j]
        total_w <- sum(weights)
        if (total_w > 0) {
          gradients[j, ] <- colSums(sweep(X, 2, A[j, ], `-`) * weights) / total_w
        }
      }

      # Update metadata states
      repertoire$metadata$state[zones == "stable"]   <- "memory"
      repertoire$metadata$state[zones == "explore"]  <- "activated"
      repertoire$metadata$state[zones == "boundary"] <- "activated"

      self$last_densities <- densities
      self$last_zones     <- zones

      list(
        densities           = densities,
        zones               = zones,
        mutation_modifiers  = mutation_modifiers,
        gradients           = gradients
      )
    },

    #' @description Print summary.
    print = function(...) {
      cat("<Microenvironment>\n")
      cat(sprintf("  Density thresholds: [%.2f, %.2f] percentile\n",
                  self$low_density_threshold, self$high_density_threshold))
      cat(sprintf("  Stabilization: %.2fx  Exploration: %.2fx\n",
                  self$stabilization_factor, self$exploration_factor))
      if (!is.null(self$last_zones)) {
        cat("  Last assessment:", paste(table(self$last_zones), collapse=", "), "\n")
      }
      invisible(self)
    }
  )
)
