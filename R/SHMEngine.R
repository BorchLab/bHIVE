#' @title SHMEngine
#' @description Somatic Hypermutation Engine implementing five biologically-
#' inspired mutation strategies for the bHIVE artificial immune system.
#'
#' @details
#' The five strategies are:
#' \describe{
#'   \item{uniform}{Classic random Gaussian noise. Mutation rate = (1-affinity) *
#'     decay^(iter-1). This is the original bHIVE behavior.}
#'   \item{airs}{AIRS-style affinity-proportional mutation. Rate = c * exp(-affinity / T).
#'     From Watkins & Timmis (AIRS2), achieving 50\% better data reduction than uniform.}
#'   \item{hotspot}{Feature-importance-weighted mutation. Features with higher gradient
#'     magnitude mutate more, analogous to AID targeting WRCY motifs in real SHM.}
#'   \item{energy}{Energy-budget-constrained mutation. Total mutation magnitude bounded
#'     by E = E_0 * (1-affinity)^2, inspired by Kleinstein's E_SHM ~ N_Mut^2 model.}
#'   \item{adaptive}{Per-feature adaptive mutation rate with moment tracking, directly
#'     implementing the Reddy (2026) insight that SHM and Adam optimizer are the
#'     same algorithm discovered independently by evolution and machine learning.}
#' }
#'
#' @examples
#' # Create different SHM engines
#' shm_uniform  <- SHMEngine$new(method = "uniform")
#' shm_adaptive <- SHMEngine$new(method = "adaptive", base_rate = 0.1)
#' shm_airs     <- SHMEngine$new(method = "airs", temperature = 0.3)
#' print(shm_adaptive)
#'
#' @param method Character. Mutation strategy.
#' @param decay Numeric. Per-iteration mutation rate decay (uniform method).
#' @param mutationMin Numeric. Minimum mutation rate floor.
#' @param c_rate Numeric. Scaling constant (airs method).
#' @param temperature Numeric. Temperature parameter (airs method).
#' @param E_0 Numeric. Energy budget base (energy method).
#' @param base_rate Numeric. Base mutation rate (hotspot, adaptive methods).
#' @param beta1 Numeric. First moment decay (adaptive method, like Adam).
#' @param beta2 Numeric. Second moment decay (adaptive method, like Adam).
#' @param adam_epsilon Numeric. Numerical stability (adaptive method).
#'
#' @importFrom R6 R6Class
#' @export
SHMEngine <- R6::R6Class(
  "SHMEngine",
  public = list(

    #' @field method Character. One of "uniform", "airs", "hotspot", "energy", "adaptive".
    method = NULL,

    #' @field params Named list of method-specific parameters.
    params = NULL,

    #' @field m1_state First moment state matrix (for adaptive method).
    m1_state = NULL,

    #' @field m2_state Second moment state matrix (for adaptive method).
    m2_state = NULL,

    #' @description Create a new SHMEngine.
    #' @param method Character. Mutation strategy.
    #' @param decay Numeric. Per-iteration mutation rate decay (uniform method).
    #' @param mutationMin Numeric. Minimum mutation rate floor.
    #' @param c_rate Numeric. Scaling constant (airs method).
    #' @param temperature Numeric. Temperature parameter (airs method).
    #' @param E_0 Numeric. Energy budget base (energy method).
    #' @param base_rate Numeric. Base mutation rate (hotspot, adaptive methods).
    #' @param beta1 Numeric. First moment decay (adaptive method, like Adam).
    #' @param beta2 Numeric. Second moment decay (adaptive method, like Adam).
    #' @param adam_epsilon Numeric. Numerical stability (adaptive method).
    initialize = function(method = "uniform",
                          decay = 1.0,
                          mutationMin = 0.01,
                          c_rate = 1.0,
                          temperature = 0.5,
                          E_0 = 1.0,
                          base_rate = 0.1,
                          beta1 = 0.9,
                          beta2 = 0.999,
                          adam_epsilon = 1e-8) {
      self$method <- match.arg(method, c("uniform", "airs", "hotspot",
                                          "energy", "adaptive"))
      stopifnot(
        "decay must be in (0, 1]" = is.numeric(decay) && decay > 0 && decay <= 1,
        "mutationMin must be non-negative" = is.numeric(mutationMin) && mutationMin >= 0,
        "temperature must be positive" = is.numeric(temperature) && temperature > 0,
        "E_0 must be positive" = is.numeric(E_0) && E_0 > 0,
        "base_rate must be positive" = is.numeric(base_rate) && base_rate > 0,
        "beta1 must be in [0, 1)" = is.numeric(beta1) && beta1 >= 0 && beta1 < 1,
        "beta2 must be in [0, 1)" = is.numeric(beta2) && beta2 >= 0 && beta2 < 1,
        "adam_epsilon must be positive" = is.numeric(adam_epsilon) && adam_epsilon > 0
      )
      self$params <- list(
        decay        = decay,
        mutationMin  = mutationMin,
        c_rate       = c_rate,
        temperature  = temperature,
        E_0          = E_0,
        base_rate    = base_rate,
        beta1        = beta1,
        beta2        = beta2,
        adam_epsilon = adam_epsilon
      )
    },

    #' @description Initialize internal state for adaptive method.
    #' @param nAntibodies Integer. Number of antibodies.
    #' @param nFeatures Integer. Number of features.
    init_state = function(nAntibodies, nFeatures) {
      if (self$method == "adaptive") {
        self$m1_state <- matrix(0, nrow = nAntibodies, ncol = nFeatures)
        self$m2_state <- matrix(0, nrow = nAntibodies, ncol = nFeatures)
      }
    },

    #' @description Reset moment states (e.g., after suppression changes antibody count).
    #' @param nAntibodies Integer. New number of antibodies.
    #' @param nFeatures Integer. Number of features.
    #' @param kept_idx Integer vector. Indices of antibodies that were kept.
    reset_state = function(nAntibodies, nFeatures, kept_idx = NULL) {
      if (self$method == "adaptive") {
        if (!is.null(kept_idx) && !is.null(self$m1_state)) {
          self$m1_state <- self$m1_state[kept_idx, , drop = FALSE]
          self$m2_state <- self$m2_state[kept_idx, , drop = FALSE]
        } else {
          self$m1_state <- matrix(0, nrow = nAntibodies, ncol = nFeatures)
          self$m2_state <- matrix(0, nrow = nAntibodies, ncol = nFeatures)
        }
      }
    },

    #' @description Print summary.
    #' @param ... Not used.
    print = function(...) {
      cat(sprintf("<SHMEngine> method='%s'\n", self$method))
      key_params <- switch(self$method,
        "uniform"  = c("decay", "mutationMin"),
        "airs"     = c("c_rate", "temperature"),
        "hotspot"  = c("base_rate"),
        "energy"   = c("E_0"),
        "adaptive" = c("base_rate", "beta1", "beta2")
      )
      for (p in key_params) {
        cat(sprintf("  %s: %s\n", p, self$params[[p]]))
      }
      invisible(self)
    }
  )
)
