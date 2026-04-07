#' @title IdiotypicNetwork
#' @description Implements Jerne's idiotypic network theory for antibody
#' repertoire regulation. Replaces crude epsilon-threshold suppression with
#' principled network dynamics based on Varela & Coutinho's (1991) second-
#' generation immune network model.
#'
#' @details
#' The idiotypic network models antibody-antibody interactions where each
#' antibody's variable region can be recognized by other antibodies. This
#' creates a regulatory network with emergent properties:
#'
#' - A bell-shaped (double-threshold) activation function: too little
#'   stimulation leads to cell death, moderate stimulation to activation,
#'   and excessive stimulation to suppression.
#' - Population dynamics with source, decay, activation, and suppression terms.
#' - Self-organized repertoire structure with memory and tolerance properties.
#'
#' This is the single most novel contribution of the overhauled bHIVE package.
#' No existing AIS implementation uses idiotypic network dynamics for
#' repertoire regulation.
#'
#' @examples
#' # Create and run idiotypic regulation
#' idi <- IdiotypicNetwork$new(theta_low = 0.01, theta_high = 0.5)
#' A <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' rep <- ImmuneRepertoire$new(A)
#' idi$regulate(rep, "gaussian", list(alpha = 0.5))
#' print(idi)
#'
#' @param theta_low Lower activation threshold.
#' @param theta_high Upper activation threshold.
#' @param source_rate Basal cell production rate.
#' @param decay_rate Natural decay rate.
#' @param dt Euler integration time step.
#' @param timeSteps Number of dynamics simulation steps.
#' @param survival_threshold Minimum population to survive.
#'
#' @importFrom R6 R6Class
#' @export
IdiotypicNetwork <- R6::R6Class(
  "IdiotypicNetwork",
  public = list(

    #' @field theta_low Lower activation threshold. Below this, cells die.
    theta_low = NULL,

    #' @field theta_high Upper activation threshold. Above this, cells are suppressed.
    theta_high = NULL,

    #' @field source_rate Rate of new cell generation (basal production).
    source_rate = NULL,

    #' @field decay_rate Natural cell death rate.
    decay_rate = NULL,

    #' @field dt Time step for Euler integration.
    dt = NULL,

    #' @field timeSteps Number of simulation time steps.
    timeSteps = NULL,

    #' @field survival_threshold Minimum population level to survive.
    survival_threshold = NULL,

    #' @field last_dynamics Result from the last regulation step.
    last_dynamics = NULL,

    #' @description Create a new IdiotypicNetwork regulator.
    #' @param theta_low Lower activation threshold.
    #' @param theta_high Upper activation threshold.
    #' @param source_rate Basal cell production rate.
    #' @param decay_rate Natural decay rate.
    #' @param dt Euler integration time step.
    #' @param timeSteps Number of dynamics simulation steps.
    #' @param survival_threshold Minimum population to survive.
    initialize = function(theta_low = 0.01,
                          theta_high = 0.5,
                          source_rate = 0.5,
                          decay_rate = 0.1,
                          dt = 0.1,
                          timeSteps = 20,
                          survival_threshold = 0.5) {
      stopifnot(
        "theta_low must be non-negative" = is.numeric(theta_low) && theta_low >= 0,
        "theta_high must be > theta_low" = is.numeric(theta_high) && theta_high > theta_low,
        "source_rate must be non-negative" = is.numeric(source_rate) && source_rate >= 0,
        "decay_rate must be non-negative" = is.numeric(decay_rate) && decay_rate >= 0,
        "dt must be positive" = is.numeric(dt) && dt > 0,
        "timeSteps must be a positive integer" = is.numeric(timeSteps) && timeSteps >= 1,
        "survival_threshold must be non-negative" = is.numeric(survival_threshold) && survival_threshold >= 0
      )
      self$theta_low           <- theta_low
      self$theta_high          <- theta_high
      self$source_rate         <- source_rate
      self$decay_rate          <- decay_rate
      self$dt                  <- dt
      self$timeSteps           <- as.integer(timeSteps)
      self$survival_threshold  <- survival_threshold
    },

    #' @description Run idiotypic network dynamics on an antibody repertoire.
    #' @param repertoire An \code{\link{ImmuneRepertoire}} object.
    #' @param affinityFunc Character. Affinity function for Ab-Ab interactions.
    #' @param affinityParams List. Parameters for the affinity function.
    #' @return Invisible self. The repertoire is modified in place (dead
    #'   antibodies removed). Access \code{$last_dynamics} for full results.
    regulate = function(repertoire, affinityFunc = "gaussian",
                        affinityParams = list(alpha = 1, c = 1, p = 2)) {
      A <- repertoire$as_matrix()

      result <- idiotypic_dynamics_cpp(
        A,
        affinityFunc,
        affinityParams$alpha %||% 1,
        affinityParams$c %||% 1,
        affinityParams$p %||% 2,
        self$theta_low,
        self$theta_high,
        self$source_rate,
        self$decay_rate,
        self$dt,
        self$timeSteps,
        self$survival_threshold
      )

      self$last_dynamics <- result

      # Apply regulation: remove dead antibodies
      kept_idx <- which(result$keep)
      if (length(kept_idx) == 0) {
        warning("Idiotypic regulation removed all antibodies. ",
                "Consider adjusting thresholds (theta_low, theta_high).")
      } else {
        repertoire$subset(kept_idx)
      }

      invisible(self)
    },

    #' @description Get the Ab-Ab affinity matrix from the last regulation.
    #' @return Numeric matrix (m x m) or NULL if not yet run.
    get_network = function() {
      if (is.null(self$last_dynamics)) return(NULL)
      self$last_dynamics$Ab_Ab_affinity
    },

    #' @description Get population levels from the last regulation.
    #' @return Numeric vector or NULL if not yet run.
    get_population = function() {
      if (is.null(self$last_dynamics)) return(NULL)
      self$last_dynamics$population
    },

    #' @description Print summary.
    #' @param ... Not used.
    print = function(...) {
      cat("<IdiotypicNetwork>\n")
      cat(sprintf("  Activation window: [%.3f, %.3f]\n",
                  self$theta_low, self$theta_high))
      cat(sprintf("  Dynamics: %d steps (dt=%.3f)\n",
                  self$timeSteps, self$dt))
      cat(sprintf("  Source: %.3f  Decay: %.3f  Survival: %.3f\n",
                  self$source_rate, self$decay_rate, self$survival_threshold))
      if (!is.null(self$last_dynamics)) {
        pop <- self$last_dynamics$population
        kept <- sum(self$last_dynamics$keep)
        cat(sprintf("  Last run: %d/%d antibodies survived\n",
                    kept, length(pop)))
      }
      invisible(self)
    }
  )
)
