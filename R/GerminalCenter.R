#' @title GerminalCenter
#' @description Models T follicular helper (Tfh) cell selection pressure on
#' B cells within a germinal center reaction. Implements resource competition
#' where antibodies compete for Tfh help, and only helped antibodies survive.
#'
#' @details
#' The germinal center is where B cells undergo affinity maturation through
#' iterative cycles of mutation and selection. Tfh cells act as quality-control
#' selectors:
#'
#' \itemize{
#'   \item Each Tfh evaluates B cell (antibody) quality using a task-aware metric
#'   \item B cells compete for Tfh help (resource competition)
#'   \item Only helped B cells survive to the next round
#'   \item Selection pressure controls the stringency of the process
#' }
#'
#' @examples
#' # Germinal center selection on Iris
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' gc <- GerminalCenter$new(nTfh = 5, selectionPressure = 0.5)
#' rep <- ImmuneRepertoire$new(X[sample(150, 20), ])
#' gc$select(rep, X, iris$Species, "classification")
#' rep$size()  # fewer antibodies after selection
#'
#' @param nTfh Integer. Number of Tfh helper cells. Each helps one B cell.
#' @param selectionPressure Numeric [0,1]. Stringency of selection.
#' @param rounds Integer. Number of competition rounds.
#'
#' @importFrom R6 R6Class
#' @export
GerminalCenter <- R6::R6Class(
  "GerminalCenter",
  public = list(

    #' @field nTfh Number of Tfh selectors (determines how many antibodies survive).
    nTfh = NULL,

    #' @field selectionPressure Numeric [0,1]. 0 = no selection (all survive),
    #'   1 = only the very best survive.
    selectionPressure = NULL,

    #' @field rounds Number of selection rounds per call.
    rounds = NULL,

    #' @field last_survivors Integer vector of survivor indices (relative
    #'   to the input repertoire) from the most recent call to \code{select()}.
    last_survivors = NULL,

    #' @description Create a new GerminalCenter.
    #' @param nTfh Integer. Number of Tfh helper cells. Each helps one B cell.
    #' @param selectionPressure Numeric [0,1]. Stringency of selection.
    #' @param rounds Integer. Number of competition rounds.
    initialize = function(nTfh = 10, selectionPressure = 0.5, rounds = 1) {
      stopifnot(
        "nTfh must be a positive integer" = is.numeric(nTfh) && nTfh >= 1,
        "selectionPressure must be in [0, 1]" = is.numeric(selectionPressure) && selectionPressure >= 0 && selectionPressure <= 1,
        "rounds must be a positive integer" = is.numeric(rounds) && rounds >= 1
      )
      self$nTfh              <- as.integer(nTfh)
      self$selectionPressure <- selectionPressure
      self$rounds            <- as.integer(rounds)
    },

    #' @description Run germinal center selection on a repertoire.
    #' @param repertoire An \code{\link{ImmuneRepertoire}} object.
    #' @param X Numeric matrix of training data.
    #' @param y Factor target vector or NULL for clustering.
    #' @param task Character: "clustering" or "classification".
    #' @param affinityFunc Character. Affinity function for evaluation.
    #' @param affinityParams List. Parameters for affinity function.
    #' @return Integer vector of survivor indices relative to the input
    #'   repertoire (composed across all selection rounds). Also stored on
    #'   \code{self$last_survivors} for inspection. Repertoire is modified
    #'   in place.
    select = function(repertoire, X, y = NULL, task = "clustering",
                      affinityFunc = "gaussian",
                      affinityParams = list(alpha = 1, c = 1, p = 2)) {

      # Track survivor indices composed across rounds so callers can
      # mirror the subset onto external per-antibody state (e.g. class
      # labels, SHM moment matrices).
      current <- seq_len(repertoire$size())

      for (round in seq_len(self$rounds)) {
        m <- repertoire$size()
        if (m <= self$nTfh) break  # Nothing to select

        A <- repertoire$as_matrix()

        # Compute quality score for each antibody
        scores <- private$.compute_quality(A, X, y, task, affinityFunc, affinityParams)

        # Selection: top-scoring antibodies get Tfh help
        # The number surviving depends on selectionPressure
        n_survive <- max(1, round(m * (1 - self$selectionPressure) + self$nTfh * self$selectionPressure))
        n_survive <- min(n_survive, m)

        # Probabilistic selection weighted by scores
        scores_pos <- scores - min(scores) + 1e-10
        probs <- scores_pos / sum(scores_pos)
        survived <- sort(unique(sample.int(m, size = n_survive, replace = FALSE, prob = probs)))

        repertoire$subset(survived)
        current <- current[survived]
      }

      self$last_survivors <- current
      invisible(current)
    },

    #' @description Print summary.
    #' @param ... Not used.
    print = function(...) {
      cat("<GerminalCenter>\n")
      cat(sprintf("  nTfh: %d\n", self$nTfh))
      cat(sprintf("  Selection pressure: %.2f\n", self$selectionPressure))
      cat(sprintf("  Rounds: %d\n", self$rounds))
      invisible(self)
    }
  ),

  private = list(

    .compute_quality = function(A, X, y, task, affinityFunc, affinityParams) {
      m <- nrow(A)

      # Compute affinity matrix (n x m)
      aff <- compute_affinity_matrix(
        X, A, affinityFunc,
        alpha = affinityParams$alpha %||% 1,
        c = affinityParams$c %||% 1,
        p = affinityParams$p %||% 2
      )

      if (task == "clustering") {
        # Quality = average affinity to assigned data points
        # Each data point assigned to nearest antibody
        assignments <- apply(aff, 1, which.max)
        scores <- numeric(m)
        for (j in seq_len(m)) {
          assigned <- which(assignments == j)
          scores[j] <- if (length(assigned) > 0) mean(aff[assigned, j]) else 0
        }

      } else {
        # Classification: quality = classification accuracy of assigned points
        assignments <- apply(aff, 1, which.max)
        scores <- numeric(m)
        # Each antibody gets the most common class label
        for (j in seq_len(m)) {
          assigned <- which(assignments == j)
          if (length(assigned) > 0) {
            labels <- y[assigned]
            # Score = proportion correctly matching the majority class
            scores[j] <- max(table(labels)) / length(labels)
          }
        }
      }

      scores
    }
  )
)
