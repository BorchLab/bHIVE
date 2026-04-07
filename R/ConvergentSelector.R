#' @title ConvergentSelector
#' @description Identifies "public antibodies" shared across independent
#' repertoires, implementing the concept of convergent selection from TCR/BCR
#' immunology as a biologically-motivated ensemble method.
#'
#' @details
#' In real immunity, certain immune receptor sequences appear across multiple
#' individuals (public clones), suggesting they are driven by common selection
#' pressures. Similarly, antibodies that appear across multiple independent
#' bHIVE runs represent the most robust patterns in the data.
#'
#' @examples
#' # Find public antibodies across multiple runs
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' results <- lapply(1:3, function(i) {
#'   m <- AINet$new(nAntibodies = 15, maxIter = 5, verbose = FALSE)
#'   m$fit(X, task = "clustering")
#'   m$result
#' })
#' conv <- ConvergentSelector$new(tolerance = 1.0, min_appearances = 2)
#' public <- conv$from_results(results)
#' nrow(public)  # consensus antibodies
#'
#' @param tolerance Numeric. Maximum distance for two antibodies to be
#'   considered the same across repertoires.
#' @param min_appearances Integer. Minimum repertoires for an antibody to be public.
#'
#' @importFrom R6 R6Class
#' @export
ConvergentSelector <- R6::R6Class(
  "ConvergentSelector",
  public = list(

    #' @field tolerance Distance tolerance for matching antibodies across repertoires.
    tolerance = NULL,

    #' @field min_appearances Minimum number of repertoires an antibody must appear
    #'   in to be considered "public".
    min_appearances = NULL,

    #' @field public_antibodies The identified public antibodies.
    public_antibodies = NULL,

    #' @description Create a new ConvergentSelector.
    #' @param tolerance Numeric. Maximum distance for two antibodies to be
    #'   considered the same across repertoires.
    #' @param min_appearances Integer. Minimum repertoires for an antibody to be public.
    initialize = function(tolerance = 0.5, min_appearances = 2) {
      self$tolerance       <- tolerance
      self$min_appearances <- min_appearances
    },

    #' @description Find public antibodies shared across multiple repertoires.
    #' @param repertoires List of ImmuneRepertoire objects or list of antibody matrices.
    #' @param distFunc Character. Distance function for matching.
    #' @return Numeric matrix of public (consensus) antibodies.
    find_public = function(repertoires, distFunc = "euclidean") {
      # Convert to list of matrices
      matrices <- lapply(repertoires, function(r) {
        if (inherits(r, "ImmuneRepertoire")) r$as_matrix()
        else if (is.list(r) && !is.null(r$antibodies)) r$antibodies
        else as.matrix(r)
      })

      n_reps <- length(matrices)
      if (n_reps < 2) {
        warning("Convergent selection requires at least 2 repertoires.")
        self$public_antibodies <- matrices[[1]]
        return(self$public_antibodies)
      }

      # Use first repertoire as reference
      ref <- matrices[[1]]
      n_ref <- nrow(ref)

      # Count appearances: for each reference antibody, how many other
      # repertoires have a matching antibody within tolerance?
      appearances <- rep(1L, n_ref)  # 1 for appearing in the reference itself

      for (r in 2:n_reps) {
        other <- matrices[[r]]
        # Distance from reference to other
        D <- compute_distance_matrix(ref, other, distFunc, 2.0, matrix(0, 0, 0))
        # For each reference antibody, check if any other antibody is within tolerance
        min_dist <- apply(D, 1, min)
        appearances <- appearances + as.integer(min_dist <= self$tolerance)
      }

      # Keep antibodies appearing in >= min_appearances repertoires
      public_idx <- which(appearances >= self$min_appearances)

      if (length(public_idx) == 0) {
        # Fall back: return antibodies from first repertoire
        warning("No public antibodies found. Try increasing tolerance or ",
                "decreasing min_appearances.")
        self$public_antibodies <- ref
      } else {
        self$public_antibodies <- ref[public_idx, , drop = FALSE]
      }

      self$public_antibodies
    },

    #' @description Run convergent selection from multiple bHIVE results.
    #' @param results List of bHIVE result objects (each with $antibodies).
    #' @param distFunc Character. Distance function.
    #' @return Numeric matrix of consensus antibodies.
    from_results = function(results, distFunc = "euclidean") {
      self$find_public(results, distFunc)
    },

    #' @description Print summary.
    #' @param ... Not used.
    print = function(...) {
      cat("<ConvergentSelector>\n")
      cat(sprintf("  Tolerance: %.3f\n", self$tolerance))
      cat(sprintf("  Min appearances: %d\n", self$min_appearances))
      if (!is.null(self$public_antibodies)) {
        cat(sprintf("  Public antibodies: %d\n", nrow(self$public_antibodies)))
      }
      invisible(self)
    }
  )
)
