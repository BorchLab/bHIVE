#' @title ClassSwitcher
#' @description Implements antibody isotype/class switching, allowing antibodies
#' to change their matching breadth. Inspired by real B cell class switching
#' from IgM (broad, pentameric) to IgG (specific, monomeric) to IgA (mucosal).
#'
#' @details
#' In bHIVE, class switching modifies the effective affinity kernel width:
#' \itemize{
#'   \item \code{IgM}: Broad matching (large kernel width) -- good for initial
#'     exploration and capturing general patterns
#'   \item \code{IgG}: Specific matching (small kernel width) -- good for
#'     fine-grained discrimination after patterns are identified
#'   \item \code{IgA}: Boundary patrol (medium kernel width) -- good for
#'     maintaining coverage at decision boundaries
#' }
#'
#' @examples
#' # Switch antibody isotypes based on microenvironment zones
#' A <- matrix(rnorm(50), nrow = 10, ncol = 5)
#' rep <- ImmuneRepertoire$new(A)
#' cs <- ClassSwitcher$new(alpha_IgM = 0.1, alpha_IgG = 5.0)
#' zones <- sample(c("stable", "explore", "boundary"), 10, replace = TRUE)
#' alphas <- cs$switch_isotypes(rep, zones)
#' table(rep$metadata$isotype)  # IgM, IgG, IgA distribution
#'
#' @importFrom R6 R6Class
#' @export
ClassSwitcher <- R6::R6Class(
  "ClassSwitcher",
  public = list(

    #' @field alpha_IgM Kernel width for IgM mode (broad).
    alpha_IgM = NULL,

    #' @field alpha_IgG Kernel width for IgG mode (specific).
    alpha_IgG = NULL,

    #' @field alpha_IgA Kernel width for IgA mode (boundary).
    alpha_IgA = NULL,

    #' @description Create a new ClassSwitcher.
    #' @param alpha_IgM Numeric. Kernel width for broad matching.
    #' @param alpha_IgG Numeric. Kernel width for specific matching.
    #' @param alpha_IgA Numeric. Kernel width for boundary matching.
    initialize = function(alpha_IgM = 0.1, alpha_IgG = 5.0, alpha_IgA = 1.0) {
      self$alpha_IgM <- alpha_IgM
      self$alpha_IgG <- alpha_IgG
      self$alpha_IgA <- alpha_IgA
    },

    #' @description Determine appropriate isotype for each antibody based on
    #'   its microenvironment zone.
    #' @param repertoire An \code{\link{ImmuneRepertoire}}.
    #' @param zones Character vector from Microenvironment assessment.
    #' @return Named numeric vector of alpha values per antibody.
    switch_isotypes = function(repertoire, zones) {
      m <- repertoire$size()
      alphas <- numeric(m)

      for (i in seq_len(m)) {
        current <- repertoire$metadata$isotype[i]
        zone <- zones[i]

        # Switch logic
        new_isotype <- switch(zone,
          "stable"   = "IgG",   # high density -> switch to specific
          "explore"  = "IgM",   # low density -> switch to broad
          "boundary" = "IgA",   # boundary -> intermediate
          current                # default: keep current
        )

        repertoire$metadata$isotype[i] <- new_isotype
        alphas[i] <- switch(new_isotype,
          "IgM" = self$alpha_IgM,
          "IgG" = self$alpha_IgG,
          "IgA" = self$alpha_IgA,
          self$alpha_IgA  # default
        )
      }

      alphas
    },

    #' @description Get alpha value for a given isotype.
    #' @param isotype Character. "IgM", "IgG", or "IgA".
    #' @return Numeric.
    get_alpha = function(isotype) {
      switch(isotype,
        "IgM" = self$alpha_IgM,
        "IgG" = self$alpha_IgG,
        "IgA" = self$alpha_IgA,
        self$alpha_IgA
      )
    },

    #' @description Print summary.
    print = function(...) {
      cat("<ClassSwitcher>\n")
      cat(sprintf("  IgM (broad):    alpha = %.3f\n", self$alpha_IgM))
      cat(sprintf("  IgG (specific): alpha = %.3f\n", self$alpha_IgG))
      cat(sprintf("  IgA (boundary): alpha = %.3f\n", self$alpha_IgA))
      invisible(self)
    }
  )
)
