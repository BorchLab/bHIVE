#' @title MemoryPool
#' @description Manages long-lived memory cells that can be recalled when
#' distribution shifts are detected. Implements the immunological distinction
#' between short-lived effector cells and long-lived memory cells.
#'
#' @examples
#' # Archive and recall memory cells
#' data(iris)
#' X <- as.matrix(iris[, 1:4])
#' A <- X[sample(150, 10), ]
#' mp <- MemoryPool$new(archive_threshold = 0.01)
#' rep <- ImmuneRepertoire$new(A)
#' mp$archive(rep, X)
#' mp$size()  # number of archived memories
#' recalled <- mp$recall(X[1:5, ])
#' nrow(recalled)  # memories relevant to query
#'
#' @importFrom R6 R6Class
#' @export
MemoryPool <- R6::R6Class(
  "MemoryPool",
  public = list(

    #' @field memory_cells Numeric matrix of archived memory antibodies.
    memory_cells = NULL,

    #' @field memory_metadata Data frame of metadata for memory cells.
    memory_metadata = NULL,

    #' @field archive_threshold Affinity threshold for archiving (only high-quality
    #'   antibodies become memory).
    archive_threshold = NULL,

    #' @field max_memory Maximum number of memory cells to store.
    max_memory = NULL,

    #' @field recall_threshold Threshold for triggering memory recall.
    recall_threshold = NULL,

    #' @description Create a new MemoryPool.
    #' @param archive_threshold Numeric. Minimum average affinity to archive.
    #' @param max_memory Integer. Maximum memory cells.
    #' @param recall_threshold Numeric. Minimum affinity to recall a memory.
    initialize = function(archive_threshold = 0.5,
                          max_memory = 100,
                          recall_threshold = 0.3) {
      self$archive_threshold <- archive_threshold
      self$max_memory        <- max_memory
      self$recall_threshold  <- recall_threshold
      self$memory_cells      <- NULL
      self$memory_metadata   <- NULL
    },

    #' @description Archive high-performing antibodies as memory cells.
    #' @param repertoire An \code{\link{ImmuneRepertoire}}.
    #' @param X Training data matrix.
    #' @param affinityFunc Character. Affinity function.
    #' @param affinityParams List. Affinity parameters.
    #' @return Integer. Number of new memory cells archived.
    archive = function(repertoire, X, affinityFunc = "gaussian",
                       affinityParams = list(alpha = 1, c = 1, p = 2)) {
      A <- repertoire$as_matrix()
      aff <- compute_affinity_matrix(
        X, A, affinityFunc,
        alpha = affinityParams$alpha %||% 1,
        c = affinityParams$c %||% 1,
        p = affinityParams$p %||% 2
      )

      # Average affinity per antibody
      avg_aff <- colMeans(aff)
      candidates <- which(avg_aff >= self$archive_threshold)

      if (length(candidates) == 0) return(0L)

      new_memory <- A[candidates, , drop = FALSE]
      new_meta <- repertoire$metadata[candidates, , drop = FALSE]
      new_meta$state <- "memory"

      if (is.null(self$memory_cells)) {
        self$memory_cells    <- new_memory
        self$memory_metadata <- new_meta
      } else {
        self$memory_cells    <- rbind(self$memory_cells, new_memory)
        self$memory_metadata <- rbind(self$memory_metadata, new_meta)
      }

      # Trim to max_memory (keep most recent)
      n_mem <- nrow(self$memory_cells)
      if (n_mem > self$max_memory) {
        keep <- seq(n_mem - self$max_memory + 1, n_mem)
        self$memory_cells    <- self$memory_cells[keep, , drop = FALSE]
        self$memory_metadata <- self$memory_metadata[keep, , drop = FALSE]
      }

      length(candidates)
    },

    #' @description Recall memory cells relevant to current data.
    #' @param X Data matrix to match against memory.
    #' @param affinityFunc Character. Affinity function.
    #' @param affinityParams List. Affinity parameters.
    #' @return Numeric matrix of recalled memory cells (may be empty).
    recall = function(X, affinityFunc = "gaussian",
                      affinityParams = list(alpha = 1, c = 1, p = 2)) {
      if (is.null(self$memory_cells) || nrow(self$memory_cells) == 0) {
        return(matrix(0, nrow = 0, ncol = ncol(as.matrix(X))))
      }

      aff <- compute_affinity_matrix(
        X, self$memory_cells, affinityFunc,
        alpha = affinityParams$alpha %||% 1,
        c = affinityParams$c %||% 1,
        p = affinityParams$p %||% 2
      )

      # Recall memories that are relevant (high average affinity to current data)
      avg_aff <- colMeans(aff)
      relevant <- which(avg_aff >= self$recall_threshold)

      if (length(relevant) == 0) {
        return(matrix(0, nrow = 0, ncol = ncol(self$memory_cells)))
      }

      self$memory_cells[relevant, , drop = FALSE]
    },

    #' @description Get current memory pool size.
    #' @return Integer.
    size = function() {
      if (is.null(self$memory_cells)) 0L else nrow(self$memory_cells)
    },

    #' @description Print summary.
    print = function(...) {
      cat(sprintf("<MemoryPool> %d cells (max %d)\n", self$size(), self$max_memory))
      cat(sprintf("  Archive threshold: %.3f\n", self$archive_threshold))
      cat(sprintf("  Recall threshold:  %.3f\n", self$recall_threshold))
      invisible(self)
    }
  )
)
