# ============================================================================
# Matrix implementation (defined first: the ANY method points at it)
# ============================================================================

#' Matrix implementation of bHIVE()
#'
#' bHIVE is a thin functional wrapper over the AINet R6 engine. This is the
#' single code path: the C++ clonal-selection and suppression backends, Lloyd
#' consolidation, scaling, target-K and the composable immunology modules all
#' live in AINet, so swarmbHIVE(), honeycombHIVE() and bHIVEParam inherit
#' every one of them. Internal callers use this function directly rather than
#' the S4 generic.
#' @keywords internal
#' @noRd
.bHIVE_matrix <- function(X, y = NULL, task = NULL, nAntibodies = 20,
                          beta = 5, epsilon = 0.01, maxIter = 50,
                          affinityFunc = "gaussian", distFunc = "euclidean",
                          affinityParams = list(alpha = 1, c = 1, p = 2,
                                                Sigma = NULL),
                          mutationDecay = 1.0, mutationMin = 0.01,
                          maxClones = Inf, stopTolerance = 0.0,
                          noImprovementLimit = Inf,
                          initMethod = c("sample", "random",
                                         "random_uniform", "kmeans++"),
                          k = 3,
                          scale = c("none", "zscore", "robust", "arcsinh"),
                          targetK = NULL, epsilonQuantile = NULL,
                          verbose = TRUE, ...) {
  .validate_bHIVE_input(X, y)
  if (is.null(task)) {
    task <- if (is.null(y)) "clustering" else "classification"
  }
  task       <- match.arg(task, c("clustering", "classification"))
  initMethod <- match.arg(initMethod)
  scale      <- match.arg(scale)

  model <- AINet$new(
    nAntibodies = nAntibodies, beta = beta, epsilon = epsilon,
    maxIter = maxIter, k = k, affinityFunc = affinityFunc,
    distFunc = distFunc, affinityParams = affinityParams,
    mutationDecay = mutationDecay, mutationMin = mutationMin,
    maxClones = maxClones, stopTolerance = stopTolerance,
    noImprovementLimit = noImprovementLimit, initMethod = initMethod,
    scale = scale, targetK = targetK, epsilonQuantile = epsilonQuantile,
    verbose = verbose, ...)
  model$fit(X, y = y, task = task)

  # A plain list with the prototypes and labels, plus the fitted model so
  # callers can predict on new data.
  res <- model$result
  res$model <- model
  res
}

#' bHIVE: B-cell Hybrid Immune Variant Engine
#'
#' Implements an artificial immune network algorithm for clustering and
#' classification tasks. The algorithm evolves a population of "antibodies"
#' via clonal selection and mutation, applies network suppression to maintain
#' diversity, and assigns data points based on affinity or distance metrics.
#'
#' \code{bHIVE()} accepts a plain matrix or data frame (observations in rows)
#' or a \linkS4class{SummarizedExperiment} /
#' \linkS4class{SingleCellExperiment} (features in rows, cells in columns, the
#' Bioconductor convention). For the latter, the cells are clustered or
#' classified using an assay or, for a SingleCellExperiment, a
#' \code{reducedDim()} embedding, and the labels are stored in
#' \code{colData()}. To use bHIVE through \pkg{bluster}, see
#' \code{\link{bHIVEParam}}.
#'
#' @param X A numeric matrix or data frame of input features, with rows as
#' observations and columns as features. Or a SummarizedExperiment or
#' SingleCellExperiment with features in rows and cells in columns.
#' @param y Optional. A factor target vector for classification. If NULL,
#' clustering will be performed. For SummarizedExperiment input, \code{y} may
#' also be the name of a \code{colData()} column holding the labels.
#' @param task Character. Specifies the task to perform: \code{"clustering"} or
#' \code{"classification"}. If NULL, it is inferred based on \code{y}.
#' @param nAntibodies Integer. The initial population size of antibodies.
#' @param beta Numeric. Clone multiplier (controls how many clones are
#' generated for top-matching antibodies).
#' @param epsilon Numeric. Similarity threshold used in network suppression;
#' antibodies closer than \code{epsilon} are considered redundant.
#' @param maxIter Integer. Maximum number of iterations to run the AI-Net
#' algorithm.
#' @param affinityFunc Character. Specifies the affinity (similarity) function
#' to use for antibody-data matching. One of \code{"gaussian"}, \code{"laplace"},
#' \code{"polynomial"}, \code{"cosine"}, or \code{"hamming"}.
#' @param distFunc Character. Specifies the distance function for clustering
#' and suppression. One of \code{"euclidean"}, \code{"manhattan"},
#' \code{"minkowski"}, \code{"cosine"}, \code{"mahalanobis"},
#' or \code{"hamming"}.
#' @param affinityParams A list of optional parameters for the chosen affinity
#' or distance function.
#'   \itemize{
#'     \item \code{alpha} (for RBF or Laplace kernel). \code{"auto"} sets the
#'     bandwidth by the median heuristic, which is recommended for
#'     high-dimensional data such as mass cytometry.
#'     \item \code{c}, \code{p} (for polynomial kernel or Minkowski distance).
#'     The polynomial affinity is the bounded form
#'     \eqn{[(x \cdot a + c) / \max(x \cdot x + c, a \cdot a + c)]^p}, which
#'     lies in [0, 1] and peaks at \eqn{a = x}.
#'     \item \code{Sigma} (for Mahalanobis distance).
#'   }
#' @param mutationDecay Numeric. Factor by which the mutation rate decays each
#' iteration (should be \eqn{\le 1.0}). Default is 1.0 (no decay).
#' @param mutationMin Numeric. Minimum mutation rate, preventing the mutation
#' scale from shrinking to zero.
#' @param maxClones Numeric. Maximum number of clones per top-matching antibody;
#'  defaults to \code{Inf}.
#' @param stopTolerance Numeric. If the change in the number of antibodies
#' (repertoire size) is \eqn{\le stopTolerance} for consecutive iterations,
#' this may trigger the \code{noImprovementLimit}.
#' @param noImprovementLimit Integer. Stops the algorithm early if there is no
#' further improvement in antibody count (beyond \code{stopTolerance}) for this
#' many consecutive iterations. Default is \code{Inf}, meaning no early stop
#' based on improvement.
#' @param initMethod Character. Method for initializing antibodies. Can be:
#'   \itemize{
#'     \item \code{"sample"} - randomly selects rows from \code{X} as initial
#'     antibodies.
#'     \item \code{"random"} - samples Gaussian noise using \code{X}'s column
#'     means/sds.
#'     \item \code{"random_uniform"} - samples uniformly in [min, max] of each
#'     column.
#'     \item \code{"kmeans++"} - tries a kmeans++-like initialization for
#'     coverage.
#'   }
#' @param k Integer. Number of top-matching antibodies (by affinity) to
#' consider cloning for each data point.
#' @param scale Character. Per-feature input scaling: \code{"none"} (default),
#' \code{"zscore"}, \code{"robust"} (median/IQR), or \code{"arcsinh"} (CyTOF).
#' Passed to \code{\link{AINet}}; makes \code{epsilon} and distances behave
#' consistently across datasets of different magnitude.
#' @param targetK Integer or NULL. If set, force the clustering result to
#' exactly \code{targetK} clusters via a seeded K-means refinement (the immune
#' network supplies the seeds). NULL (default) keeps the emergent cluster count.
#' @param epsilonQuantile Numeric in (0, 1) or NULL. If set, the suppression
#' threshold adapts each iteration to this quantile of pairwise antibody
#' distances instead of the fixed \code{epsilon}.
#' @param verbose Logical. If \code{TRUE}, reports progress each iteration
#' with \code{message()}.
#' @param assay.type String or integer. The assay holding the features to
#' cluster, for SummarizedExperiment input. Ignored when \code{use.dimred} is
#' set.
#' @param use.dimred String or integer. The \code{reducedDim()} embedding to
#' cluster, for SingleCellExperiment input. Takes precedence over
#' \code{assay.type}. Recommended for sparse single-cell RNA-seq, since an
#' assay is converted to a dense matrix.
#' @param name String. The \code{colData()} column (and \code{metadata()}
#' entry) the results are stored under.
#' @param ... Additional arguments forwarded to \code{\link{AINet}}, including
#' immunology modules (\code{shm}, \code{idiotypic}, \code{germinalCenter},
#' \code{microenvironment}, \code{activation}, \code{memory},
#' \code{classSwitcher}, \code{init}), \code{consolidate} /
#' \code{consolidationSteps}, \code{coverageBoost} and, for classification,
#' \code{minPerClass} (antibodies seeded per class, default 5). For the
#' SummarizedExperiment methods, all arguments of the matrix method.
#'
#' @return For matrix or data frame input, a list:
#'   \itemize{
#'     \item \code{antibodies}: Final prototype matrix (prototypes x features),
#'     in the scaled space the model was trained in.
#'     \item \code{assignments}: For clustering, integer cluster IDs, where
#'     cluster \code{i} is row \code{i} of \code{antibodies}. For
#'     classification, the predicted class labels.
#'     \item \code{antibody_classes}: Classification only. The class label of
#'     each antibody.
#'     \item \code{antibody_index}: Classification only. The row of
#'     \code{antibodies} that labeled each observation.
#'     \item \code{antibodies_unscaled}: The prototypes back-transformed to
#'     the original feature units.
#'     \item \code{task}: The chosen task.
#'     \item \code{model}: The fitted \code{\link{AINet}} object. Use
#'     \code{res$model$predict(newdata)} to label new observations.
#'   }
#'   For SummarizedExperiment or SingleCellExperiment input, \code{X} with the
#'   labels added as a factor column \code{colData(X)[[name]]} and the full
#'   result list above in \code{metadata(X)[[name]]}.
#'
#' @seealso \code{\link{bHIVEParam}} for \pkg{bluster} integration,
#'   \code{\link{AINet}} for the underlying R6 engine.
#'
#' @examples
#' # Example 1: Clustering with the Iris dataset
#' data(iris)
#' X <- as.matrix(iris[, 1:4])  # Numeric features only
#' res <- bHIVE(X = X,
#'              task = "clustering",
#'              nAntibodies = 30,
#'              beta = 5,
#'              epsilon = 0.01,
#'              maxIter = 20,
#'              k = 3,
#'              verbose = FALSE)
#' table(res$assignments)
#'
#' # Example 2: Classification with Iris species
#' y <- iris$Species
#' res <- bHIVE(X = X,
#'               y = y,
#'               task = "classification",
#'               nAntibodies = 30,
#'               beta = 5,
#'               epsilon = 0.01,
#'               maxIter = 20,
#'               k = 3,
#'               verbose = FALSE)
#' table(res$assignments, y)
#'
#' # Example 3: SingleCellExperiment input (features x cells)
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(logcounts = t(X)),
#'                             colData = DataFrame(species = iris$Species))
#' sce <- bHIVE(sce, targetK = 3, maxIter = 10, verbose = FALSE)
#' table(sce$bHIVE, sce$species)
#'
#' # Cluster on an embedding, or classify with a colData column as labels
#' reducedDim(sce, "PCA") <- prcomp(X)$x[, 1:2]
#' sce <- bHIVE(sce, use.dimred = "PCA", name = "bHIVE_pca", verbose = FALSE)
#' sce <- bHIVE(sce, y = "species", name = "bHIVE_class", verbose = FALSE)
#' mean(sce$bHIVE_class == sce$species)
#'
#' @importFrom methods setGeneric setMethod is
#' @export
setGeneric("bHIVE", function(X, ...) standardGeneric("bHIVE"))

#' @rdname bHIVE
#' @export
setMethod("bHIVE", "ANY", .bHIVE_matrix)

#' @rdname bHIVE
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("bHIVE", "SummarizedExperiment",
          function(X, y = NULL, ..., assay.type = "logcounts",
                   name = "bHIVE") {
  mat <- t(as.matrix(assay(.bhive_check_assay(X, assay.type), assay.type)))
  .bHIVE_se_finish(X, mat, y, name, ...)
})

#' @rdname bHIVE
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDim reducedDimNames
#' @export
setMethod("bHIVE", "SingleCellExperiment",
          function(X, y = NULL, ..., assay.type = "logcounts",
                   use.dimred = NULL, name = "bHIVE") {
  if (is.null(use.dimred)) {
    mat <- t(as.matrix(assay(.bhive_check_assay(X, assay.type), assay.type)))
  } else {
    if (is.character(use.dimred) && !use.dimred %in% reducedDimNames(X)) {
      stop(sprintf("reducedDim '%s' not found. Available: %s", use.dimred,
                   paste(reducedDimNames(X), collapse = ", ")), call. = FALSE)
    }
    mat <- as.matrix(reducedDim(X, use.dimred))
  }
  .bHIVE_se_finish(X, mat, y, name, ...)
})


# ============================================================================
# Internal helpers
# ============================================================================

#' Names of every argument bHIVE can tune
#'
#' The union of the matrix method's arguments and AINet$initialize()'s, minus
#' the data and bookkeeping arguments. swarmbHIVE() uses this to decide which
#' grid columns to forward.
#' @keywords internal
#' @noRd
.bhive_tunable_args <- function() {
  setdiff(union(names(formals(.bHIVE_matrix)),
                names(formals(AINet$public_methods$initialize))),
          c("X", "y", "task", "verbose", "..."))
}

#' Check that an assay exists, with a readable error
#' @keywords internal
#' @noRd
.bhive_check_assay <- function(se, assay.type) {
  nms <- SummarizedExperiment::assayNames(se)
  if (is.character(assay.type) && !assay.type %in% nms) {
    stop(sprintf("assay '%s' not found. Available: %s", assay.type,
                 paste(nms, collapse = ", ")), call. = FALSE)
  }
  se
}

#' Resolve y for SummarizedExperiment input
#'
#' Accepts NULL (clustering), a vector of length ncol(se), or the name of a
#' colData() column. Returns a factor or NULL.
#' @keywords internal
#' @noRd
.bhive_resolve_y <- function(se, y) {
  if (is.null(y)) return(NULL)
  cd <- SummarizedExperiment::colData(se)
  if (is.character(y) && length(y) == 1L && y %in% colnames(cd)) {
    y <- cd[[y]]
  }
  if (length(y) != ncol(se)) {
    stop("y must be a colData() column name or have one value per column ",
         "(cell) of X.", call. = FALSE)
  }
  if (anyNA(y)) stop("y contains missing labels.", call. = FALSE)
  if (is.factor(y)) y else factor(y)
}

#' Fit bHIVE on cells extracted from a SummarizedExperiment and store results
#' @param mat Cells x features matrix.
#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors metadata metadata<-
.bHIVE_se_finish <- function(se, mat, y, name, ...) {
  y   <- .bhive_resolve_y(se, y)
  res <- .bHIVE_matrix(mat, y = y, ...)
  labels <- if (res$task == "classification") {
    factor(res$assignments, levels = levels(y))
  } else {
    factor(res$assignments)
  }
  colData(se)[[name]] <- labels
  metadata(se)[[name]] <- res
  se
}
