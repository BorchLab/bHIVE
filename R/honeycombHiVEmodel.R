#' Mulilayered honeycombHIVE for caret
#'
#' A \code{caret} wrapper for the \code{\link{honeycombHIVE}} function,
#' enabling seamless integration with the \code{caret} package for
#' hyperparameter tuning, cross-validation, and performance evaluation.
#'
#' @name honeycombHIVEmodel
#' @section Parameters:
#'   \itemize{
#'     \item \code{nAntibodies}: Number of initial antibodies in the network.
#'     \item \code{beta}: Clone multiplier controlling the number of clones per
#'           antibody.
#'     \item \code{epsilon}: Threshold for network suppression to remove
#'           redundant antibodies.
#'     \item \code{layers}: Number of hierarchical layers for iterative
#'           refinement.
#'     \item \code{refineOptimizer}: Optimizer for gradient-based refinement
#'           (e.g. "sgd", "momentum", "adagrad", "adam", "rmsprop").
#'     \item \code{refineSteps}: Number of gradient update steps in refinement.
#'     \item \code{refineLR}: Learning rate for refinement.
#'   }
#'
#' @section Supported Tasks:
#'   \itemize{
#'     \item \code{"Classification"}: Assigns class labels to input observations.
#'     \item \code{"Clustering"}: Groups data points based on similarity
#'           (though typically \code{caret} is used for supervised tasks).
#'   }
#'
#' @return A \code{caret} model definition list. Pass it to
#'   \code{\link[caret]{train}} for model training and evaluation.
#'
#' @examples
#' \dontrun{
#'   library(caret)
#'   # Example: Classification with Iris
#'   data(iris)
#'   X <- as.matrix(iris[, 1:4])
#'   y <- iris$Species
#'
#'   train_control <- trainControl(method = "cv", number = 5)
#'   set.seed(42)
#'   model <- train(
#'     x = X,
#'     y = y,
#'     method = honeycombHIVEmodel,
#'     trControl = train_control,
#'     tuneGrid = expand.grid(
#'       nAntibodies = c(10, 20),
#'       beta = c(3, 5),
#'       epsilon = c(0.01, 0.05),
#'       layers = c(1, 2),
#'       refineOptimizer = "adam",
#'       refineSteps = 5,
#'       refineLR = 0.01
#'     )
#'   )
#'   print(model)
#' }
#'
#' @seealso \code{\link[caret]{train}}, \code{\link[caret]{trainControl}}
#' @export
honeycombHIVEmodel <- list(
  label = "Multilayered Artificial Immune Network (honeycombHIVE)",
  library = c("bHIVE"),
  type = c("Classification", "Clustering"),

  # Tuning parameters, including new refinement settings
  parameters = data.frame(
    parameter = c("nAntibodies", "beta", "epsilon", "layers",
                  "refineOptimizer", "refineSteps", "refineLR"),
    class = c("numeric", "numeric", "numeric", "integer",
              "character", "numeric", "numeric"),
    label = c(
      "Number of Antibodies",
      "Clone Multiplier (beta)",
      "Epsilon Threshold",
      "Number of Layers",
      "Refinement Optimizer",
      "Refinement Steps",
      "Refinement Learning Rate"
    )
  ),

  # Default tuning grid: additional refinement parameters are set to fixed defaults
  grid = function(x, y, len = NULL) {
    expand.grid(
      nAntibodies = seq(10, 50, length.out = len),
      beta = seq(1, 10, length.out = len),
      epsilon = seq(0.01, 0.1, length.out = len),
      layers = seq(1, 3, length.out = len),
      refineOptimizer = "adam",
      refineSteps = 5,
      refineLR = 0.01
    )
  },

  # The core "fit" function: trains honeycombHIVE on the given data
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    task <- if (is.factor(y)) "classification" else "clustering"

    # Train the honeycombHIVE model with refinement parameters passed in from the tuning grid
    results <- honeycombHIVE(
      X = x,
      y = y,
      task = task,
      nAntibodies = param$nAntibodies,
      beta = param$beta,
      epsilon = param$epsilon,
      layers = param$layers,
      refine = TRUE,
      refineLoss = "categorical_crossentropy",
      refineSteps = param$refineSteps,
      refineLR = param$refineLR,
      refinePushAway = (task == "classification"),
      refineOptimizer = param$refineOptimizer,
      # Fixed defaults for additional hyperparameters (could be tuned in an extended grid)
      refineMomentumCoef = 0.9,
      refineBeta1 = 0.9,
      refineBeta2 = 0.999,
      refineRmspropDecay = 0.9,
      refineEpsilon = 1e-8,
      ...
    )

    # The model object returned to caret ("modelFit") contains:
    # 1) The entire results list from honeycombHIVE
    # 2) The final layer's output (for prediction purposes)
    final_layer <- results[[length(results)]]

    modelFit <- list(
      task = task,
      results = results,      # all layers
      finalLayer = final_layer,  # convenience
      obsLevels = if (task == "classification") lev else NULL
    )

    if (task == "classification") {
      modelFit$levels <- lev
    } else {
      modelFit$levels <- NULL
    }

    class(modelFit) <- "honeycombHIVEmodel"
    return(modelFit)
  },

  # Prediction function for new data
  predict = function(modelFit, newdata, submodels = NULL) {
    task <- modelFit$task
    finalLayer <- modelFit$finalLayer

    # For clustering, prediction is not defined in this wrapper.
    if (task == "clustering") {
      stop("Prediction not defined for clustering in this caret wrapper.")
    }

    # 1) Get final antibodies from the final layer
    finalAb <- finalLayer$antibodies
    if (is.null(finalAb)) {
      stop("No antibodies found in final layer. Check your honeycombHIVE output.")
    }

    # 2) Get the predicted label for each antibody from the final layer's predictions
    antibodyLabels <- finalLayer$predictions
    if (is.null(antibodyLabels)) {
      stop("No final antibody labels found. Check how your honeycombHIVE stores them.")
    }

    # 3) Compute distances from newdata to finalAb (using Euclidean squared distance)
    newdata <- as.matrix(newdata)
    finalAb <- as.matrix(finalAb)
    preds <- vector("list", nrow(newdata))
    for (i in seq_len(nrow(newdata))) {
      diffs <- sweep(finalAb, 2, newdata[i, ], FUN = "-")
      dists <- rowSums(diffs^2)
      idx_min <- which.min(dists)
      preds[[i]] <- idx_min
    }
    nearest_idx <- unlist(preds)

    # 4) Return the corresponding label for each observation
    out <- antibodyLabels[nearest_idx]
    out <- factor(out, levels = modelFit$levels)
    return(out)
  },

  # Probability predictions for classification are not supported in this wrapper.
  prob = function(modelFit, newdata, submodels = NULL) {
    stop("Probability predictions are not supported for honeycombHIVE in this wrapper.")
  }
)
