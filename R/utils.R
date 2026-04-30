utils::globalVariables(c("PC1", "PC2", "Group", "Feature", "Prototype", "Layer"))

#' Null-coalesce operator
#'
#' Returns \code{x} if it is not \code{NULL}, otherwise returns \code{y}.
#' Used internally by bHIVE R6 classes for parameter defaults.
#'
#' @param x Value to test.
#' @param y Default value if \code{x} is NULL.
#' @return \code{x} if not NULL, otherwise \code{y}.
#' @name null-coalesce
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

.validate_bHIVE_input <- function(X, 
                                  y = NULL) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("Input X must be a matrix or data frame.")
  }
  if (!is.null(y)) {
    if (!is.factor(y)) {
      stop("y must be a factor (classification).")
    }
    if (nrow(X) != length(y)) {
      stop("X and y must have the same number of rows.")
    }
  }
  if (anyNA(X)) {
    stop("Input X contains missing values. Please handle missing values before running bHIVE.")
  }
  invisible(TRUE)
}

#Determines how to move a single prototype 'ab_vec' in feature space
.update_prototype <- function(ab_vec,
                              x_i,
                              same_label = NA,
                              task = c("clustering","classification"),
                              loss = c("categorical_crossentropy", "binary_crossentropy",
                                       "kullback_leibler", "cosine", "mae"),
                              push_away = TRUE) {
  task <- match.arg(task)
  loss <- match.arg(loss)

  # must return a numeric vector of length = length(ab_vec)
  d <- length(ab_vec)
  zero_vec <- numeric(d)          # a zero vector for quick returns
  pull_vec <- x_i - ab_vec        # pulling ab_vec toward x_i
  push_vec <- -pull_vec           # pushing ab_vec away from x_i

  # =================
  # 1) CLUSTERING
  # =================
  if (task == "clustering") {
    return(pull_vec)
  }

  # =================
  # 2) CLASSIFICATION
  # =================
  if (task == "classification") {
    # If no assigned label => no movement
    if (is.na(same_label)) {
      return(zero_vec)
    }
    is_same <- isTRUE(same_label)

    # (A) If 'is_same' => we PULL
    # (B) Else if push_away => we PUSH
    # (C) else => zero
    if (loss %in% c("categorical_crossentropy","binary_crossentropy","kullback_leibler","cosine")) {
      if (is_same) {
        return(pull_vec)
      } else if (push_away) {
        return(push_vec)
      } else {
        return(zero_vec)
      }
    }

    else if (loss == "mae") {
      # sign-based approach
      if (is_same) {
        return(sign(pull_vec))
      } else if (push_away) {
        return(sign(push_vec))
      } else {
        return(zero_vec)
      }
    }

    return(zero_vec)
  }

  # If something unexpected
  return(zero_vec)
}
