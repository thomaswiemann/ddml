#' Compute ensemble models.
#'
#' Compute ensemble models.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param type A string indicating the type of ensemble. Multiple types may be
#'     passed in form of a vector of strings.
#' @param models A list of lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}}
#'     }
#' @param cv_folds The number for cross-validation folds.
#' @param cv_res Optional output from \code{\link{crossval}}. If availbale, this
#'     avoids unnecessary computation.
#' @param setup_parallel An list containing two named elements:
#'     \itemize{
#'         \item{\code{type} A string of value \code{"static"} or
#'             \code{"dynamic"}, indicating whether job scheduling should be
#'             static (via \code{parSapply}) or dynamic (via \code{foreach}).}
#'         \item{\code{cores} The number of processor units utilized. Note that
#'             if \code{cores == 1}, \code{crossval} will not be computed in
#'             parallel.
#'         }
#'     }
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{ensemble} returns an object of S3 class "\code{ensemble}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class.
#'
#' An object of class "\code{ensemble}" is a list containig the following
#'     components:
#' \describe{
#' \item{\code{mdl_fits}}{A list containing the fitted model objects
#'     corresponding to the estimators passed via \code{models}.}
#' \item{\code{weights}}{A matrix where each column gives the ensemble weights
#'     for the corresponding ensemble as specified in \code{type}.}
#' \item{\code{models}}{Passthrough of the input \code{models}.}
#' \item{\code{cv_res}}{Passed output from \code{\link{crossval}}.}
#' \item{\code{mdl_w_iv}}{A vector of model indices corresponding to models that
#'     selected at least one instrument. }
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ensemble
ensemble <- function(y, X, Z = NULL,
                     type = c("average"),
                     models,
                     cv_folds = 5,
                     cv_res = NULL,
                     setup_parallel = list(type = 'dynamic', cores = 1),
                     silent = F) {
  # Data parameters
  nmodels <- length(models)
  # Compute ensemble weights
  ens_w_res <- ensemble_weights(y, X, Z,
                                type, models,
                                cv_folds, cv_res, setup_parallel, silent)
  weights <- ens_w_res$weights
  cv_res <- ens_w_res$cv_res
  # Check for excluded models
  mdl_include <- which(rowSums(abs(weights)) > 0)
  # Compute fit for each included model
  mdl_fits <- rep(list(NULL), nmodels)
  ens_Z <- rep(FALSE, nmodels)
  for (m in 1:nmodels) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Else fit on data. Begin by selecting the model constructor and the
    #     variable assignment.
    mdl_fun <- list(what = models[[m]]$fun, args = models[[m]]$args)
    assign_X <- models[[m]]$assign_X
    assign_Z <- models[[m]]$assign_Z
    # Then fit the model
    mdl_fun$args$y <- y
    mdl_fun$args$X <- cbind(X[, assign_X],
                            Z[, assign_Z])
    mdl_fits[[m]] <- do.call(do.call, mdl_fun)
    # Check whether instruments were selected

    if (!is.null(Z)) {
      index_iv <- (length(assign_X) + 1):length(c(assign_X, assign_Z))
      ens_Z[m] <- any_iv(obj = mdl_fits[[m]],
                         index_iv = index_iv,
                         names_iv = colnames(Z[1:2, assign_Z]))
    } else {
      ens_Z[m] <- TRUE
    }#IFELSE
  }#FOR
  # Check whether (further) models should be excluded based on IV selection
  mdl_w_iv <- setdiff(cv_res$cv_Z, setdiff(c(1:nmodels), which(ens_Z)))
  if (length(mdl_w_iv) != length(cv_res$cv_Z)) {
    # Recompute weights on smaller set of admissable models
    cv_res_adj <- cv_res
    cv_res_adj$cv_Z <- mdl_w_iv
    weights <- ensemble_weights(y, X, Z,
                                type, models,
                                cv_folds, cv_res_adj,
                                setup_parallel, silent)$weights
  }#IF
  # Organize and return output
  output <- list(mdl_fits = mdl_fits, weights = weights,
                 models = models, cv_res = cv_res, mdl_w_iv = mdl_w_iv)
  class(output) <- "ensemble"
  return(output)
}#ENSEMBLE

# Complementary methods ========================================================
#' Predict method for ensemble models.
#'
#' Predict method for ensemble models.
#'
#' @export predict.ensemble
#' @export
predict.ensemble <- function(obj, newX = NULL, newZ = NULL){
  # Data parameters
  nmodels <- length(obj$mdl_fits)
  # Check for excluded models
  mdl_include <- which(rowSums(abs(obj$weights)) > 0)
  # Calculate fitted values for each model
  first_fit <- T
  for (m in 1:nmodels) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Get assign_X and assing_Z
    assign_X <- obj$models[[m]]$assign_X
    assign_Z <- obj$models[[m]]$assign_Z
    # Compute predictions
    fitted <- predict(obj$mdl_fits[[m]],
                      newdata = cbind(newX[, assign_X, drop = F],
                                      newZ[, assign_Z, drop = F]))

    # Initialize matrix of fitted values
    if (first_fit) {
      fitted_mat <- matrix(0, length(fitted), nmodels)
      first_fit <- F
    }#IF
    fitted_mat[, m] <- as(fitted, "matrix")
  }#FOR
  # Compute matrix of fitted values by ensemble type and return
  fitted_ens <- fitted_mat %*% obj$weights
  return(fitted_ens)
}#PREDICT.ENSEMBLE

# Complementary functions ======================================================
#' Compute ensemble weights.
#'
#' Compute ensemble weights.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param type A string indicating the type of ensemble. Multiple types may be
#'     passed in form of a vector of strings.
#' @param models A list of lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}}
#'     }
#' @param cv_folds The number for cross-validation folds.
#' @param cv_res Optional output from \code{\link{crossval}}. If availbale, this
#'     avoids unnecessary computation.
#' @param setup_parallel An list containing two named elements:
#'     \itemize{
#'         \item{\code{type} A string of value \code{"static"} or
#'             \code{"dynamic"}, indicating whether job scheduling should be
#'             static (via \code{parSapply}) or dynamic (via \code{foreach}).}
#'         \item{\code{cores} The number of processor units utilized. Note that
#'             if \code{cores == 1}, \code{crossval} will not be computed in
#'             parallel.
#'         }
#'     }
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{ensemble_weights} returns a list containing the following
#'     components:
#' \describe{
#' \item{\code{weights}}{A matrix where each column gives the ensemble weights
#'     for the corresponding ensemble as specified in \code{type}.}
#' \item{\code{cv_res}}{Passed output from \code{\link{crossval}}.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ensemble_weights
ensemble_weights <- function(y, X, Z = NULL,
                             type = c("average"),
                             models,
                             cv_folds = 5,
                             cv_res = NULL,
                             setup_parallel = list(type = 'dynamic', cores = 1),
                             silent = F) {
  # Data parameters
  nmodels <- length(models)
  ntype <- length(type)
  # Check whether out-of-sample residuals should be calculated to inform the
  #     ensemble weights, and whether previous results are available.
  if (any(c("stacking", "stacking_nn", "stacking_01", "cv") %in% type) &&
      (is.null(cv_res))) {
    # Run crossvalidation procedure
    cv_res <- crossval(y, X, Z, models, cv_folds, setup_parallel, silent)
  }#IF
  # Compute weights for each ensemble type
  weights <- matrix(0, length(models), length(type))
  for (k in 1:ntype) {
    # Check if multiple models were selected. If not, w is a unit vector.
    if (length(cv_res$cv_Z) == 1 & "average" != type[k]) {
      weights[cv_res$cv_Z, k] <- 1
      next
    }#IF
    if (type[k] == "average") {
      # Assign 1 to all included models and normalize
      weights[, k] <- 1
      weights[, k] <- weights[, k] / sum(weights[, k])
    } else if (type[k] == "stacking_01") {
      # For stacking with weights constrained between 0 and 1: |w|_1 = 1, solve
      # the quadratic programming problem.
      sq_resid <- Matrix::crossprod(cv_res$oos_resid[, cv_res$cv_Z])
      ncv_Z <- length(cv_res$cv_Z)
      A <- matrix(1, 1, ncv_Z) # |w|_1 = 1 constraint
      # Define sink so LowRankQP output is not printed
      sink(file=file())
      # Calculate solution
      r <- LowRankQP::LowRankQP(Vmat = sq_resid,
                                  dvec = rep(0, ncv_Z),
                                  Amat = A,
                                  bvec = 1,
                                  uvec = rep(1, ncv_Z),
                                  method = 'LU',
                                  verbose = F)
      weights[cv_res$cv_Z, k] <- r$alpha
      # Remove sink so output is no longer surpressed
      sink()
    } else if (type[k] == "stacking_nn") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_res$oos_resid[, cv_res$cv_Z, drop = F]
      # For non-negative stacking, calculate the non-negatuve ols coefficients
      weights[cv_res$cv_Z, k] <- nnls::nnls(oos_fitted, y)$x
    } else if (type[k] == "stacking") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_res$oos_resid[, cv_res$cv_Z, drop = F]
      # For unconstrained stacking, simply calculate the ols coefficients
      weights[cv_res$cv_Z, k] <- ols(y, oos_fitted)$coef
    } else if (type[k] == "cv") {
      # Find MSPE-minimizing model
      mdl_min <- which.min(Matrix::colMeans(cv_res$oos_resid^2)[cv_res$cv_Z,
                                                                drop = F])
      mdl_min <- c(1:nmodels)[cv_res$cv_Z][mdl_min]
      # Assign unit weight to the best model
      weights[mdl_min, k] <- 1
    }#IFELSE
  }#FOR
  # Organize and return output
  output <- list(weights = weights, cv_res = cv_res)
  return(output)
}#ENSEMBLE_WEIGHTS
