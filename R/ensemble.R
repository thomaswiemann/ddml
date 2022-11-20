#' Compute ensemble learners.
#'
#' Compute ensemble learners.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param type A string indicating the type of ensemble. Multiple types may be
#'     passed in form of a vector of strings.
#' @param learners A list of lists, each containing four named elements:
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
#' @param cv_results Optional output from \code{\link{crossval}}. If availbale, this
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
#' @param silent A boolean indicating whether current learners and folds should be
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
#'     corresponding to the estimators passed via \code{learners}.}
#' \item{\code{weights}}{A matrix where each column gives the ensemble weights
#'     for the corresponding ensemble as specified in \code{type}.}
#' \item{\code{learners}}{Passthrough of the input \code{learners}.}
#' \item{\code{cv_results}}{Passed output from \code{\link{crossval}}.}
#' \item{\code{mdl_w_iv}}{A vector of model indices corresponding to learners that
#'     selected at least one instrument. }
#' }
#'
#' @export ensemble
ensemble <- function(y, X, Z = NULL,
                     type = c("average"),
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     cv_results = NULL,
                     silent = F, progress = NULL) {
  # Data parameters
  nlearners <- length(learners)
  # Compute ensemble weights
  ens_w_res <- ensemble_weights(y, X, Z,
                                type = type, learners = learners,
                                cv_folds = cv_folds,
                                cv_subsamples = cv_subsamples,
                                cv_results = cv_results,
                                silent = silent, progress = progress)
  weights <- ens_w_res$weights
  cv_results <- ens_w_res$cv_results
  # Check for excluded learners
  mdl_include <- which(rowSums(abs(weights)) > 0)
  # Compute fit for each included model
  mdl_fits <- rep(list(NULL), nlearners)
  for (m in 1:nlearners) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Check whether X, Z assignment has been specified. If not, include all.
    if (is.null(learners[[m]]$assign_X))
      learners[[m]]$assign_X <- c(1:ncol(X))
    if (is.null(learners[[m]]$assign_Z) & !is.null(Z))
      learners[[m]]$assign_Z <- c(1:ncol(Z))
    # Else fit on data. Begin by selecting the model constructor and the
    #     variable assignment.
    mdl_fun <- list(what = learners[[m]]$fun, args = learners[[m]]$args)
    assign_X <- learners[[m]]$assign_X
    assign_Z <- learners[[m]]$assign_Z
    # Then fit the model
    mdl_fun$args$y <- y
    mdl_fun$args$X <- cbind(X[, assign_X],
                            Z[, assign_Z])
    mdl_fits[[m]] <- do.call(do.call, mdl_fun)
  }#FOR
  # Organize and return output
  output <- list(mdl_fits = mdl_fits, weights = weights,
                 learners = learners, cv_results = cv_results)
  class(output) <- "ensemble"
  return(output)
}#ENSEMBLE

# Complementary methods ========================================================
#' Predict method for ensemble learners.
#'
#' Predict method for ensemble learners.
#'
#' @export predict.ensemble
#' @export
predict.ensemble <- function(obj, newX = NULL, newZ = NULL){
  # Data parameters
  nlearners <- length(obj$mdl_fits)
  # Check for excluded learners
  mdl_include <- which(rowSums(abs(obj$weights)) > 0)
  # Calculate fitted values for each model
  first_fit <- T
  for (m in 1:nlearners) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Get assign_X and assing_Z
    assign_X <- obj$learners[[m]]$assign_X
    assign_Z <- obj$learners[[m]]$assign_Z
    # Compute predictions
    fitted <- predict(obj$mdl_fits[[m]],
                      newdata = cbind(newX[, assign_X],
                                      newZ[, assign_Z]))

    # Initialize matrix of fitted values
    if (first_fit) {
      fitted_mat <- matrix(0, length(fitted), nlearners)
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
#' @param learners A list of lists, each containing four named elements:
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
#' @param cv_results Optional output from \code{\link{crossval}}. If availbale, this
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
#' @param silent A boolean indicating whether current learners and folds should be
#'     printed to the console.
#'
#' @return \code{ensemble_weights} returns a list containing the following
#'     components:
#' \describe{
#' \item{\code{weights}}{A matrix where each column gives the ensemble weights
#'     for the corresponding ensemble as specified in \code{type}.}
#' \item{\code{cv_results}}{Passed output from \code{\link{crossval}}.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ensemble_weights
ensemble_weights <- function(y, X, Z = NULL,
                             type = c("average"),
                             learners,
                             cv_folds = 5,
                             cv_subsamples = NULL,
                             cv_results = NULL,
                             silent = F, progress = NULL) {
  # Data parameters
  nlearners <- length(learners)
  ntype <- length(type)
  # Check whether out-of-sample residuals should be calculated to inform the
  #     ensemble weights, and whether previous results are available.
  cv_stacking <- c("stacking", "stacking_nn", "stacking_01", "stacking_best")
  if (any(cv_stacking %in% type) & is.null(cv_results)) {
    # Run crossvalidation procedure
    cv_results <- crossval(y, X, Z,
                           learners = learners,
                           cv_folds = cv_folds,
                           cv_subsamples = cv_subsamples,
                           silent = silent, progress = progress)
  }#IF
  # Compute weights for each ensemble type
  weights <- matrix(0, length(learners), length(type))
  for (k in 1:ntype) {
    if (type[k] == "average") {
      # Assign 1 to all included learners and normalize
      weights[, k] <- 1
      weights[, k] <- weights[, k] / sum(weights[, k])
    } else if (type[k] == "stacking_01") {
      # For stacking with weights constrained between 0 and 1: |w|_1 = 1, solve
      # the quadratic programming problem.
      sq_resid <- Matrix::crossprod(cv_results$oos_resid)
      A <- matrix(1, 1, nlearners) # |w|_1 = 1 constraint
      # Define sink so LowRankQP output is not printed
      sink(file=nullfile())
      # Calculate solution
      r <- LowRankQP::LowRankQP(Vmat = sq_resid,
                                dvec = rep(0, nlearners),
                                Amat = A,
                                bvec = 1,
                                uvec = rep(1, nlearners),
                                method = 'LU',
                                verbose = F)
      weights[, k] <- r$alpha
      # Remove sink so output is no longer surpressed
      sink()
    } else if (type[k] == "stacking_nn") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_results$oos_resid
      # For non-negative stacking, calculate the non-negatuve ols coefficients
      weights[, k] <- nnls::nnls(oos_fitted, y)$x
    } else if (type[k] == "stacking") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_results$oos_resid
      # For unconstrained stacking, simply calculate the ols coefficients
      weights[, k] <- ols(y, oos_fitted)$coef
    } else if (type[k] == "stacking_best") {
      # Find MSPE-minimizing model
      mdl_min <- which.min(Matrix::colMeans(cv_results$oos_resid^2)[, drop = F])
      mdl_min <- c(1:nlearners)[mdl_min]
      # Assign unit weight to the best model
      weights[mdl_min, k] <- 1
    }#IFELSE
  }#FOR
  # Organize and return output
  output <- list(weights = weights, cv_results = cv_results)
  return(output)
}#ENSEMBLE_WEIGHTS
