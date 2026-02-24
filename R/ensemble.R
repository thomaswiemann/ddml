#' Stacking estimator using combinations of base learners
#'
#' @family utilities
#'
#' @description Computes an ensemble of learners based on the specified
#'     aggregation type and computes cross-validated out-of-sample
#'     predictions to inform the weights.
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param Z Optional instrumental variables matrix.
#' @param type A character string indicating the type of ensemble to compute.
#'     Default is \code{"average"}.
#' @param learners A list of base learners.
#' @param cv_folds Number of cross-validation folds.
#' @param cv_subsamples Optional list of subsamples for cross-validation.
#' @param cv_results Optional pre-computed cross-validation results.
#' @param custom_weights Optional custom weights matrix.
#' @param silent A boolean indicating whether to suppress progress messages.
#'
#' @return An object of class \code{ensemble} containing:
#'     \item{mdl_fits}{List of fitted base learners.}
#'     \item{weights}{Computed ensemble weights.}
#'     \item{learners}{The base learners used.}
#'     \item{cv_results}{Cross-validation results if computed.}
#'     \item{mean_y}{Mean of the outcome variable.}
#'     \item{constant_y}{Boolean indicating if y is constant.}
#' @export
#'
#' @examples
#' \donttest{
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#'
#' # Fit an ensemble of ols, lasso, and ridge
#' ens_fit = ensemble(y, X,
#'                    type = "nnls",
#'                    learners = list(list(fun = ols),
#'                                   list(fun = mdl_glmnet),
#'                                   list(fun = mdl_glmnet,
#'                                        args = list(alpha = 0))),
#'                    cv_folds = 5,
#'                    silent = TRUE)
#' ens_fit$weights
#' predict(ens_fit, newdata = X)[1:5]
#' }
ensemble <- function(y, X, Z = NULL,
                     type = "average",
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     cv_results = NULL,
                     custom_weights = NULL,
                     silent = FALSE) {
  # Data parameters
  nlearners <- length(learners)
  # Check if y is constant
  if (length(unique(y)) == 1) {
    warning(paste("Outcome variable y is constant. Ensemble will return",
                   "mean(y) for all predictions."))
    # Return minimal output needed for predictions
    output <- list(
      mdl_fits = NULL,
      weights = NULL,
      learners = learners,
      cv_results = NULL,
      mean_y = mean(y),
      constant_y = TRUE
    )
    class(output) <- "ensemble"
    return(output)
  }#IF
  # Compute ensemble weights
  ens_w_res <- ensemble_weights(y, X, Z,
                                type = type, learners = learners,
                                cv_folds = cv_folds,
                                cv_subsamples = cv_subsamples,
                                cv_results = cv_results,
                                custom_weights = custom_weights,
                                silent = silent)
  weights <- ens_w_res$weights
  cv_results <- ens_w_res$cv_results
  # Check for excluded learners
  mdl_include <- which(rowSums(abs(weights)) > 0)
  if (length(mdl_include) == 0) {
    warning("None of the learners are assigned positive stacking weights.")
  }#IF
  # Compute fit for each included model
  mdl_fits <- rep(list(NULL), nlearners)
  for (m in seq_len(nlearners)) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Check whether X, Z assignment has been specified. If not, include all.
    if (is.null(learners[[m]]$assign_X))
      learners[[m]]$assign_X <- 1:ncol(X)
    if (is.null(learners[[m]]$assign_Z) & !is.null(Z))
      learners[[m]]$assign_Z <- 1:ncol(Z)
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
                 learners = learners, cv_results = cv_results,
                 mean_y = mean(y), constant_y = FALSE)
  class(output) <- "ensemble"
  return(output)
}#ENSEMBLE

# Complementary methods ========================================================

#' @exportS3Method
predict.ensemble <- function(object, newdata, newZ = NULL, ...){
  # Data parameters
  nlearners <- length(object$learners)
  # If y was constant, return mean_y for all observations
  if (!is.null(object$constant_y) && object$constant_y) {
    return(matrix(object$mean_y, nrow(newdata), nlearners))
  }#IF
  # Check for excluded learners
  mdl_include <- which(rowSums(abs(object$weights)) > 0)
  if (length(mdl_include) == 0) {
    fitted_ens <- matrix(object$mean_y, nrow(newdata), 1)
    #warning("None of the learners are assigned positive stacking weights.")
  }#IF
  # Calculate fitted values for each model
  first_fit <- TRUE
  for (m in seq_len(nlearners)) {
    # Skip model if not assigned positive weight
    if (!(m %in% mdl_include)) next
    # Get assign_X and assing_Z
    assign_X <- object$learners[[m]]$assign_X
    assign_Z <- object$learners[[m]]$assign_Z
    # Compute predictions
    fitted <- stats::predict(object$mdl_fits[[m]],
                             newdata = cbind(newdata[, assign_X],
                                             newZ[, assign_Z]))

    # Initialize matrix of fitted values
    if (first_fit) {
      fitted_mat <- matrix(0, length(fitted), nlearners)
      first_fit <- FALSE
    }#IF
    fitted_mat[, m] <- methods::as(fitted, "matrix")
  }#FOR
  # Compute matrix of fitted values by ensemble type and return
  if (length(mdl_include) > 0) fitted_ens <- fitted_mat %*% object$weights
  return(fitted_ens)
}#PREDICT.ENSEMBLE

# Complementary functions ======================================================

#' Compute Stacking Weights for Base Learners
#'
#' @family utilities
#'
#' @description Computes the stacking weights for an ensemble of base learners
#'     using cross-validated out-of-sample predictions.
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param Z Optional instrumental variables matrix.
#' @param type A character string or vector indicating the type(s) of ensemble
#'     weights to compute. Default is \code{"average"}.
#' @param learners A list of base learners.
#' @param cv_folds Number of cross-validation folds.
#' @param cv_subsamples Optional list of subsamples for cross-validation.
#' @param cv_results Optional pre-computed cross-validation results.
#' @param custom_weights Optional custom weights matrix.
#' @param silent A boolean indicating whether to suppress progress messages.
#'
#' @return A list containing:
#'     \item{weights}{A matrix of computed ensemble weights.}
#'     \item{cv_results}{Cross-validation results used for computing weights.}
#' @export
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#'
#' # Compute stacking weights via NNLS
#' ew = ensemble_weights(y, X,
#'                       type = "nnls",
#'                       learners = list(list(fun = ols),
#'                                      list(fun = mdl_glmnet)),
#'                       cv_folds = 5,
#'                       silent = TRUE)
#' ew$weights
#' }
ensemble_weights <- function(y, X, Z = NULL,
                             type = "average",
                             learners,
                             cv_folds = 5,
                             cv_subsamples = NULL,
                             cv_results = NULL,
                             custom_weights = NULL,
                             silent = FALSE) {
  # Data parameters
  nlearners <- length(learners)
  ncustom <- ncol(custom_weights)
  ncustom <- ifelse(is.null(ncustom), 0, ncustom)
  ntype <- length(type)

  # Check whether out-of-sample residuals should be calculated to inform the
  #     ensemble weights, and whether previous results are available.
  cv_stacking <- c("ols", "nnls", "nnls1", "singlebest")
  if (any(cv_stacking %in% type) & is.null(cv_results)) {
    # Run crossvalidation procedure
    cv_results <- crossval(y, X, Z,
                           learners = learners,
                           cv_folds = cv_folds,
                           cv_subsamples = cv_subsamples,
                           silent = silent)
  }#IF
  # Compute weights for each ensemble type
  weights <- matrix(0, nlearners, ntype + ncustom)
  for (k in seq_len(ntype)) {
    if (type[k] == "average") {
      # Assign 1 to all included learners and normalize
      weights[, k] <- 1
      weights[, k] <- weights[, k] / sum(weights[, k])
    } else if (type[k] == "nnls1") {
      # For stacking with weights constrained between 0 and 1: |w|_1 = 1, solve
      # the quadratic programming problem.
      sq_resid <- Matrix::crossprod(cv_results$oos_resid)
      A <- cbind(matrix(1, nlearners, 1), diag(1, nlearners))
      # Calculate solution
      # Note: quadprog only solves for pos.def matrices. nearPD finds nearest
      #     pos.def matrix as a workaround.
      r <- quadprog::solve.QP(Dmat = Matrix::nearPD(sq_resid)$mat,
                             dvec = matrix(0, nlearners, 1),
                             Amat = A,
                             bvec = c(1, rep(0, nlearners)))
      weights[, k] <- r$solution
    } else if (type[k] == "nnls") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_results$oos_resid
      # For non-negative stacking, calculate the non-negatuve ols coefficients
      weights[, k] <- nnls::nnls(oos_fitted, y)$x
    } else if (type[k] == "ols") {
      # Reconstruct out of sample fitted values
      oos_fitted <- as.numeric(y) - cv_results$oos_resid
      # For unconstrained stacking, simply calculate the ols coefficients
      weights[, k] <- ols(y, oos_fitted, const = FALSE)$coef
    } else if (type[k] == "singlebest") {
      # Find MSPE-minimizing model
      mdl_min <- which.min(Matrix::colMeans(cv_results$oos_resid^2)[, drop = FALSE])
      mdl_min <- (seq_len(nlearners))[mdl_min]
      # Assign unit weight to the best model
      weights[mdl_min, k] <- 1
    }#IFELSE
  }#FOR
  # Append weights with custom weights
  if (!(ncustom == 0)) {
    weights[, (ntype + 1):(ntype + ncustom)] <- custom_weights
  }#IF
  # Assign ensemble types to columns
  if (!(ncustom == 0) && is.null(colnames(custom_weights))) {
    colnames(custom_weights) <- paste0("custom_", 1:ncustom)
  }#IF
  colnames(weights) <- c(type, colnames(custom_weights))
  # Organize and return output
  output <- list(weights = weights, cv_results = cv_results)
  return(output)
}#ENSEMBLE_WEIGHTS
