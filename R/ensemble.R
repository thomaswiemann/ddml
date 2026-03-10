#' Stacking Estimator Using Combinations of Base Learners
#'
#' @family utilities
#'
#' @description Computes an ensemble of learners based on the specified
#'     aggregation type and computes cross-validated out-of-sample
#'     predictions to inform the weights.
#'
#' @param y The outcome variable.
#' @param X The feature matrix.
#' @param type A character string indicating the type of ensemble to compute.
#'     Default is \code{"average"}.
#' @param learners A list of base learners. See
#'     \code{\link{ddml-intro}} for the full specification.
#' @param cv_folds Number of cross-validation folds.
#' @param cv_subsamples Optional list of subsamples for cross-validation.
#' @param cv_results Optional pre-computed cross-validation results.
#' @param custom_weights Optional custom weights matrix.
#' @param silent A boolean indicating whether to suppress progress messages.
#'
#' @return An object of class \code{ensemble} containing:
#'     \describe{
#'         \item{\code{mdl_fits}}{List of fitted base learners.}
#'         \item{\code{weights}}{Computed ensemble weights.}
#'         \item{\code{learners}}{The base learners used.}
#'         \item{\code{cv_results}}{Cross-validation results if
#'             computed.}
#'         \item{\code{mean_y}}{Mean of the outcome variable.}
#'         \item{\code{constant_y}}{Boolean indicating if y is
#'             constant.}
#'     }
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
#'                    learners = list(list(what = ols),
#'                                   list(what = mdl_glmnet),
#'                                   list(what = mdl_glmnet,
#'                                        args = list(alpha = 0))),
#'                    cv_folds = 5,
#'                    silent = TRUE)
#' ens_fit$weights
#' predict(ens_fit, newdata = X)[1:5]
#' }
ensemble <- function(y, X,
                     type = "average",
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     cv_results = NULL,
                     custom_weights = NULL,
                     silent = FALSE) {
  # Normalize learner specs
  learners <- normalize_learners(learners)

  # Wrap single-learner spec into a 1-element list
  if (is_single_learner(learners)) {
    learners <- list(learners)
  }#IF

  # Data parameters
  nlearners <- length(learners)

  # Constant-y: return trivial average weights of correct shape
  if (length(unique(y)) == 1) {
    warning(paste("Outcome variable y is constant. Ensemble will",
                  "return mean(y) for all predictions."),
            call. = FALSE)
    nensb <- length(type) + n_custom(custom_weights)
    weights <- matrix(1 / nlearners, nlearners, nensb)
    colnames(weights) <- c(type, colnames(custom_weights))
    output <- list(
      mdl_fits = NULL,
      weights = weights,
      learners = learners,
      cv_results = NULL,
      mean_y = mean(y),
      constant_y = TRUE)
    class(output) <- "ensemble"
    return(output)
  }#IF

  # Compute ensemble weights
  ens_w_res <- ensemble_weights(y, X,
                                type = type, learners = learners,
                                cv_folds = cv_folds,
                                cv_subsamples = cv_subsamples,
                                cv_results = cv_results,
                                custom_weights = custom_weights,
                                silent = silent)
  weights <- ens_w_res$weights
  cv_results <- ens_w_res$cv_results
  # Warn if all learner weights are zero across ensemble columns
  if (!any(rowSums(abs(weights)) > 0)) {
    warning("None of the learners are assigned positive stacking ",
            "weights.", call. = FALSE)
  }#IF

  # Fit all base learners
  mdl_fits <- rep(list(NULL), nlearners)
  for (m in seq_len(nlearners)) {
    if (is.null(learners[[m]]$assign_X))
      learners[[m]]$assign_X <- seq_len(ncol(X))
    mdl_fits[[m]] <- fit_learner(learners[[m]], y, X)
  }#FOR

  # Organize and return output
  output <- list(mdl_fits = mdl_fits, weights = weights,
                 learners = learners, cv_results = cv_results,
                 mean_y = mean(y), constant_y = FALSE)
  class(output) <- "ensemble"
  return(output)
}#ENSEMBLE

# Complementary methods ========================================================

#' Predict Method for \code{ensemble} Objects
#'
#' @param object A fitted \code{ensemble} object.
#' @param newdata A feature matrix for prediction.
#' @param ... Currently unused.
#' @param type Character; \code{"ensemble"} (default) returns
#'     weighted ensemble predictions, \code{"bylearner"} returns
#'     the raw per-learner prediction matrix.
#'
#' @return A matrix of predictions. When \code{type = "ensemble"},
#'     one column per ensemble type; when \code{type = "bylearner"},
#'     one column per base learner.
#'
#' @exportS3Method
predict.ensemble <- function(object, newdata, ...,
                             type = "ensemble") {
  type <- match.arg(type, c("ensemble", "bylearner"))
  nlearners <- length(object$learners)
  # Constant-y: return mean_y with the correct number of columns
  if (!is.null(object$constant_y) && object$constant_y) {
    ncols <- if (type == "bylearner") nlearners
      else ncol(object$weights)
    return(matrix(object$mean_y, nrow(newdata), ncols))
  }#IF
  # Calculate fitted values for each learner
  fitted_mat <- matrix(0, nrow(newdata), nlearners)
  for (m in seq_len(nlearners)) {
    assign_X <- object$learners[[m]]$assign_X
    fitted <- stats::predict(object$mdl_fits[[m]],
                             newdata = newdata[, assign_X,
                                               drop = FALSE])
    fitted_mat[, m] <- methods::as(fitted, "matrix")
  }#FOR
  if (type == "bylearner") return(fitted_mat)
  fitted_mat %*% object$weights
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
#' @param type A character string or vector indicating the type(s) of ensemble
#'     weights to compute. Default is \code{"average"}.
#' @param learners Optional list of base learners.
#'     Required when \code{cv_results} is not supplied (learners are
#'     needed to run cross-validation). When \code{cv_results} is
#'     supplied, \code{learners} may be omitted; the number of
#'     learners is inferred from the cross-validation residuals.
#'     See \code{\link{ddml-intro}} for the full specification.
#' @param cv_folds Number of cross-validation folds.
#' @param cv_subsamples Optional list of subsamples for cross-validation.
#' @param cv_results Optional pre-computed cross-validation results.
#' @param custom_weights Optional custom weights matrix.
#' @param silent A boolean indicating whether to suppress progress messages.
#'
#' @return A list containing:
#'     \describe{
#'         \item{\code{weights}}{A matrix of computed ensemble
#'             weights.}
#'         \item{\code{cv_results}}{Cross-validation results used
#'             for computing weights.}
#'     }
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
#'                       learners = list(list(what = ols),
#'                                      list(what = mdl_glmnet)),
#'                       cv_folds = 5,
#'                       silent = TRUE)
#' ew$weights
#' }
ensemble_weights <- function(y, X,
                             type = "average",
                             learners = NULL,
                             cv_folds = 5,
                             cv_subsamples = NULL,
                             cv_results = NULL,
                             custom_weights = NULL,
                             silent = FALSE) {
  # Data parameters
  nlearners <- if (!is.null(cv_results)) {
    ncol(cv_results$cv_resid)
  } else {
    length(learners)
  }#IFELSE
  ncustom <- n_custom(custom_weights)
  ntype <- length(type)

  # Single learner: trivial weight of 1, skip cross-validation
  if (nlearners == 1) {
    weights <- matrix(1, 1, ntype + ncustom)
    if (ncustom > 0)
      weights[, (ntype + 1):(ntype + ncustom)] <- custom_weights
    if (ncustom > 0 && is.null(colnames(custom_weights)))
      colnames(custom_weights) <- paste0("custom_", seq_len(ncustom))
    colnames(weights) <- c(type, colnames(custom_weights))
    return(list(weights = weights, cv_results = NULL))
  }#IF

  # Check whether cross-validation is needed for data-driven weights
  cv_stacking <- c("ols", "nnls", "nnls1", "singlebest")
  if (any(cv_stacking %in% type) && is.null(cv_results)) {
    cv_results <- crossval(y, X,
                           learners = learners,
                           cv_folds = cv_folds,
                           cv_subsamples = cv_subsamples,
                           silent = silent)
  }#IF

  # Compute weights for each ensemble type via dispatch
  weights <- matrix(0, nlearners, ntype + ncustom)
  for (k in seq_len(ntype)) {
    weights[, k] <- compute_w(type[k], nlearners,
                              cv_results, y)
  }#FOR

  # Append custom weights
  if (ncustom > 0) {
    weights[, (ntype + 1):(ntype + ncustom)] <- custom_weights
  }#IF
  if (ncustom > 0 && is.null(colnames(custom_weights))) {
    colnames(custom_weights) <- paste0("custom_", seq_len(ncustom))
  }#IF
  colnames(weights) <- c(type, colnames(custom_weights))

  list(weights = weights, cv_results = cv_results)
}#ENSEMBLE_WEIGHTS

# Weight-type dispatch: maps a type string to a weight vector.
compute_w <- function(type, nlearners, cv_results, y) {
  switch(type,
    average = rep(1 / nlearners, nlearners),
    nnls1   = compute_w_nnls1(nlearners, cv_results),
    nnls    = compute_w_nnls(cv_results, y),
    ols     = compute_w_ols(cv_results, y),
    singlebest = compute_w_singlebest(nlearners, cv_results))
}#COMPUTE_W

compute_w_nnls1 <- function(nlearners, cv_results) {
  # QP solver requires a positive-definite matrix. Empirical
  # cross-products of CV residuals can be PSD (not PD) when
  # learners are collinear. nearPD finds the nearest PD matrix.
  sq_resid <- Matrix::crossprod(cv_results$cv_resid)
  A <- cbind(matrix(1, nlearners, 1), diag(1, nlearners))
  r <- tryCatch(
    quadprog::solve.QP(Dmat = Matrix::nearPD(sq_resid)$mat,
                       dvec = matrix(0, nlearners, 1),
                       Amat = A,
                       bvec = c(1, rep(0, nlearners))),
    error = function(e) {
      warning("nnls1 weight optimization failed: ",
              conditionMessage(e),
              ". Falling back to equal weights.",
              call. = FALSE)
      NULL
    })
  if (!is.null(r)) r$solution
  else rep(1 / nlearners, nlearners)
}#COMPUTE_W_NNLS1

# Unconstrained non-negative least squares (Wolpert-style).
# Unlike nnls1, weights are NOT normalized to sum to 1.
compute_w_nnls <- function(cv_results, y) {
  cv_fitted <- as.numeric(y) - cv_results$cv_resid
  nnls::nnls(cv_fitted, y)$x
}#COMPUTE_W_NNLS

compute_w_ols <- function(cv_results, y) {
  cv_fitted <- as.numeric(y) - cv_results$cv_resid
  ols(y, cv_fitted, const = FALSE)$coef
}#COMPUTE_W_OLS

compute_w_singlebest <- function(nlearners, cv_results) {
  mdl_min <- which.min(Matrix::colMeans(cv_results$cv_resid^2))
  mdl_min <- (seq_len(nlearners))[mdl_min]
  w <- rep(0, nlearners)
  w[mdl_min] <- 1
  w
}#COMPUTE_W_SINGLEBEST
