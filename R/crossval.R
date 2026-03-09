#' Estimator of the Mean Squared Prediction Error using Cross-Validation.
#'
#' @family utilities
#'
#' @description Estimator of the mean squared prediction error of
#'     different learners using cross-validation.
#'
#' @details \code{crossval} estimates the mean squared prediction error
#'     (MSPE) of \eqn{J} base learners via \eqn{K}-fold
#'     cross-validation. It is the inner workhorse of the stacking
#'     machinery used by \code{\link{ensemble_weights}} to determine
#'     ensemble weights.
#'
#' Given a generic conditional expectation function \eqn{f_0(\cdot)}
#'     (e.g., \eqn{E[Y\vert X]}, \eqn{E[D\vert X]}), let
#'     \eqn{\{I_1, \ldots, I_K\}} be a \eqn{K}-fold partition of
#'     \eqn{\{1, \ldots, n\}} and let \eqn{\hat{f}_j^{(-k)}} denote
#'     learner \eqn{j} trained on all observations outside fold
#'     \eqn{I_k}. The out-of-sample residual for observation
#'     \eqn{i \in I_k} is
#'
#' \eqn{\hat{e}_{i,j} = y_i - \hat{f}_j^{(-k)}(X_i).}
#'
#' Since every observation belongs to exactly one fold, this yields a
#'     complete \eqn{n \times J} residual matrix. The cross-validated
#'     MSPE for learner \eqn{j} is
#'
#' \eqn{\widehat{\textrm{MSPE}}_j = n^{-1} \sum_{i=1}^{n} \hat{e}_{i,j}^2,}
#'
#' and the cross-validated \eqn{R^2} is
#'
#' \eqn{\hat{R}^2_j = 1 - \widehat{\textrm{MSPE}}_j \,/\, \hat{\sigma}^2_y,}
#'
#' where \eqn{\hat{\sigma}^2_y} is the sample variance of \eqn{y}.
#'
#' @inheritParams ddml_plm
#' @param y The outcome variable.
#' @param X A (sparse) matrix of predictive variables.
#' @param learners \code{learners} is a list of lists, each containing three
#'     named elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to variables in \code{X} that are passed to
#'             the base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{what}. Omission of \code{assign_X}
#'     results in inclusion of all predictive variables in \code{X}.
#' @param cv_folds Number of folds used for cross-validation.
#' @param cv_subsamples List of vectors with sample indices for
#'     cross-validation.
#' @param parallel An optional named list with parallel processing
#'     options. When \code{NULL} (the default), computation is
#'     sequential. Supported fields:
#'     \describe{
#'         \item{\code{cores}}{Number of cores to use.}
#'         \item{\code{export}}{Character vector of object names to
#'             export to parallel workers (for custom learners that
#'             reference global objects).}
#'         \item{\code{packages}}{Character vector of additional
#'             package names to load on workers (for custom learners
#'             that use packages not imported by \code{ddml}).}
#'     }
#'
#' @return \code{crossval} returns a list containing the following components:
#'     \describe{
#'         \item{\code{mspe}}{A vector of MSPE estimates,
#'             each corresponding to a base learner (in chronological
#'             order).}
#'         \item{\code{r2}}{A vector of cross-validated \eqn{R^2}
#'             values, each corresponding to a base learner (in
#'             chronological order).}
#'         \item{\code{cv_resid}}{A matrix of out-of-sample residuals,
#'             each column corresponding to a base learner (in
#'             chronological order).}
#'         \item{\code{cv_subsamples}}{Pass-through of
#'             \code{cv_subsamples}. See above.}
#'     }
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]
#'
#' # Compare ols, lasso, and ridge using 4-fold cross-validation
#' cv_res <- crossval(y, X,
#'                    learners = list(list(what = ols),
#'                                    list(what = mdl_glmnet),
#'                                    list(what = mdl_glmnet,
#'                                         args = list(alpha = 0))),
#'                    cv_folds = 4,
#'                    silent = TRUE)
#' cv_res$mspe
crossval <- function(y, X,
                     learners,
                     cv_folds = 10,
                     cluster_variable = seq_along(y),
                     cv_subsamples = NULL,
                     silent = FALSE,
                     parallel = NULL) {
  # Unpack parallel options
  p <- parse_parallel(parallel)
  num_cores <- p$num_cores
  parallel_export <- p$export
  parallel_packages <- p$packages

  # Normalize learner specs before parallel dispatch
  learners <- normalize_learners(learners)

  # Validate inputs
  validate_inputs(y = y, X = X, learners = learners,
                  cv_folds = cv_folds)

  # Data parameters
  nobs <- length(y)
  nlearners <- length(learners)

  # Get cv subsample tuple
  indx <- get_crossfit_indices(cluster_variable,
                               sample_folds = cv_folds,
                               subsamples = cv_subsamples)
  cv_subsamples <- indx$subsamples
  cv_folds <- length(cv_subsamples)
  nobs <- length(unlist(cv_subsamples)) # In case subsamples are user-provided

  # Define the computation function
  cv_fun <- function(x) {
    j <- ceiling(x / cv_folds)
    i <- x - cv_folds * (ceiling(x / cv_folds) - 1)
    fold_x <- cv_subsamples[[i]]
    crossval_compute(test_sample = fold_x,
                     learner = learners[[j]],
                     y, X)
  }#CV_FUN

  # Compute out-of-sample errors
  njobs <- cv_folds * nlearners
  cl <- NULL
  if (num_cores > 1) {
    cl <- tryCatch(
      setup_parallel_cluster(num_cores, parallel_export,
                             parallel_packages),
      error = function(e) {
        warning("Parallel setup failed: ",
                conditionMessage(e),
                ". Falling back to sequential.",
                call. = FALSE)
        NULL
      }
    )
    if (!is.null(cl))
      on.exit(parallel::stopCluster(cl), add = TRUE)
  }#IF

  if (silent) {
    op <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(op), add = TRUE)
  }#IF
  cv_res <- pbapply::pbsapply(seq_len(njobs), cv_fun,
                              cl = cl)

  # Compile residual matrix
  cv_resid <- unlist(cv_res)
  cv_resid <- matrix(cv_resid, nobs, nlearners)
  cv_resid <- cv_resid[order(unlist(cv_subsamples)), , drop = FALSE]

  # Compute MSPE and R-squared by learner
  mspe <- colMeans(cv_resid^2)
  y_var <- as.numeric(stats::var(y))
  r2 <- if (y_var > 0) 1 - mspe / y_var else rep(NA_real_, length(mspe))

  # Organize and return output
  output <- list(mspe = mspe, r2 = r2,
                 cv_resid = cv_resid,
                 cv_subsamples = cv_subsamples)
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
crossval_compute <- function(test_sample, learner,
                             y, X) {
  if (is.null(learner$assign_X)) learner$assign_X <- seq_len(ncol(X))

  mdl_fun <- list(what = learner$what, args = learner$args)
  assign_X <- learner$assign_X

  mdl_fun$args$y <- y[-test_sample]
  mdl_fun$args$X <- X[-test_sample, assign_X, drop = FALSE]

  mdl_fit <- tryCatch(
    do.call(do.call, mdl_fun),
    error = function(e) {
      stop("Learner fitting failed: ", conditionMessage(e),
           call. = FALSE)
    }
  )

  cv_fitted <- stats::predict(mdl_fit,
                              X[test_sample, assign_X,
                                drop = FALSE])
  if (!is.matrix(cv_fitted)) cv_fitted <- as.matrix(cv_fitted)
  cv_resid <- y[test_sample] - cv_fitted
  return(cv_resid)
}#CROSSVAL_COMPUTE
