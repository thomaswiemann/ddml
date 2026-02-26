#' Estimator of the Mean Squared Prediction Error using Cross-Validation.
#'
#' @family utilities
#'
#' @description Estimator of the mean squared prediction error of
#'     different learners using cross-validation.
#'
#' @inheritParams ddml_plm
#' @param y The outcome variable.
#' @param X A (sparse) matrix of predictive variables.
#' @param Z Optional additional (sparse) matrix of predictive variables.
#' @param learners \code{learners} is a list of lists, each containing four
#'     named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to variables in \code{Z} that are passed to the
#'             base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all predictive variables in \code{X} (and/or
#'     \code{Z}).
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
#'             each corresponding to a base learners (in chronological order).}
#'         \item{\code{oos_resid}}{A matrix of out-of-sample prediction errors,
#'             each column corresponding to a base learners (in chronological
#'             order).}
#'         \item{\code{cv_subsamples}}{Pass-through of \code{cv_subsamples}.
#'             See above.}
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
crossval <- function(y, X, Z = NULL,
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

  # Normalize learner specs: resolve fun/what before parallel dispatch
  learners <- normalize_learners(learners)

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
                     y, X, Z)
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
  oos_resid <- unlist(cv_res)
  oos_resid <- matrix(oos_resid, nobs, nlearners)
  oos_resid <- oos_resid[order(unlist(cv_subsamples)), , drop = FALSE]

  # Compute MSPE and R-squared by learner
  mspe <- colMeans(oos_resid^2)
  y_var <- as.numeric(stats::var(y))
  r2 <- if (y_var > 0) 1 - mspe / y_var else rep(NA_real_,
                                                   length(mspe))

  # Organize and return output
  output <- list(mspe = mspe, r2 = r2,
                 oos_resid = oos_resid,
                 cv_subsamples = cv_subsamples)
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
crossval_compute <- function(test_sample, learner,
                             y, X, Z = NULL) {
  if (is.null(learner$assign_X)) learner$assign_X <- seq_len(ncol(X))
  if (is.null(learner$assign_Z) && !is.null(Z))
    learner$assign_Z <- seq_len(ncol(Z))

  mdl_fun <- list(what = learner$what, args = learner$args)
  assign_X <- learner$assign_X
  assign_Z <- learner$assign_Z

  mdl_fun$args$y <- y[-test_sample]
  mdl_fun$args$X <- cbind(X[-test_sample, assign_X, drop = FALSE],
                          Z[-test_sample, assign_Z, drop = FALSE])

  mdl_fit <- tryCatch(
    do.call(do.call, mdl_fun),
    error = function(e) {
      stop("Learner fitting failed: ", conditionMessage(e),
           call. = FALSE)
    }
  )

  oos_fitted <- stats::predict(mdl_fit,
                               cbind(X[test_sample, assign_X,
                                       drop = FALSE],
                                     Z[test_sample, assign_Z,
                                       drop = FALSE]))
  if (!is.matrix(oos_fitted)) oos_fitted <- as.matrix(oos_fitted)
  oos_resid <- y[test_sample] - oos_fitted
  return(oos_resid)
}#CROSSVAL_COMPUTE
