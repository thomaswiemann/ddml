#' Compute out-of-sample residuals via cross-validation.
#'
#' Compute out-of-sample residuals via cross-validation.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param Z An optional instrument matrix.
#' @param models A list of lists, each containing three named elements:
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
#' @param cv_subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
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
#' @return \code{crossval} a list containig the following components:
#' \describe{
#' \item{\code{oos_resid}}{A matrix of out-of-sample residuals, each column
#'     corresponding to a model in \code{models}.}
#' \item{\code{cv_Z}}{A vector of model indices corresponding to models that
#'     selected at least one instrument in each cross-validation fold.}
#' }
#'
#' @examples
#'
#' # Simple example w/o instruments ============================================
#' # Simulate small dataset
#' X <- cbind(1, matrix(rnorm(100*99), 100, 99)) # Simulate features
#' y <- X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
#' # Define arguments
#' models <- list(list(fun = rlasso,
#'                     args = list(include = NULL,
#'                                 iter_resid = 1, d = 5),
#'                     assign_X = c(1:ncol(X))), # rlasso
#'                list(fun = ols,
#'                     args = list(),
#'                     assign_X = c(1:ncol(X))), # ols w/ all features
#'                list(fun = ols,
#'                     args = list(),
#'                     assign_X = c(1:50)) # ols w/ some arbitrary features
#' # Compute cross-validation residuals
#' cv_res <- crossval(y, X,
#'                    models = models,
#'                    cv_folds = 3,
#'                    silent = T)
#' # Compute MSPE by model
#' colMeans(cv_res$oos_resid)
#'
#' @export crossval
crossval <- function(y, X, Z = NULL,
                     models,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     setup_parallel = list(type = 'dynamic',
                                           cores = 1),
                     silent = F) {
  # Data parameters
  nobs <- length(y)
  nmodels <- length(models)

  # Create CV fold tuple
  if (is.null(cv_subsamples)) {
    cv_subsamples <- split(c(1:nobs), sample(rep(c(1:cv_folds),
                                              ceiling(nobs / cv_folds)))[1:nobs])
  }#IF
  cv_folds <- length(cv_subsamples)
  nobs <- length(unlist(cv_subsamples))

  # Run cross-validation depending on parallelization specification
  if (setup_parallel$cores == 1) {
    # If only one core is available, don't parallelize
    cv_res <- sapply(c(1:(cv_folds * nmodels)), function(x) {
      # Select model and cv-fold for this job
      j <- ceiling(x / cv_folds) # jth model
      i <- x - cv_folds * (ceiling(x / cv_folds) - 1) # ith CV fold
      fold_x <- cv_subsamples[[i]]
      # Print fold and lambda
      if(!silent) print(paste0(paste0("model = ", j, paste0(": Fold ", i))))
      # Compute model for this fold
      crossval_compute(test_sample = fold_x,
                       model = models[[j]],
                       y, X, Z)
    })#SAPPLY
  } else if (setup_parallel$cores > 1) {
    # If two or more cores are available: parallelize. Begin by setting up the
    #    clusters. Load data to each cluster. Then load the ddml package.
    cl <- parallel::makeCluster(setup_parallel$cores, outfile="")
    parallel::clusterExport(cl,
                            c("y", "X", "Z",
                              "models", "cv_folds",
                              "nobs", "nmodels", "cv_folds",
                              "silent"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(ddml))

    # Use different parallelization approaches depending on parallel$type.
    if (setup_parallel$type == 'static') {
      # For static job scheduling, use parSapply.
      cv_res <- parallel::parSapply(cl, c(1:(cv_folds * nmodels)), function(x){
        # Select model and cv-fold for this job
        j <- ceiling(x / cv_folds) # jth model
        i <- x - cv_folds * (ceiling(x / cv_folds) - 1) # ith CV fold
        fold_x <- cv_subsamples[[i]]
        # Print fold and lambda
        if(!silent) print(paste0(paste0("model = ", j, paste0(": Fold ", i))))
        # Compute model for this fold
        ddml::crossval_compute(test_sample = fold_x,
                         model = models[[j]],
                         y, X, Z)
      })#PARSAPPLY
    } else if (setup_parallel$type == 'dynamic') {
      # For dynamic job scheduling, use foreach
      doParallel::registerDoParallel(cl, cores = setup_parallel$cores)
      `%dopar%` <- foreach::`%dopar%` # workaround
      cv_res <- foreach::foreach(x = 1:(cv_folds * nmodels),
                        .combine = "rbind") %dopar% {
        # Select model and cv-fold for this job
        j <- ceiling(x / cv_folds)
        i <- x - cv_folds * (ceiling(x / cv_folds) - 1)
        fold_x <- cv_subsamples[[i]]
        # Print fold and lambda
        if(!silent) print(paste0(paste0("model = ", j, paste0("; Fold ", i))))
        # Compute model for this fold
        ddml::crossval_compute(test_sample = fold_x,
                         model = models[[j]],
                         y, X, Z)
      }#FOREACH
    }#IFESLE
    # Terminate all parrallel clusters
    parallel::stopCluster(cl)
  }#IFELSE

  # Compile residual matrix and instrument variable selection list. This is a
  #     bit messy due as the crossvalidation procedure returns a list of lists.
  #     The solution differs by a transpose depending on whether foreach was
  #     used to compute the results.

  if (setup_parallel$cores == 1 || setup_parallel$type == "static"){
    cv_res <- t(cv_res)
  }#IF

  oos_resid <- unlist(sapply(c(1:(cv_folds * nmodels)),
                             function(x){cv_res[[x, 1]]}))
  cv_Z <- unlist(sapply(c(1:(cv_folds * nmodels)), function(x){cv_res[[x, 2]]}))
  # Then compile to final results in matrix format
  oos_resid <- matrix(oos_resid, nobs, nmodels)
  cv_Z <- which(apply(matrix(cv_Z, cv_folds, nmodels), 2, all))

  # Organize and return output
  output <- list(oos_resid = oos_resid,
                 cv_Z = cv_Z)
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
#' Compute results for a single cross-validation instance.
#'
#' Compute results for a single cross-validation instance.
#'
#' @export crossval_compute
crossval_compute <- function(test_sample, model,
                             y, X, Z = NULL) {
  # Check whether X, Z assignment has been specified. If not, include all.
  if (is.null(model$assign_X)) model$assign_X <- c(1:ncol(X))
  if (is.null(model$assign_Z) & !is.null(Z)) model$assign_Z <- c(1:ncol(Z))

  # Extract model arguments
  mdl_fun <- list(what = model$fun, args = model$args)
  assign_X <- model$assign_X
  assign_Z <- model$assign_Z

  # Compute model for this fold
  #     Note: this is effectively copying the data -- fix needed.
  mdl_fun$args$y <- y[-test_sample];
  mdl_fun$args$X <- cbind(X[-test_sample, assign_X, drop = F],
                          Z[-test_sample, assign_Z, drop = F])
  mdl_fit <- do.call(do.call, mdl_fun)

  # Compute out of sample residuals
  oos_fitted <- predict(mdl_fit, cbind(X[test_sample, assign_X, drop = F],
                                       Z[test_sample, assign_Z, drop = F]))
  oos_resid <- y[test_sample] - as(oos_fitted, "matrix")

  # Check whether instruments were selected (optional) and return fold
  #     results.
  if (!is.null(Z)) {
    index_iv <- (length(assign_X) + 1):length(c(assign_X, assign_Z))
    cv_Z <- any_iv(obj = mdl_fit,
                   index_iv = index_iv,
                   names_iv = colnames(Z[1:2, assign_Z, drop = F]))
  } else {
    cv_Z <- TRUE
  }#IFELSE
  # Return residuals and cv_Z
  return(list(oos_resid = oos_resid,
              cv_Z = cv_Z))
}#CROSSVAL_COMPUTE
