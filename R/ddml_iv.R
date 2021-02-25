#' Compute LIE-conform DDML IV estimators.
#'
#' Compute LIE-conform DDML IV estimators.
#'
#' @param y A response vector.
#' @param D An endogeneous variable vector.
#' @param Z An instrumental variable vector or matrix.
#' @param X A control matrix.
#' @param models May take one of two forms, depending on whether a single model
#'     should be used for residualization, or whether an ensemble procedure
#'     should be employed.
#'     If a single model should be used, \code{models} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If an ensemble should be used, \code{models} is a list of lists, each
#'     containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The function used to trained the model. The
#'             function must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{assign_X} A vector of indices corresponding to features
#'             in \code{X} that should be used for training.}
#'         \item{\code{assign_Z} An optional vector of indices corresponding to
#'             instruments in \code{Z} that should be used for training.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}}
#'     }
#' @param models_FS Same as \code{models}. May be used to consider an seperate
#'     set of models to be used for the estimation of E[D|X,Z] (first stage).
#' @param ens_type A string indicating the type of ensemble. Multiple types may
#'     be passed in form of a vector of strings.
#' @param cv_folds The number for cross-validation folds.
#' @param sample_folds The number of split sample folds used for calculation of
#'     the out-of-sample predictions.
#' @param subsamples An optional list of vectors, each containing indices of
#'     a test-sample. If not used-provided, the split sample folds are randomly
#'     drawn.
#' @param setup_parallel An list containing details on the parallelization of
#'    \code{crossval}.
#' @param silent A boolean indicating whether current models and folds should be
#'     printed to the console.
#'
#' @return \code{ddml_iv} returns an object of S3 class
#'     \code{cddml_iv}.
#'
#' An object of class \code{ddml_iv}is a list containig the
#'     following components:
#' @return \code{crosspred} returns a list containig the following
#'     components:
#' \describe{
#' \item{\code{coef}}{A vector with the LIE-conform DDML IV coefficent in the
#'     first entry.}
#'     \item{\code{...}}{...}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export ddml_iv
ddml_iv <- function(y, D, Z, X = matrix(1, nobs(y)),
                    models,
                    models_FS = models,
                    ens_type = c("average"),
                    cv_folds = 5,
                    sample_folds = 2,
                    subsamples = NULL,
                    setup_parallel = list(type = 'dynamic', cores = 1),
                    silent = F) {
  # Data parameters
  nobs <- length(y)
  # Draw samples if not user-supplied
  if (is.null(subsamples)) {
    subsamples <- split(c(1:nobs),
                        sample(rep(c(1:sample_folds),
                                   ceiling(nobs / sample_folds)))[1:nobs])
  }#IF
  sample_folds <- length(subsamples)

  # Compute number of required model computations
  cmp_cv <- length(setdiff(ens_type, "average")) > 0
  nmodels <- length(models); nensb <- length(ens_type)
  if(!silent & cmp_cv)  {
    num_comp <- nmodels * cv_folds * sample_folds *
      (2 + nensb - ("average" %in% ens_type)) +
      ("average" %in% ens_type) * nmodels * sample_folds
    print(paste0("# required model computations: ", num_comp))
  }#IF

  # Compute estimates of E[D|X,Z]
  D_XZ_res <- crosspred(D, X, Z,
                        models_FS, ens_type, cv_folds,
                        sample_folds, subsamples, compute_is_predictions = T,
                        setup_parallel, silent)

  # Compute estimates of E[y|X]
  y_X_res <- crosspred(y, X, Z = NULL,
                       models, ens_type, cv_folds,
                       sample_folds, subsamples, compute_is_predictions = F,
                       setup_parallel, silent)

  # Check whether multiple ensembles are computed simultaneously
  sim_ens <- "list" %in% class(D_XZ_res$is_fitted[[1]]) # use construction

  # If a single ensemble is calculated, no loops are required.
  if (!sim_ens) {
    # Compute LIE-conform estimates of E[D|X]
    D_X_res <- crosspred(D_XZ_res$is_fitted, X, Z = NULL,
                         models, ens_type, cv_folds,
                         sample_folds, subsamples, compute_is_predictions = F,
                         setup_parallel, silent)

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - D_X_res$oos_fitted
    V_r <- D_XZ_res$oos_fitted - D_X_res$oos_fitted

    # Compute IV estimate with constructed variables
    iv_fit <- tsls(y_r, D_r, V_r)

    # Organize complementary ensemble output
    coef <- iv_fit$coef[1]
    weights <- list(D_XZ = D_XZ_res$weights,
                    D_X = D_X_res$weights,
                    y_X = y_X_res$weights)
    mspe <- list(D_XZ = D_XZ_res$mspe,
                 D_X = D_X_res$mspe,
                 y_X = y_X_res$mspe)
    anyiv_cv <- list(D_XZ = D_XZ_res$anyiv_cv,
                     D_X = D_X_res$anyiv_cv,
                     y_X = y_X_res$anyiv_cv)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (sim_ens) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    nensb <- length(ens_type)
    coef <- matrix(0, 1, nensb)
    mspe <- anyiv_cv <- iv_fit <- rep(list(1), nensb)
    weights <- replicate(3, array(0, dim = c(nmodels, nensb, sample_folds)),
                         simplify = F)
    weights[[1]] <- D_XZ_res$weights; weights[[3]] <- y_X_res$weights
    names(weights) <- c("D_XZ", "D_X", "y_X")
    for (j in 1:nensb) {
      # Compute LIE-conform estimates of E[D|X]
      D_X_res <- crosspred(D_XZ_res$is_fitted[[j]], X, Z = NULL,
                           models, ens_type[j], cv_folds,
                           sample_folds, subsamples, compute_is_predictions = F,
                           setup_parallel, silent)

      # Residualize
      y_r <- y - y_X_res$oos_fitted[, j]
      D_r <- D - D_X_res$oos_fitted
      V_r <- D_XZ_res$oos_fitted[, j] - D_X_res$oos_fitted

      # Compute IV estimate with constructed variables
      iv_fit_j <- tsls(y_r, D_r, V_r)

      # Organize complementary ensemble output
      coef[j] <- iv_fit_j$coef[1]
      mspe[[j]] <- list(D_XZ = D_XZ_res$mspe,
                        D_X = D_X_res$mspe,
                        y_X = y_X_res$mspe)
      anyiv_cv[[j]] <- list(D_XZ = D_XZ_res$anyiv_cv,
                            D_X = D_X_res$anyiv_cv,
                            y_X = y_X_res$anyiv_cv)
      iv_fit[[j]] <- iv_fit_j
      weights[[2]][, j, ] <- D_X_res$weights
    }#FOR
    # Name output appropriately by ensemble type
    names(mspe) <- names(anyiv_cv) <- names(iv_fit) <-
      ens_type
  }#IF

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights,
                   mspe = mspe, anyiv_cv = anyiv_cv,
                   models = models, models_FS = models_FS,
                   iv_fit = iv_fit,
                   subsamples = subsamples,
                   ens_type = ens_type,
                   nobs = nobs,
                   y = y,
                   D = D)

  # Amend class and return
  class(ddml_fit) <- c("ddml_iv")
  return(ddml_fit)
}#DDML_IV

#' Function for exporting ddml_iv fits.
#'
#' Function for exporting ddml_iv fits.
#'
#' @export export_ddml.ddml_iv
#' @export
export_ddml.ddml_iv <- function(obj, filename) {

  # Data parameters
  nensemble <- length(obj$ens_type)

  # Assign fold-id
  fid <- rep(0, obj$nobs)
  sample_folds <- length(obj$subsamples)
  for (k in 1:sample_folds) {
    fid[obj$subsamples[[k]]] <- k
  }#FOR

  # Check whether a) multiple ensembles, or b) single ensemble, or c) single
  #     model were used for orthogonalization.
  y_mat <- D_mat <- Z_mat <- matrix(0, obj$nobs, nensemble)
  if (nensemble > 1) {
    # Multiple ensembles
    # Compile variables
    for (e in 1:nensemble) {
      y_mat[, e] <- obj$iv_fit[[e]]$y
      D_mat[, e] <- obj$iv_fit[[e]]$X_[, 1]
      Z_mat[, e] <- obj$iv_fit[[e]]$Z_[, 1]
    }#FOR

    # Assign names
    cnames <- paste0("_", names(obj$iv_fit))
    colnames(y_mat) <- paste0("y", cnames)
    colnames(D_mat) <- paste0("D", cnames)
    colnames(Z_mat) <- paste0("Z", cnames)
  } else {
    # For single ensemble or single model, only the naming differs. First,
    #     compile variables.
    y_mat <- obj$iv_fit$y
    D_mat <- as.matrix(obj$iv_fit$X_[, 1])
    Z_mat <- as.matrix(obj$iv_fit$Z_[, 1])
    # Then assign names.
        calc_ensemble <- !("what" %in% names(obj$models)) |
      !("what" %in% names(obj$models_FS))
    if (calc_ensemble) {
      # Single ensemble
      cnames <- paste0("_", obj$ens_type)
      colnames(y_mat) <- paste0("y", cnames)
      colnames(D_mat) <- paste0("D", cnames)
      colnames(Z_mat) <- paste0("Z", cnames)

    } else {
      # Single model
      # Single ensemble
      cnames <- paste0("_", "mdl")
      colnames(y_mat) <- paste0("y", cnames)
      colnames(D_mat) <- paste0("D", cnames)
      colnames(Z_mat) <- paste0("Z", cnames)
    }#IFELSE
  }#IFELSE

  # Print orthogonalized variables to csv
  ids <- cbind(c(1:obj$nobs), fid); colnames(ids) <- c("id", "fid")
  v_mat <- cbind(y_mat, D_mat, Z_mat)
  write.table(cbind(ids, v_mat),
            paste0(filename, ".csv"),
            row.names = FALSE,
            sep = ",")
  # Print dict
  y_name <- colnames(obj$y); if (is.null(y_name)) y_name <- "y"
  D_name <- colnames(obj$D); if (is.null(D_name)) D_name <- "D"
  dict <- rbind(colnames(v_mat),
                c(rep(y_name, nensemble), rep(D_name, 2 * nensemble)),
                c(rep("yeq", nensemble), rep("deq", nensemble),
                  rep("zeq", nensemble)))
  write.table(dict, paste0(filename, "_dict.csv"),
            col.names = FALSE, row.names = FALSE,
            sep = ",")
}#EXPORT_DDML.DDML_IV
