#' Estimator for the Flexible Partially Linear IV Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()], [ddml::coef.ddml()],
#'     [ddml::confint.ddml()], [ddml::tidy.ddml()],
#'     [ddml::glance.ddml()], [ddml::diagnostics()],
#'     [AER::ivreg()]
#'
#' @description Estimator for the flexible partially linear IV model.
#'
#' @details \code{ddml_fpliv} provides a double/debiased machine learning
#'     estimator for the parameter of interest \eqn{\theta_0} in the partially
#'     linear IV model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, Z, U)} is a random vector such that
#'     \eqn{E[U\vert X, Z] = 0} and \eqn{E[Var(E[D\vert X, Z]\vert X)] \neq 0},
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' @inheritParams ddml_pliv
#' @param Z A (sparse) matrix of instruments.
#' @param learners_DXZ,learners_DX Optional arguments to allow for different
#'     estimators of \eqn{E[D \vert X, Z]}, \eqn{E[D \vert X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DXZ,custom_ensemble_weights_DX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DXZ},\code{learners_DX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DXZ},\code{custom_ensemble_weights_DX} must
#'     have the same number of columns.
#' @param enforce_LIE Indicator equal to 1 if the law of iterated expectations
#'     is enforced in the first stage.
#' @param fitted An optional named list of per-equation cross-fitted
#'     predictions, typically obtained via \code{fit$fitted}. Not
#'     supported when \code{enforce_LIE = TRUE}. See
#'     \code{\link{ddml_plm}} for details.
#' @param save_crossval Logical; store inner cross-validation
#'     residuals for exact weight recomputation on pass-through.
#'     See \code{\link{ddml_plm}} for details.
#'
#' @return \code{ddml_fpliv} returns an object of S3 class
#'     \code{ddml_fpliv}. An object of class \code{ddml_fpliv} is a list
#'     containing the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the \eqn{\theta_0} estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{iv_fit}}{Object of class \code{ivreg} from the IV
#'             regression of \eqn{Y - \hat{E}[Y\vert X]} on
#'             \eqn{D - \hat{E}[D\vert X]} using
#'             \eqn{\hat{E}[D\vert X,Z] - \hat{E}[D\vert X]} as the instrument.}
#'         \item{\code{learners},\code{learners_DX},\code{learners_DXZ},
#'             \code{cluster_variable},\code{subsamples},
#'             \code{cv_subsamples},\code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
#' Ahrens A, Hansen C B, Schaffer M E, Wiemann T (2024). "Model Averaging and 
#'     Double Machine Learning." Journal of Applied Econometrics, 40(3): 249-269.
#'
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B, Newey W,
#'     Robins J (2018). "Double/debiased machine learning for treatment and
#'     structural parameters." The Econometrics Journal, 21(1), C1-C68.
#'
#' Wolpert D H (1992). "Stacked generalization." Neural Networks, 5(2), 241-259.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex", drop = FALSE]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear IV model using a single base learner: Ridge.
#' fpliv_fit <- ddml_fpliv(y, D, Z, X,
#'                         learners = list(what = mdl_glmnet,
#'                                         args = list(alpha = 0)),
#'                         sample_folds = 2,
#'                         silent = TRUE)
#' summary(fpliv_fit)
ddml_fpliv <- function(y, D, Z, X,
                       learners,
                       learners_DXZ = learners,
                       learners_DX = learners,
                       sample_folds = 10,
                       ensemble_type = "nnls",
                       shortstack = FALSE,
                       cv_folds = 10,
                       enforce_LIE = TRUE,
                       custom_ensemble_weights = NULL,
                       custom_ensemble_weights_DXZ = custom_ensemble_weights,
                       custom_ensemble_weights_DX = custom_ensemble_weights,
                       cluster_variable = seq_along(y),
                       silent = FALSE,
                       parallel = NULL,
                       fitted = NULL,
                       splits = NULL,
                       save_crossval = TRUE,
                       ...) {
  # Validate inputs
  validate_inputs(y = y, D = D, X = X, Z = Z, learners = learners,
                  sample_folds = sample_folds, cv_folds = cv_folds,
                  ensemble_type = ensemble_type)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DXZ, learners_DXZ)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  if (!is.null(fitted) && enforce_LIE) {
    stop(paste("Pass-through via the 'fitted' argument is not",
               "currently supported when enforce_LIE = TRUE."))
  }#IF

  # Data parameters
  nobs <- length(y)
  nensb_raw <- length(ensemble_type) # number of ensembles w/o custom weights

  # Check for multivariate endogenous variables
  D <- as.matrix(D)
  nD <- ncol(D)

  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- !shortstack &
    any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (class(learners[[1]]) != "function")

  # Normalize deprecated split arguments into a single splits object
  splits <- normalize_splits(splits = splits, ...)
  validate_fitted_splits_pair(fitted, splits, w_cv)

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = if (w_cv) cv_folds,
                             subsamples = splits$subsamples,
                             cv_subsamples = splits$cv_subsamples)
  check_subsamples(indxs$subsamples, NULL, stratify = FALSE)

  # Estimation start
  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }
  info_msg("ddml_fpliv: estimating (", mode_str, ")",
           silent = silent)

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X, Z = NULL,
                     learners = learners, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     compute_insample_predictions = FALSE,
                     silent = silent, label = "E[Y|X]",
                     parallel = parallel,
                     fitted = fitted$y_X)

  # Compute estimates of E[D|X,Z]. Also calculate in-sample predictions when
  #     the LIE is enforced.
  D_XZ_res_list <- compute_CEF_list(
    D, X, Z = Z, learners = learners_DXZ,
    ensemble_type = ensemble_type,
    shortstack = shortstack,
    custom_ensemble_weights = custom_ensemble_weights_DXZ,
    subsamples = indxs$subsamples,
    cv_subsamples = indxs$cv_subsamples,
    compute_insample_predictions = enforce_LIE,
    silent = silent,
    label_prefix = "E[D", label_suffix = "|X,Z]",
    parallel = parallel,
    fitted = fitted$D_XZ)

  # When the LIE is not enforced, estimating E[D|X] is straightforward.
  if (!enforce_LIE) {
    D_X_res_list <- compute_CEF_list(
      D, X, Z = NULL, learners = learners_DX,
      ensemble_type = ensemble_type,
      shortstack = shortstack,
      custom_ensemble_weights = custom_ensemble_weights_DX,
      subsamples = indxs$subsamples,
      cv_subsamples = indxs$cv_subsamples,
      compute_insample_predictions = FALSE,
      silent = silent,
      label_prefix = "E[D", label_suffix = "|X]",
      parallel = parallel)
  }#IF

  # Update ensemble type to account for (optional) custom weights
  ensb_info <- update_ensemble_info(y_X_res$weights)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb
  multiple_ensembles <- ensb_info$multiple_ensembles

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {
    # Check whether the law of iterated expectations (LIE) should be enforced.
    #     When the LIE is enforced (recommended), the estimates of E[D|X,Z] are
    #     used for the calculation of the estimates of E[D|X].
    if (enforce_LIE) {
      D_X_res_list <- list()
      for (k in seq_len(nD)) {
        D_X_res_list[[k]] <-
          get_CEF(D_XZ_res_list[[k]]$is_fitted, X, Z = NULL,
                  learners = learners_DX,
                  ensemble_type = ensemble_type,
                  shortstack = shortstack,
                  subsamples = indxs$subsamples,
                  cv_subsamples = indxs$cv_subsamples,
                  compute_insample_predictions = FALSE,
                  silent = silent,
                  label = paste0("E[D", k, "|X]"),
                  shortstack_y = D_XZ_res_list[[k]]$oos_fitted,
                  parallel = parallel)
      }#FOR
    }#IFELSE

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - get_oosfitted(D_X_res_list)
    V_r <- get_oosfitted(D_XZ_res_list) - get_oosfitted(D_X_res_list)

    # Compute IV estimate with constructed variables
    iv_fit <- AER::ivreg(y_r ~ D_r | V_r)

    # Organize complementary ensemble output
    coef_vec <- stats::coef(iv_fit)[-1]
    coef <- matrix(coef_vec, nrow = nD, ncol = 1)
    colnames(coef) <- ensemble_type

    # Compute scores and Jacobian
    D_r_mat <- as.matrix(D_r)
    V_r_mat <- as.matrix(V_r)
    D_hat <- stats::lm.fit(V_r_mat, D_r_mat)$fitted.values
    e <- as.vector(y_r - D_r_mat %*% coef)
    scores <- list(D_hat * e)
    J_list <- list(-crossprod(D_hat, D_r_mat) / nobs)
    coef_names <- names(coef_vec)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML IV estimate for each.
    coef <- matrix(0, nD, nensb)
    iv_fit <- rep(list(1), nensb)
    scores <- vector("list", nensb)
    J_list <- vector("list", nensb)
    nlearners <- length(learners)
    nlearners_DX <- length(learners_DX); nlearners_DXZ <- length(learners_DXZ)
    # Assign names for more legible output
    colnames(coef) <- names(iv_fit) <- ensemble_type

    # Create intermediate weight matrices for enforce_LIE = TRUE
    if (enforce_LIE) {
      # weights
      weights_DX <- array(0, dim = c(nlearners_DX, nensb, sample_folds))
      dimnames(weights_DX) <- dimnames(y_X_res$weights)
      weights_DX <- replicate(nD, weights_DX, simplify = FALSE)
    }#IF

    # Compute coefficients for each ensemble
    for (j in seq_len(nensb)) {
      # When the LIE is enforced, compute LIE-conform estimates of E[D|X].
      #     Otherwise use the previously calculated estimates of E[D|X].
      if (enforce_LIE) {
        D_X_res_list <- list()
        for (k in seq_len(nD)) {
          label_jk <- paste0("E[D", k, "|X] (",
                             ensemble_type[j], ")")

          # Check whether j is a custom ensemble specification. Necessary to
          #     assign correct corresponding custom_weights vector.
          if (j <= nensb_raw) { # j is not a custom specification
            D_X_res_list[[k]] <-
              get_CEF(D_XZ_res_list[[k]]$is_fitted[[j]], X, Z = NULL,
                      learners = learners_DX,
                      ensemble_type = ensemble_type[j],
                      shortstack = shortstack,
                      subsamples = indxs$subsamples,
                      cv_subsamples = indxs$cv_subsamples,
                      compute_insample_predictions = FALSE,
                      silent = silent,
                      label = label_jk,
                      shortstack_y = D_XZ_res_list[[k]]$oos_fitted[, j],
                      parallel = parallel)
          } else { # j is a custom specification
            D_X_res_list[[k]] <-
              get_CEF(D_XZ_res_list[[k]]$is_fitted[[j]], X, Z = NULL,
                      learners = learners_DX,
                      ensemble_type = "average",
                      shortstack = shortstack,
                      custom_ensemble_weights =
                        custom_ensemble_weights_DX[, j - nensb_raw, drop = FALSE],
                      subsamples = indxs$subsamples,
                      cv_subsamples = indxs$cv_subsamples,
                      compute_insample_predictions = FALSE,
                      silent = silent,
                      label = label_jk,
                      shortstack_y = D_XZ_res_list[[k]]$oos_fitted[, j],
                      parallel = parallel)
            # Remove "average" oos_fitted and weights
            D_X_res_list[[k]]$oos_fitted <- D_X_res_list[[k]]$oos_fitted[, -1]
            D_X_res_list[[k]]$weights <- D_X_res_list[[k]]$weights[, -1, ,
                                                                   drop = FALSE]
        }#IFELSE
      }#FOR
      }#IF

      # Residualize
      if (enforce_LIE) {
        D_r <- D - get_oosfitted(D_X_res_list)
        V_r <- get_oosfitted(D_XZ_res_list, j) - get_oosfitted(D_X_res_list)
      } else {
        D_r <- D - get_oosfitted(D_X_res_list, j)
        V_r <- get_oosfitted(D_XZ_res_list, j) - get_oosfitted(D_X_res_list, j)
      }#IFELSE

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute IV estimate with constructed variables
      iv_fit_j <- AER::ivreg(y_r ~ D_r | V_r)

      # Organize complementary ensemble output
      coef[, j] <- stats::coef(iv_fit_j)[-1]
      iv_fit[[j]] <- iv_fit_j

      # Compute scores and Jacobian
      D_r_mat <- as.matrix(D_r)
      V_r_mat <- as.matrix(V_r)
      D_hat <- stats::lm.fit(V_r_mat, D_r_mat)$fitted.values
      e_j <- as.vector(y_r - D_r_mat %*% coef[, j])
      scores[[j]] <- D_hat * e_j
      J_list[[j]] <- -crossprod(D_hat, D_r_mat) / nobs

      if (enforce_LIE) {
        for (k in seq_len(nD)) weights_DX[[k]][, j, ] <- D_X_res_list[[k]]$weights
      }#IF
    }#FOR
    rownames(coef) <- names(iv_fit[[1]]$coefficients)[-1]
    coef_names <- rownames(coef)
  }#IF

  # Ensemble metrics
  weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  r2 <- list(y_X = y_X_res$r2)
  for (k in seq_len(nD)) {
    if (enforce_LIE & multiple_ensembles) {
      weights[[paste0("D", k, "_X")]] <- weights_DX[[k]]
    } else {
      weights[[paste0("D", k, "_X")]] <-
        D_X_res_list[[k]]$weights
    }#IFELSE
  }#FOR
  for (k in seq_len(nD)) {
    eq <- paste0("D", k, "_XZ")
    weights[[eq]] <- D_XZ_res_list[[k]]$weights
    mspe[[eq]] <- D_XZ_res_list[[k]]$mspe
    r2[[eq]] <- D_XZ_res_list[[k]]$r2
  }#FOR

  fitted_export <- list(
    y_X = build_fitted_entry(y_X_res, save_crossval),
    D_XZ = build_fitted_from_list(D_XZ_res_list, save_crossval)
  )
  splits_export <- list(subsamples = indxs$subsamples,
                        cv_subsamples = indxs$cv_subsamples)

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_DX = learners_DX,
                   iv_fit = iv_fit,
                   cluster_variable = cluster_variable,
                   subsamples = indxs$subsamples,
                   cv_subsamples = indxs$cv_subsamples,
                   ensemble_type = ensemble_type,
                   enforce_LIE = enforce_LIE,
                   coefficients = coef,
                   scores = scores,
                   J = J_list,
                   coef_names = coef_names,
                   nobs = nobs,
                   sample_folds = sample_folds,
                   cv_folds = if (shortstack) NULL
                     else cv_folds,
                   shortstack = shortstack,
                   fitted = fitted_export,
                   splits = splits_export,
                   r2 = r2)

  # Print estimation completion
  elapsed <- round(proc.time()[3] - t0, 1)
  info_msg("ddml_fpliv: completed in ", elapsed, "s",
           silent = silent)

  # Amend class and return
  class(ddml_fit) <- c("ddml_fpliv", "ddml")
  return(ddml_fit)
}#DDML_FPLIV
