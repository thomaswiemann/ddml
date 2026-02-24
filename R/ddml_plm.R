#' Estimator for the Partially Linear Model.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()]
#'
#' @description Estimator for the partially linear model.
#'
#' @details \code{ddml_plm} provides a double/debiased machine learning
#'     estimator for the parameter of interest \eqn{\theta_0} in the partially
#'     linear model given by
#'
#' \eqn{Y = \theta_0D + g_0(X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{E[Cov(U, D\vert X)] = 0} and \eqn{E[Var(D\vert X)] \neq 0}, and
#'     \eqn{g_0} is an unknown nuisance function.
#'
#' @param y The outcome variable.
#' @param D A matrix of endogenous variables.
#' @param X A (sparse) matrix of control variables.
#' @param learners May take one of two forms, depending on whether a single
#'     learner or stacking with multiple learners is used for estimation of the
#'     conditional expectation functions.
#'     If a single learner is used, \code{learners} is a list with two named
#'     elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is a list of
#'     lists, each containing four named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to control variables in \code{X} that are passed to
#'             the base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} results in inclusion of
#'     all variables in \code{X}.
#' @param learners_DX Optional argument to allow for different estimators of
#'     \eqn{E[D|X]}. Setup is identical to \code{learners}.
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into final
#'     estimate of the conditional expectation functions. Possible values are:
#'     \itemize{
#'         \item{\code{"nnls"} Non-negative least squares.}
#'         \item{\code{"nnls1"} Non-negative least squares with the constraint
#'             that all weights sum to one.}
#'         \item{\code{"singlebest"} Select base learner with minimum MSPE.}
#'         \item{\code{"ols"} Ordinary least squares.}
#'         \item{\code{"average"} Simple average over base learners.}
#'     }
#'     Multiple ensemble types may be passed as a vector of strings.
#' @param shortstack Boolean to use short-stacking.
#' @param cv_folds Number of folds used for cross-validation in ensemble
#'     construction.
#' @param custom_ensemble_weights A numerical matrix with user-specified
#'     ensemble weights. Each column corresponds to a custom ensemble
#'     specification, each row corresponds to a base learner in \code{learners}
#'     (in chronological order). Optional column names are used to name the
#'     estimation results corresponding the custom ensemble specification.
#' @param custom_ensemble_weights_DX Optional argument to allow for different
#'     custom ensemble weights for \code{learners_DX}. Setup is identical to
#'     \code{custom_ensemble_weights}. Note: \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DX} must have the same number of columns.
#' @param cluster_variable A vector of cluster indices.
#' @param subsamples List of vectors with sample indices for cross-fitting.
#' @param cv_subsamples List of lists, each corresponding to a subsample
#'     containing vectors with subsample indices for cross-validation.
#' @param cv_subsamples_list Deprecated; use \code{cv_subsamples} instead.
#' @param silent Boolean to silence estimation updates.
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
#' @return \code{ddml_plm} returns an object of S3 class
#'     \code{ddml_plm}. An object of class \code{ddml_plm} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{coef}}{A vector with the \eqn{\theta_0} estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{ols_fit}}{Object of class \code{lm} from the second
#'             stage regression of \eqn{Y - \hat{E}[Y|X]} on
#'             \eqn{D - \hat{E}[D|X]}.}
#'         \item{\code{learners},\code{learners_DX},\code{cluster_variable},
#'             \code{subsamples}, \code{cv_subsamples},
#'             \code{ensemble_type}}{Pass-through of selected user-provided
#'             arguments. See above.}
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
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the partially linear model using a single base learner, ridge.
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
#'
#' # Estimate the partially linear model using short-stacking with base learners
#' #     ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' plm_fit <- ddml_plm(y, D, X,
#'                     learners = list(list(fun = ols),
#'                                     list(fun = mdl_glmnet),
#'                                     list(fun = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     custom_ensemble_weights = weights_everylearner,
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(plm_fit)
ddml_plm <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     cluster_variable = seq_along(y),
                     subsamples = NULL,
                     cv_subsamples = NULL,
                     cv_subsamples_list = NULL,
                     silent = FALSE,
                     parallel = NULL) {
  # Validate inputs
  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds, cv_folds = cv_folds,
                  ensemble_type = ensemble_type)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  # Backward compatibility for renamed parameter
  if (!is.null(cv_subsamples_list)) {
    if (!is.null(cv_subsamples))
      stop("Specify cv_subsamples or cv_subsamples_list, not both.")
    message("Note: cv_subsamples_list has been renamed to cv_subsamples.")
    cv_subsamples <- cv_subsamples_list
  }#IF

  # Data parameters
  nobs <- length(y)

  # Check for multivariate endogenous variables
  D <- as.matrix(D)
  nD <- ncol(D)

  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- !shortstack &
    any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (class(learners[[1]]) != "function" |
     class(learners_DX[[1]]) != "function")

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = if (w_cv) cv_folds,
                             subsamples = subsamples,
                             cv_subsamples = cv_subsamples)
  check_subsamples(indxs$subsamples, NULL, stratify = FALSE)

  # Estimation start
  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }
  info_msg("ddml_plm: estimating (", mode_str, ")",
           silent = silent)

  # Compute estimates of E[y|X]
  y_X_res <- get_CEF(y, X,
                     learners = learners,
                     ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = "E[Y|X]",
                     parallel = parallel)

  # Compute estimates of E[D|X], loop through endogenous variables
  D_X_res_list <- compute_CEF_list(
    D, X, learners = learners_DX,
    ensemble_type = ensemble_type,
    shortstack = shortstack,
    custom_ensemble_weights = custom_ensemble_weights_DX,
    subsamples = indxs$subsamples,
    cv_subsamples = indxs$cv_subsamples,
    silent = silent,
    label_prefix = "E[D", label_suffix = "|X]",
    parallel = parallel)

  # Update ensemble type to account for (optional) custom weights
  ensb_info <- update_ensemble_info(y_X_res$weights)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb
  multiple_ensembles <- ensb_info$multiple_ensembles

  # If a single ensemble is calculated, no loops are required.
  if (!multiple_ensembles) {

    # Residualize
    y_r <- y - y_X_res$oos_fitted
    D_r <- D - get_oosfitted(D_X_res_list)

    # Compute OLS estimate with constructed variables
    ols_fit <- stats::lm(y_r ~ D_r)

    # Organize complementary ensemble output
    coef <- stats::coef(ols_fit)[-1]

    # Compute scores and Jacobian
    D_r_mat <- as.matrix(D_r)
    e <- as.vector(y_r - D_r_mat %*% coef)
    scores <- list(D_r_mat * e)
    J_list <- list(-crossprod(D_r_mat) / nobs)
    coef_names <- names(coef)
  }#IF

  # If multiple ensembles are calculated, iterate over each type.
  if (multiple_ensembles) {
    # Iterate over ensemble type. Compute DDML estimate for each.
    coef <- matrix(0, nD, nensb)
    mspe <- ols_fit <- rep(list(1), nensb)
    scores <- vector("list", nensb)
    J_list <- vector("list", nensb)
    nlearners <- length(learners); nlearners_DX <- length(learners_DX)

    # Compute coefficients for each ensemble
    for (j in seq_len(nensb)) {
      # Residualize
      D_r <- D - get_oosfitted(D_X_res_list, j)

      # Residualize y
      y_r <- y - y_X_res$oos_fitted[, j]

      # Compute OLS estimate with constructed variables
      ols_fit_j <- stats::lm(y_r ~ D_r)

      # Organize complementary ensemble output
      coef[, j] <- stats::coef(ols_fit_j)[-1]
      ols_fit[[j]] <- ols_fit_j

      # Compute scores and Jacobian
      D_r_mat <- as.matrix(D_r)
      e_j <- as.vector(y_r - D_r_mat %*% coef[, j])
      scores[[j]] <- D_r_mat * e_j
      J_list[[j]] <- -crossprod(D_r_mat) / nobs
    }#FOR

    # Assign names for more legible output
    colnames(coef) <- names(ols_fit) <- dimnames(y_X_res$weights)[[2]]
    rownames(coef) <- names(ols_fit_j$coefficients)[-1]
    coef_names <- rownames(coef)
  }#IF

  # Ensemble metrics and per-learner residuals
  weights <- list(y_X = y_X_res$weights)
  mspe <- list(y_X = y_X_res$mspe)
  r2 <- list(y_X = y_X_res$r2)
  oos_resid_bylearner <- list(
    y_X = y_X_res$oos_resid_bylearner)
  for (k in seq_len(nD)){
    weights[[paste0("D", k, "_X")]] <-
      D_X_res_list[[k]]$weights
    mspe[[paste0("D", k, "_X")]] <-
      D_X_res_list[[k]]$mspe
    r2[[paste0("D", k, "_X")]] <-
      D_X_res_list[[k]]$r2
    oos_resid_bylearner[[paste0("D", k, "_X")]] <-
      D_X_res_list[[k]]$oos_resid_bylearner
  }#FOR

  # Organize output
  ddml_fit <- list(coef = coef, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DX = learners_DX,
                   ols_fit = ols_fit,
                   cluster_variable = cluster_variable,
                   subsamples = indxs$subsamples,
                   cv_subsamples = indxs$cv_subsamples,
                   ensemble_type = ensemble_type,
                   coefficients = coef,
                   scores = scores,
                   J = J_list,
                   coef_names = coef_names,
                   nobs = nobs,
                   sample_folds = sample_folds,
                   cv_folds = if (shortstack) NULL
                     else cv_folds,
                   shortstack = shortstack,
                   oos_resid_bylearner =
                     oos_resid_bylearner,
                   r2 = r2)

  # Print estimation completion
  elapsed <- round(proc.time()[3] - t0, 1)
  info_msg("ddml_plm: completed in ", elapsed, "s",
           silent = silent)

  # Amend class and return
  class(ddml_fit) <- c("ddml_plm", "ddml")
  return(ddml_fit)
}#DDML_PLM

