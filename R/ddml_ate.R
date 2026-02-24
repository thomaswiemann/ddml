#' Estimators of Average Treatment Effects.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml()]
#'
#' @description Estimators of the average treatment effect and the average
#'     treatment effect on the treated.
#'
#' @details \code{ddml_ate} and \code{ddml_att} provide double/debiased machine
#'     learning  estimators for the average treatment effect and the average
#'     treatment effect on the treated, respectively, in the interactive model
#'     given by
#'
#' \eqn{Y = g_0(D, X) + U,}
#'
#' where \eqn{(Y, D, X, U)} is a random vector such that
#'     \eqn{\operatorname{supp} D = \{0,1\}}, \eqn{E[U\vert D, X] = 0}, and
#'     \eqn{\Pr(D=1\vert X) \in (0, 1)} with probability 1,
#'     and \eqn{g_0} is an unknown nuisance function.
#'
#' In this model, the average treatment effect is defined as
#'
#' \eqn{\theta_0^{\textrm{ATE}} \equiv E[g_0(1, X) - g_0(0, X)]}.
#'
#' and the average treatment effect on the treated is defined as
#'
#' \eqn{\theta_0^{\textrm{ATT}} \equiv E[g_0(1, X) - g_0(0, X)\vert D = 1]}.
#'
#' @inheritParams ddml_plm
#' @param D The binary endogenous variable of interest.
#' @param subsamples_byD List of two lists corresponding to the two treatment
#'     levels. Each list contains vectors with sample indices for
#'     cross-fitting.
#' @param cv_subsamples_byD List of two lists, each corresponding to one of the
#'     two treatment levels. Each of the two lists contains lists, each
#'     corresponding to a subsample and contains vectors with subsample indices
#'     for cross-validation.
#' @param stratify Boolean for stratified cross-fitting: if \code{TRUE},
#'     subsamples are constructed to be balanced across treatment levels.
#' @param trim Number in (0, 1) for trimming the estimated propensity scores at
#'     \code{trim} and \code{1-trim}.
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
#' @return \code{ddml_ate} and \code{ddml_att} return an object of S3 class
#'     \code{ddml_ate} and \code{ddml_att}, respectively. An object of class
#'     \code{ddml_ate} or \code{ddml_att} is a list containing
#'     the following components:
#'     \describe{
#'         \item{\code{ate} / \code{att}}{A vector with the average treatment
#'             effect / average treatment effect on the treated estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the computation
#'             of scores. Used in [ddml::summary.ddml()].}
#'         \item{\code{oos_pred}}{List of matrices, providing the reduced form
#'             predicted values.}
#'         \item{\code{learners},\code{learners_DX},\code{cluster_variable},
#'             \code{subsamples_byD},\code{cv_subsamples_byD},
#'             \code{ensemble_type}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
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
#' # Estimate the average treatment effect using a single base learner, ridge.
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(what = mdl_glmnet,
#'                                     args = list(alpha = 0)),
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
#'
#' # Estimate the average treatment effect using short-stacking with base
#' #     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the ATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' ate_fit <- ddml_ate(y, D, X,
#'                     learners = list(list(fun = ols),
#'                                     list(fun = mdl_glmnet),
#'                                     list(fun = mdl_glmnet,
#'                                          args = list(alpha = 0))),
#'                     ensemble_type = 'nnls',
#'                     custom_ensemble_weights = weights_everylearner,
#'                     shortstack = TRUE,
#'                     sample_folds = 2,
#'                     silent = TRUE)
#' summary(ate_fit)
ddml_ate <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 10,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 10,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     cluster_variable = seq_along(y),
                     stratify = TRUE,
                     subsamples = NULL,
                     subsamples_byD = NULL,
                     cv_subsamples = NULL,
                     cv_subsamples_byD = NULL,
                     trim = 0.01,
                     silent = FALSE,
                     parallel = NULL) {
  # Validate inputs
  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds, cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)

  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- !shortstack &
    any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (class(learners[[1]]) != "function")

  # Create crossfitting and cv tuples
  indxs <- get_sample_splits(cluster_variable = cluster_variable,
                             sample_folds = sample_folds,
                             cv_folds = if (w_cv) cv_folds,
                             D = D, stratify = stratify,
                             subsamples = subsamples,
                             subsamples_byD = subsamples_byD,
                             cv_subsamples = cv_subsamples,
                             cv_subsamples_byD = cv_subsamples_byD)
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D)

  # Estimation start
  t0 <- proc.time()[3]
  mode_str <- if (!is.null(parallel)) {
    p <- parse_parallel(parallel)
    paste0("parallel, ", p$num_cores, " cores")
  } else {
    "sequential"
  }
  info_msg("ddml_ate: estimating (", mode_str, ")",
           silent = silent)

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = FALSE],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = indxs$subsamples_byD[[1]],
                        cv_subsamples = indxs$cv_subsamples_byD[[1]],
                        silent = silent, label = "E[Y|D=0,X]",
                        auxiliary_X = get_auxiliary_X(indxs$aux_indx[[1]], X),
                        parallel = parallel)

  # Compute estimates of E[y|D=1,X]
  y_X_D1_res <- get_CEF(y[-is_D0], X[-is_D0, , drop = FALSE],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = indxs$subsamples_byD[[2]],
                        cv_subsamples = indxs$cv_subsamples_byD[[2]],
                        silent = silent, label = "E[Y|D=1,X]",
                        auxiliary_X = get_auxiliary_X(indxs$aux_indx[[2]], X),
                        parallel = parallel)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, label = "E[D|X]",
                     parallel = parallel)

  # Update ensemble type to account for (optional) custom weights
  ensb_info <- update_ensemble_info(y_X_D0_res$weights)
  ensemble_type <- ensb_info$ensemble_type
  nensb <- ensb_info$nensb
  multiple_ensembles <- ensb_info$multiple_ensembles

  # Construct reduced form variables
  g_X_byD <- extrapolate_CEF(D = D,
                             CEF_res_byD = list(list(y_X_D0_res, d=0),
                                                list(y_X_D1_res, d=1)),
                             aux_indx = indxs$aux_indx)
  m_X <- D_X_res$oos_fitted

  # Trim propensity scores, return warnings
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)

  # Compute the ATE using the constructed variables
  if (!multiple_ensembles) {
    g0 <- g_X_byD[, , 1]
    g1 <- g_X_byD[, , 2]
    m <- as.vector(m_X_tr)
    psi_b <- matrix(
      D * (y - g1) / m - (1 - D) * (y - g0) / (1 - m) + g1 - g0,
      nobs, 1)
    ate <- mean(psi_b)
    names(ate) <- ensemble_type
    psi_a <- matrix(-1, nobs, 1)
  } else {
    y_copy <- matrix(rep(y, nensb), nobs, nensb)
    D_copy <- matrix(rep(D, nensb), nobs, nensb)
    psi_b <- D_copy * (y_copy - g_X_byD[, , 2]) / m_X_tr -
      (1 - D_copy) * (y_copy - g_X_byD[, , 1]) / (1 - m_X_tr) +
      g_X_byD[, , 2] - g_X_byD[, , 1]
    ate <- colMeans(psi_b)
    names(ate) <- ensemble_type
    psi_a <- matrix(-1, nobs, nensb)
  }#IFELSE

  # Compute scores and Jacobian from psi_a/psi_b
  scores <- lapply(seq_len(nensb), function(j) {
    as.matrix(psi_a[, j] * ate[j] + psi_b[, j])
  })
  J_list <- lapply(seq_len(nensb), function(j) {
    as.matrix(mean(psi_a[, j]))
  })
  coef_names <- "ATE"

  # Organize complementary ensemble output
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  y_X_D1 = y_X_D1_res$weights,
                  D_X = D_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               y_X_D1 = y_X_D1_res$mspe,
               D_X = D_X_res$mspe)

  # Organize reduced form predicted values
  oos_pred <- list(EY_D0_X = g_X_byD[, , 1],
                   EY_D1_X = g_X_byD[, , 2],
                   ED_X = m_X)

  # Organize output
  ddml_fit <- list(ate = ate, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   cluster_variable = cluster_variable,
                   subsamples_byD = indxs$subsamples_byD,
                   cv_subsamples_byD =
                     indxs$cv_subsamples_byD,
                   ensemble_type = ensemble_type,
                   coefficients = ate,
                   scores = scores,
                   J = J_list,
                   coef_names = coef_names,
                   nobs = nobs,
                   sample_folds = sample_folds,
                   cv_folds = if (shortstack) NULL
                     else cv_folds,
                   shortstack = shortstack)

  # Print estimation completion
  elapsed <- round(proc.time()[3] - t0, 1)
  info_msg("ddml_ate: completed in ", elapsed, "s",
           silent = silent)

  # Amend class and return
  class(ddml_fit) <- c("ddml_ate", "ddml")
  return(ddml_fit)
}#DDML_ATE
