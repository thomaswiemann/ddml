#' @rdname ddml_ate
#'
#' @export
ddml_att <- function(y, D, X,
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
                     balance_on = "clusters",
                     subsamples = NULL,
                     subsamples_byD = NULL,
                     cv_subsamples = NULL,
                     cv_subsamples_byD = NULL,
                     trim = 0.01,
                     silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)

  # Check whether ddml uses conventional stacking w/ data driven weights
  w_cv <- !shortstack &
    any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
    (class(learners[[1]]) != "function")

  # Create crossfitting and cv tuples
  indxs <- get_all_indx(cluster_variable = cluster_variable,
                        sample_folds = sample_folds, cv_folds = cv_folds,
                        D = D, by_D = TRUE, stratify = stratify,
                        balance_on = balance_on,
                        subsamples = subsamples,
                        subsamples_byD = subsamples_byD,
                        cv_subsamples = cv_subsamples,
                        cv_subsamples_byD = cv_subsamples_byD,
                        compute_cv_indices = w_cv, compute_aux_X_indices = TRUE)

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        subsamples = indxs$subsamples_byD[[1]],
                        cv_subsamples = indxs$cv_subsamples_byD[[1]],
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxiliary_X = get_auxiliary_X(indxs$aux_indx[[1]], X))

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     subsamples = indxs$subsamples,
                     cv_subsamples = indxs$cv_subsamples,
                     silent = silent, progress = "E[D|X]: ")

  # Compute estimates of E[D] -- simple computation of averages here
  D_res <- get_CEF(D, matrix(1, nobs, 1),
                   learners = list(what = ols),
                   ensemble_type = "average",
                   shortstack = FALSE,
                   cv_subsamples = NULL,
                   subsamples = indxs$subsamples,
                   silent = TRUE)

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_D0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_X_D0<- extrapolate_CEF(D = D,
                           CEF_res_byD = list(list(fit = y_X_D0_res, d = 0)),
                           aux_indx = indxs$aux_indx)[, , 1]
  m_X <- D_X_res$oos_fitted

  # Trim propensity scores, return warnings
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)

  # Compute the ATT using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  p_copy <- matrix(rep(D_res$oos_fitted, nensb), nobs, nensb)
  psi_b <- D_copy * (y_copy - g_X_D0) / p_copy -
    m_X_tr * (1 - D_copy) * (y_copy - g_X_D0) / (p_copy * (1 - m_X_tr))
  psi_a <- -D_copy / p_copy
  att <- -colMeans(psi_b) / colMeans(psi_a)
  names(att) <- ensemble_type

  # Organize complementary ensemble output
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  D_X = D_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               D_X = D_X_res$mspe)

  # Organize reduced form predicted values
  oos_pred <- list(EY_D0_X = g_X_D0, ED_X = m_X, ED = D_res$oos_fitted)

  # Organize output
  ddml_fit <- list(att = att, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   cluster_variable = cluster_variable,
                   subsamples_byD = indxs$subsamples_byD,
                   cv_subsamples_byD = indxs$cv_subsamples_byD,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_att"
  return(ddml_fit)
}#DDML_ATT

#' @rdname summary.ddml_ate
#'
#' @export
summary.ddml_att <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- " "
  # Compute and print inference results
  coefficients <- organize_interactive_inf_results(coef = object$att,
                                                   psi_a = object$psi_a,
                                                   psi_b = object$psi_b,
                                                   ensemble_type =
                                                     object$ensemble_type,
                                                   cluster_variable =
                                                     object$cluster_variable)
  class(coefficients) <- c("summary.ddml_att", class(coefficients))
  coefficients
}#SUMMARY.DDML_ATT

#' @rdname print.summary.ddml_ate
#'
#' @export
print.summary.ddml_att <- function(x, digits = 3, ...) {
  cat("ATT estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_ATT
