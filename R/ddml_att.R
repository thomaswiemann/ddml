#' @rdname ddml_ate
#'
#' @export
ddml_att <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 2,
                     ensemble_type = "nnls",
                     shortstack = FALSE,
                     cv_folds = 5,
                     custom_ensemble_weights = NULL,
                     custom_ensemble_weights_DX = custom_ensemble_weights,
                     subsamples_D0 = NULL,
                     subsamples_D1 = NULL,
                     cv_subsamples_list_D0 = NULL,
                     cv_subsamples_list_D1 = NULL,
                     trim = 0.01,
                     silent = FALSE) {
  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)
  nobs_D0 <- length(is_D0)
  nobs_D1 <- nobs - nobs_D0

  # Create sample fold tuple by treatment
  if (is.null(subsamples_D0) | is.null(subsamples_D1)) {
    subsamples_D0 <- generate_subsamples(nobs_D0, sample_folds)
    subsamples_D1 <- generate_subsamples(nobs_D1, sample_folds)
  }#IF
  sample_folds <- length(subsamples_D0)

  # Create cv-subsamples tuple by treatment
  if (is.null(cv_subsamples_list_D0) | is.null(cv_subsamples_list_D1)) {
    cv_subsamples_list_D0 <- rep(list(NULL), sample_folds)
    cv_subsamples_list_D1 <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_D0_k <- nobs_D0 - length(subsamples_D0[[k]])
      nobs_D1_k <- nobs_D1 - length(subsamples_D1[[k]])
      cv_subsamples_list_D0[[k]] <- generate_subsamples(nobs_D0_k, cv_folds)
      cv_subsamples_list_D1[[k]] <- generate_subsamples(nobs_D1_k, cv_folds)
    }# FOR
  }#IF

  # Merge subsamples across treatment and create auxilliary control matrix
  subsamples <- subsamples_D0
  cv_subsamples_list <- cv_subsamples_list_D0
  auxilliary_X_D0 <- rep(list(NULL), sample_folds)
  auxilliary_X_D1 <- rep(list(NULL), sample_folds)
  for (k in 1:sample_folds) {
    # Sample folds
    subsamples[[k]] <- sort(c((1:nobs)[is_D0][subsamples_D0[[k]]],
                              (1:nobs)[-is_D0][subsamples_D1[[k]]]))
    # CV folds
    nobs_k <- nobs - length(subsamples[[k]])
    is_D0_k <- which(D[-subsamples[[k]]] == 0)
    is_D1_k <- which(D[-subsamples[[k]]] == 1)
    for (j in 1:cv_folds) {
      indx_D0 <- is_D0_k[cv_subsamples_list_D0[[k]][[j]]]
      indx_D1 <- is_D1_k[cv_subsamples_list_D1[[k]][[j]]]
      cv_subsamples_list[[k]][[j]] <- sort(c(indx_D0, indx_D1))
    }#FOR

    # Auxilliary X (only need treated observations)
    auxilliary_X_D1[[k]] <- X[-is_D0, , drop=F][subsamples_D1[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        custom_ensemble_weights = custom_ensemble_weights,
                        cv_subsamples_list = cv_subsamples_list_D0,
                        subsamples = subsamples_D0,
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxilliary_X = auxilliary_X_D1)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     custom_ensemble_weights = custom_ensemble_weights_DX,
                     cv_subsamples_list = cv_subsamples_list,
                     subsamples = subsamples,
                     silent = silent, progress = "E[D|X]: ")

  # Compute estimates of E[D] -- simple computation of averages here
  D_res <- get_CEF(D, matrix(1, nobs, 1),
                   learners = list(what = ols),
                   ensemble_type = "average",
                   shortstack = FALSE,
                   cv_subsamples_list = NULL,
                   subsamples = subsamples,
                   silent = TRUE)

  # Update ensemble type to account for (optional) custom weights
  ensemble_type <- dimnames(y_X_D0_res$weights)[[2]]
  nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_D0 <- matrix(0, nobs, nensb)
  g_D0[is_D0, ] <- y_X_D0_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      g_D0[-is_D0][subsamples_D1[[k]]] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      g_D0[-is_D0, ][subsamples_D1[[k]], ] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  m_X <- D_X_res$oos_fitted

  # Trim propensity scores, return warnings
  m_X_tr <- trim_propensity_scores(m_X, trim, ensemble_type)

  # Compute the ATT using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  p_copy <- matrix(rep(D_res$oos_fitted, nensb), nobs, nensb)
  psi_b <- D_copy * (y_copy - g_D0) / p_copy -
    m_X_tr * (1 - D_copy) * (y_copy - g_D0) / (p_copy * (1 - m_X_tr))
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
  oos_pred <- list(EY_D0_X = g_D0, ED_X = m_X, ED = D_res$oos_fitted)

  # Organize output
  ddml_fit <- list(att = att, weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   subsamples_D0 = subsamples_D0,
                   subsamples_D1 = subsamples_D1,
                   cv_subsamples_list_D0 = cv_subsamples_list_D0,
                   cv_subsamples_list_D1 = cv_subsamples_list_D1,
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
                                                     object$ensemble_type)
  class(coefficients) <- c("summary.ddml_att", class(coefficients))
  coefficients
}#SUMMARY.DDML_ATT

#' @rdname print.summary.ddml_ate
#'
#' @export
print.summary.ddml_att <- function(x, ...) {
  cat("ATT estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x)
}#PRINT.SUMMARY.DDML_ATT
