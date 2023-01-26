#' Title
#'
#' @param y abc
#' @param D abc
#' @param X abc
#' @param learners abc
#' @param learners_DX abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param shortstack abc
#' @param cv_folds abc
#' @param subsamples_D0 abc
#' @param subsamples_D1 abc
#' @param cv_subsamples_list_D0 abc
#' @param cv_subsamples_list_D1 abc
#' @param silent abc
#'
#' @return abc
#' @export
#'
#' @examples
#' 1 + 1
ddml_ate <- function(y, D, X,
                     learners,
                     learners_DX = learners,
                     sample_folds = 2,
                     ensemble_type = "average",
                     shortstack = FALSE,
                     cv_folds = 5,
                     subsamples_D0 = NULL,
                     subsamples_D1 = NULL,
                     cv_subsamples_list_D0 = NULL,
                     cv_subsamples_list_D1 = NULL,
                     silent = F) {
  # Data parameters
  nobs <- length(y)
  is_D0 <- which(D == 0)
  nobs_D0 <- length(is_D0)
  nobs_D1 <- nobs - nobs_D0
  nensb <- length(ensemble_type)

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
    subsamples[[k]] <- sort(c(c(1:nobs)[is_D0][subsamples_D0[[k]]],
                              c(1:nobs)[-is_D0][subsamples_D1[[k]]]))
    # CV folds
    nobs_k <- nobs - length(subsamples[[k]])
    is_D0_k <- which(D[-subsamples[[k]]] == 0)
    is_D1_k <- which(D[-subsamples[[k]]] == 1)
    for (j in 1:cv_folds) {
      indx_D0 <- is_D0_k[cv_subsamples_list_D0[[k]][[j]]]
      indx_D1 <- is_D1_k[cv_subsamples_list_D1[[k]][[j]]]
      cv_subsamples_list[[k]][[j]] <- sort(c(indx_D0, indx_D1))
    }#FOR

    # Auxilliary X
    auxilliary_X_D1[[k]] <- X[-is_D0, , drop=F][subsamples_D1[[k]], , drop=F]
    auxilliary_X_D0[[k]] <- X[is_D0, , drop=F][subsamples_D0[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|D=0,X]
  y_X_D0_res <- get_CEF(y[is_D0], X[is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        cv_subsamples_list = cv_subsamples_list_D0,
                        subsamples = subsamples_D0,
                        silent = silent, progress = "E[Y|D=0,X]: ",
                        auxilliary_X = auxilliary_X_D1)

  # Compute estimates of E[y|D=1,X]
  y_X_D1_res <- get_CEF(y[-is_D0], X[-is_D0, , drop = F],
                        learners = learners, ensemble_type = ensemble_type,
                        shortstack = shortstack,
                        cv_subsamples_list = cv_subsamples_list_D1,
                        subsamples = subsamples_D1,
                        silent = silent, progress = "E[Y|D=1,X]: ",
                        auxilliary_X = auxilliary_X_D0)

  # Compute estimates of E[D|X]
  D_X_res <- get_CEF(D, X,
                     learners = learners_DX, ensemble_type = ensemble_type,
                     shortstack = shortstack,
                     cv_subsamples_list = cv_subsamples_list,
                     subsamples = subsamples,
                     silent = silent, progress = "E[D|,X]: ")

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  g_D0 <- g_D1 <- matrix(0, nobs, nensb)
  g_D0[is_D0, ] <- y_X_D0_res$oos_fitted
  g_D1[-is_D0, ] <- y_X_D1_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      g_D1[is_D0][subsamples_D0[[k]]] <- y_X_D1_res$auxilliary_fitted[[k]]
      g_D0[-is_D0][subsamples_D1[[k]]] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      g_D1[is_D0, ][subsamples_D0[[k]], ] <- y_X_D1_res$auxilliary_fitted[[k]]
      g_D0[-is_D0, ][subsamples_D1[[k]], ] <- y_X_D0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  m_X <- D_X_res$oos_fitted

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  ate <- colMeans(D_copy * (y_copy - g_D1) / m_X +
                    (1 - D_copy) * (y_copy - g_D0) / (1 - m_X) +
                    g_D1 - g_D0)
  names(ate) <- ensemble_type

  # Organize complementary ensemble output
  weights <- list(y_X_D0 = y_X_D0_res$weights,
                  y_X_D1 = y_X_D1_res$weights,
                  D_X = D_X_res$weights)

  # Store complementary ensemble output
  mspe <- list(y_X_D0 = y_X_D0_res$mspe,
               y_X_D1 = y_X_D1_res$mspe,
               D_X = D_X_res$mspe)

  # Organize output
  ddml_fit <- list(ate = ate, weights = weights, mspe = mspe,
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
  class(ddml_fit) <- c("ddml_ate")
  return(ddml_fit)
}#DDML_ATE
