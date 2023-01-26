#' Title
#'
#' @param y abc
#' @param D abc
#' @param Z abc
#' @param X abc
#' @param learners abc
#' @param learners_DXZ abc
#' @param learners_ZX abc
#' @param sample_folds abc
#' @param ensemble_type abc
#' @param cv_folds abc
#' @param subsamples_Z0 abc
#' @param subsamples_Z1 abc
#' @param cv_subsamples_list_Z0 abc
#' @param cv_subsamples_list_Z1 abc
#' @param silent abc
#'
#' @return abc
#' @export
#'
#' @examples
#' 1 + 1
ddml_late <- function(y, D, Z, X,
                     learners,
                     learners_DXZ = learners,
                     learners_ZX = learners,
                     sample_folds = 2,
                     ensemble_type = "average",
                     cv_folds = 5,
                     subsamples_Z0 = NULL,
                     subsamples_Z1 = NULL,
                     cv_subsamples_list_Z0 = NULL,
                     cv_subsamples_list_Z1 = NULL,
                     silent = F) {
  # Data parameters
  nobs <- length(y)
  is_Z0 <- which(Z == 0)
  nobs_Z0 <- length(is_Z0)
  nobs_Z1 <- nobs - nobs_Z0
  nensb <- length(ensemble_type)

  # Create sample fold tuple by treatment
  if (is.null(subsamples_Z0) | is.null(subsamples_Z1)) {
    subsamples_Z0 <- generate_subsamples(nobs_Z0, sample_folds)
    subsamples_Z1 <- generate_subsamples(nobs_Z1, sample_folds)
  }#IF
  sample_folds <- length(subsamples_Z0)

  # Create cv-subsamples tuple by treatment
  if (is.null(cv_subsamples_list_Z0) | is.null(cv_subsamples_list_Z1)) {
    cv_subsamples_list_Z0 <- rep(list(NULL), sample_folds)
    cv_subsamples_list_Z1 <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      nobs_Z0_k <- nobs_Z0 - length(subsamples_Z0[[k]])
      nobs_Z1_k <- nobs_Z1 - length(subsamples_Z1[[k]])
      cv_subsamples_list_Z0[[k]] <- generate_subsamples(nobs_Z0_k, cv_folds)
      cv_subsamples_list_Z1[[k]] <- generate_subsamples(nobs_Z1_k, cv_folds)
    }# FOR
  }#IF

  # Merge subsamples across treatment and create auxilliary control matrix
  subsamples <- subsamples_Z0
  cv_subsamples_list <- cv_subsamples_list_Z0
  auxilliary_X_Z0 <- rep(list(NULL), sample_folds)
  auxilliary_X_Z1 <- rep(list(NULL), sample_folds)
  for (k in 1:sample_folds) {
    # Sample folds
    subsamples[[k]] <- sort(c(c(1:nobs)[is_Z0][subsamples_Z0[[k]]],
                              c(1:nobs)[-is_Z0][subsamples_Z1[[k]]]))
    # CV folds
    nobs_k <- nobs - length(subsamples[[k]])
    is_Z0_k <- which(Z[-subsamples[[k]]] == 0)
    is_Z1_k <- which(Z[-subsamples[[k]]] == 1)
    for (j in 1:cv_folds) {
      indx_Z0 <- is_Z0_k[cv_subsamples_list_Z0[[k]][[j]]]
      indx_Z1 <- is_Z1_k[cv_subsamples_list_Z1[[k]][[j]]]
      cv_subsamples_list[[k]][[j]] <- sort(c(indx_Z0, indx_Z1))
    }#FOR

    # Auxilliary X
    auxilliary_X_Z1[[k]] <- X[-is_Z0, , drop=F][subsamples_Z1[[k]], , drop=F]
    auxilliary_X_Z0[[k]] <- X[is_Z0, , drop=F][subsamples_Z0[[k]], , drop=F]
  }#FOR

  # Print to progress to console
  if (!silent) cat("DDML estimation in progress. \n")

  # Compute estimates of E[y|Z=0,X]
  y_X_Z0_res <- crosspred(y[is_Z0], X[is_Z0, , drop = F],
                          learners = learners, ensemble_type = ensemble_type,
                          cv_subsamples_list = cv_subsamples_list_Z0,
                          subsamples = subsamples_Z0,
                          silent = silent, progress = "E[Y|Z=0,X]: ",
                          auxilliary_X = auxilliary_X_Z1)
  update_progress(silent)

  # Compute estimates of E[y|Z=1,X]
  y_X_Z1_res <- crosspred(y[-is_Z0], X[-is_Z0, , drop = F],
                          learners = learners, ensemble_type = ensemble_type,
                          cv_subsamples_list = cv_subsamples_list_Z1,
                          subsamples = subsamples_Z1,
                          silent = silent, progress = "E[Y|Z=1,X]: ",
                          auxilliary_X = auxilliary_X_Z0)
  update_progress(silent)

  # Compute estimates of E[D|Z=0,X]
  D_X_Z0_res <- crosspred(D[is_Z0], X[is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          cv_subsamples_list = cv_subsamples_list_Z0,
                          subsamples = subsamples_Z0,
                          silent = silent, progress = "E[D|Z=0,X]: ",
                          auxilliary_X = auxilliary_X_Z1)
  update_progress(silent)

  # Compute estimates of E[D|Z=1,X]
  D_X_Z1_res <- crosspred(D[-is_Z0], X[-is_Z0, , drop = F],
                          learners = learners_DXZ,
                          ensemble_type = ensemble_type,
                          cv_subsamples_list = cv_subsamples_list_Z1,
                          subsamples = subsamples_Z1,
                          silent = silent, progress = "E[D|Z=1,X]: ",
                          auxilliary_X = auxilliary_X_Z0)
  update_progress(silent)

  # Compute estimates of E[Z|X]
  Z_X_res <- crosspred(Z, X,
                       learners = learners_ZX, ensemble_type = ensemble_type,
                       cv_subsamples_list = cv_subsamples_list,
                       subsamples = subsamples,
                       compute_insample_predictions = F,
                       silent = silent, progress = "E[Z|X]: ")
  update_progress(silent)

  # Check whether multiple ensembles are computed simultaneously
  multiple_ensembles <- nensb > 1

  # Construct reduced form variables
  l_Z0 <- l_Z1 <- p_Z0 <- p_Z1 <- matrix(0, nobs, nensb)
  l_Z0[is_Z0, ] <- y_X_Z0_res$oos_fitted
  l_Z1[-is_Z0, ] <- y_X_Z1_res$oos_fitted
  p_Z0[is_Z0, ] <- D_X_Z0_res$oos_fitted
  p_Z1[-is_Z0, ] <- D_X_Z1_res$oos_fitted
  if (!multiple_ensembles) {
    for (k in 1:sample_folds) {
      l_Z1[is_Z0][subsamples_Z0[[k]]] <- y_X_Z1_res$auxilliary_fitted[[k]]
      l_Z0[-is_Z0][subsamples_Z1[[k]]] <- y_X_Z0_res$auxilliary_fitted[[k]]
      p_Z1[is_Z0][subsamples_Z0[[k]]] <- D_X_Z1_res$auxilliary_fitted[[k]]
      p_Z0[-is_Z0][subsamples_Z1[[k]]] <- D_X_Z0_res$auxilliary_fitted[[k]]
    }#FOR
  } else {
    for (k in 1:sample_folds) {
      l_Z1[is_Z0, ][subsamples_Z0[[k]], ] <- y_X_Z1_res$auxilliary_fitted[[k]]
      l_Z0[-is_Z0, ][subsamples_Z1[[k]], ] <- y_X_Z0_res$auxilliary_fitted[[k]]
      p_Z1[is_Z0, ][subsamples_Z0[[k]], ] <- D_X_Z1_res$auxilliary_fitted[[k]]
      p_Z0[-is_Z0, ][subsamples_Z1[[k]], ] <- D_X_Z0_res$auxilliary_fitted[[k]]
    }#FOR
  }#IF
  r_X <- Z_X_res$oos_fitted

  # Compute the ATE using the constructed variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  D_copy <- matrix(rep(D, nensb), nobs, nensb)
  Z_copy <- matrix(rep(Z, nensb), nobs, nensb)
  numerator <- colMeans(Z_copy * (y_copy - l_Z1) / r_X +
                          (1 - Z_copy) * (y_copy - l_Z0) / (1 - r_X) +
                          l_Z1 - l_Z0)
  denominator <- colMeans(Z_copy * (D_copy - p_Z1) / r_X +
                            (1 - Z_copy) * (D_copy - p_Z0) / (1 - r_X) +
                            p_Z1 - p_Z0)
  late <- numerator / denominator
  names(late) <- ensemble_type

  # Organize complementary ensemble output
  weights <- list(y_X_Z0 = y_X_Z0_res$weights,
                  y_X_Z1 = y_X_Z1_res$weights,
                  D_X_Z0 = D_X_Z0_res$weights,
                  D_X_Z1 = D_X_Z1_res$weights,
                  Z_X = Z_X_res$weights)
  if (multiple_ensembles) {
    for (j in 1:3) {
      dimnames(weights[[j]]) <- list(NULL, ensemble_type, NULL)
    }#FOR
  }#IF

  # Store complementary ensemble output
  mspe <- list(y_X_Z0 = y_X_Z0_res$mspe,
               y_X_Z1 = y_X_Z1_res$mspe,
               D_X_Z0 = D_X_Z0_res$mspe,
               D_X_Z1 = D_X_Z1_res$mspe,
               Z_X = Z_X_res$mspe)

  # Organize output
  ddml_fit <- list(late = late, weights = weights, mspe = mspe,
                   learners = learners,
                   learners_DXZ = learners_DXZ,
                   learners_ZX = learners_ZX,
                   subsamples_Z0 = subsamples_Z0,
                   subsamples_Z1 = subsamples_Z1,
                   cv_subsamples_list_Z0 = cv_subsamples_list_Z0,
                   cv_subsamples_list_Z1 = cv_subsamples_list_Z1,
                   ensemble_type = ensemble_type)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- c("ddml_late")
  return(ddml_fit)
}#DDML_LATE
