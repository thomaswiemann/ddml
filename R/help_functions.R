# Collection of small internal functions

# Simple generalized inverse wrapper.
csolve <- function(X) {
  # Attempt inversion
  X_inv <- tryCatch(solve(X), error = function(e) NA)
  # If inversion failed, calculate generalized inverse
  if (any(is.na(X_inv))) {
    X_inv <- MASS::ginv(X)
  }#IF
  # Return (generalized) inverse
  X_inv
}#CSOLVE

# Function to pull oosresid from get_CEF results
get_oosfitted <- function(res_list, j = NULL) {
  if (is.null(j)) {
    vapply(res_list, function (x) x$oos_fitted,
           FUN.VALUE = c(res_list[[1]]$oos_fitted))
  } else {
    vapply(res_list, function (x) x$oos_fitted[, j],
           FUN.VALUE = res_list[[1]]$oos_fitted[, 1])
  }#IFELSE
}#GET_OOSRESID

# Function to trim propensity scores and warn user
trim_propensity_scores <- function(m_X, trim, ensemble_type) {
  # Data parameter
  nensb <- length(ensemble_type)
  # Trim by ensemble type
  for (j in seq_len(nensb)) {
    indx_trim_0 <- which(m_X[, j] <= trim)
    indx_trim_1 <- which(m_X[, j] >= 1 - trim)
    ntrim <- length(c(indx_trim_0, indx_trim_1))
    if (ntrim > 0) {
      # Warn user
      if (nensb == 1) {
        warning(paste0(ntrim, " propensity scores were trimmed."))
      } else {
        warning(paste0(ensemble_type[j], ": ", ntrim,
                       " propensity scores were trimmed."))
      }#IFELSE
      # Replace scores by constant
      m_X[indx_trim_0, j] <- trim
      m_X[indx_trim_1, j] <- 1 - trim
    }#IF
  }#FOR
  # Return trimmed scores
  m_X
}#TRIM_PROPENSITY_SCORES

# Input validation checks for DDML estimators
validate_inputs <- function(y = NULL, D = NULL, X = NULL, Z = NULL, learners = NULL,
                            sample_folds = NULL, cv_folds = NULL,
                            ensemble_type = NULL, trim = NULL,
                            require_binary_D = FALSE) {
  nobs <- length(y)
  if (!is.null(y)) {
    if (!is.numeric(y) || anyNA(y) || nobs == 0) {
      stop("y must be a numeric vector with no NAs.")
    }
  }

  if (!is.null(D) && !is.null(y)) {
    D_mat <- as.matrix(D)
    if (!is.numeric(D_mat) || anyNA(D_mat)) {
      stop("D must be numeric with no NAs.")
    }
    if (nrow(D_mat) != nobs) {
      stop("Length/number of rows of D must match length of y.")
    }
    if (require_binary_D) {
      if (!all(D_mat %in% c(0, 1))) {
        stop("D must be binary (0 or 1).")
      }
    }
  }

  if (!is.null(X) && !is.null(y)) {
    if (NROW(X) != nobs) {
      stop("Number of rows of X must match length of y.")
    }
  }

  if (!is.null(Z) && !is.null(y)) {
    if (NROW(Z) != nobs) {
      stop("Number of rows of Z must match length of y.")
    }
  }

  if (!is.null(learners)) {
    if (!is.list(learners)) {
      stop("learners must be a list.")
    } else {
      is_single <- "what" %in% names(learners)
      if (!is_single) {
        for (l in learners) {
          if (!is.list(l) || !"fun" %in% names(l)) {
            stop("learners structure is invalid.")
          }
        }
      }
    }
  }

  if (!is.null(sample_folds)) {
    if (!is.numeric(sample_folds) || length(sample_folds) > 1 || sample_folds < 1 || sample_folds %% 1 != 0) {
      stop("sample_folds must be a positive integer.")
    }
  }

  if (!is.null(cv_folds)) {
    if (!is.numeric(cv_folds) || length(cv_folds) > 1 || cv_folds < 1 || cv_folds %% 1 != 0) {
      stop("cv_folds must be a positive integer.")
    }
  }

  if (!is.null(ensemble_type)) {
    allowed_types <- c("nnls", "nnls1", "singlebest", "ols", "average")
    if (!is.character(ensemble_type) || any(!ensemble_type %in% allowed_types)) {
      stop("ensemble_type must be one or more of: nnls, nnls1, singlebest, ols, average.")
    }
  }

  if (!is.null(trim)) {
    if (!is.numeric(trim) || length(trim) > 1 || trim <= 0 || trim >= 0.5) {
      stop("trim must be a numeric value strictly between 0 and 0.5.")
    }
  }
}#VALIDATE_INPUTS

validate_custom_weights <- function(custom_weights, learners) {
  if (is.null(custom_weights)) return(invisible(NULL))
  if (!is.numeric(custom_weights)) {
    stop("custom_ensemble_weights must be numeric.")
  }
  custom_weights <- as.matrix(custom_weights)
  n_learners <- if ("what" %in% names(learners)) 1 else length(learners)
  if (nrow(custom_weights) != n_learners) {
    stop("Number of rows in custom_ensemble_weights must match the number of base learners.")
  }
}#VALIDATE_CUSTOM_WEIGHTS

# Compute ncustom and nensb from ensemble config.
compute_ncustom_nensb <- function(ensemble_type,
                                  custom_ensemble_weights) {
  ncustom <- ncol(custom_ensemble_weights)
  ncustom <- if (is.null(ncustom)) 0L else ncustom
  nensb <- length(ensemble_type) + ncustom
  list(ncustom = ncustom, nensb = nensb)
}#COMPUTE_NCUSTOM_NENSB

# Update ensemble info from CEF result weights.
update_ensemble_info <- function(res_weights) {
  ensemble_type <- dimnames(res_weights)[[2]]
  nensb <- if (is.null(ensemble_type)) 1L else length(ensemble_type)
  multiple_ensembles <- nensb > 1
  list(ensemble_type = ensemble_type, nensb = nensb,
       multiple_ensembles = multiple_ensembles)
}#UPDATE_ENSEMBLE_INFO

# Compute CEF for each column of M, collecting results in a list.
compute_CEF_list <- function(M, X, Z = NULL,
                             learners, ensemble_type,
                             shortstack,
                             custom_ensemble_weights,
                             subsamples, cv_subsamples,
                             compute_insample_predictions = FALSE,
                             silent = FALSE,
                             label_prefix, label_suffix,
                             parallel = NULL) {
  nM <- ncol(M)
  res_list <- vector("list", nM)
  for (k in seq_len(nM)) {
    res_list[[k]] <- get_CEF(
      M[, k, drop = FALSE], X, Z = Z,
      learners = learners,
      ensemble_type = ensemble_type,
      shortstack = shortstack,
      custom_ensemble_weights = custom_ensemble_weights,
      subsamples = subsamples,
      cv_subsamples = cv_subsamples,
      compute_insample_predictions =
        compute_insample_predictions,
      silent = silent,
      label = paste0(label_prefix, k, label_suffix),
      parallel = parallel)
  }#FOR
  res_list
}#COMPUTE_CEF_LIST
