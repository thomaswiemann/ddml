# Collection of small internal functions

# Resolve estimator progress messages.
#
# Merges user-supplied message overrides (passed via dots by
# internal callers like ddml_ate) with estimator defaults.
# Standard start/finish templates are derived from `name`.
#
# @param dots The `list(...)` captured in the estimator.
# @param name Estimator name, e.g. "ddml_apo".
# @param labels Named list of equation-specific labels,
#   e.g. list(y_X = "E[Y|D=1,X]", D_X = "P(D=1|X)").
# @return A named list of message strings.
resolve_messages <- function(dots, name, labels = list()) {
  defaults <- c(
    list(start = paste0(name, ": estimating (%s)"),
         finish = paste0(name, ": completed in %s s")),
    labels)
  user <- dots[["messages"]]
  if (is.null(user)) return(defaults)
  c(user, defaults[setdiff(names(defaults), names(user))])
}#RESOLVE_MESSAGES

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
trim_propensity_scores <- function(m_X, trim, ensemble_type,
                                   silent = FALSE) {
  nensb <- length(ensemble_type)
  for (j in seq_len(nensb)) {
    indx_trim_0 <- which(m_X[, j] <= trim)
    indx_trim_1 <- which(m_X[, j] >= 1 - trim)
    ntrim <- length(c(indx_trim_0, indx_trim_1))
    if (ntrim > 0) {
      if (!silent) {
        if (nensb == 1) {
          warning(paste0(ntrim,
                         " propensity scores were trimmed."))
        } else {
          warning(paste0(ensemble_type[j], ": ", ntrim,
                         " propensity scores were trimmed."))
        }#IFELSE
      }#IF
      m_X[indx_trim_0, j] <- trim
      m_X[indx_trim_1, j] <- 1 - trim
    }#IF
  }#FOR
  # Return trimmed scores
  m_X
}#TRIM_PROPENSITY_SCORES

# Detect whether learners is a single-learner spec or stacking.
# Single: list(what = fn, args = ...) — first element is not a list.
# Stacking: list(list(what = fn), list(what = fn2)) — first element is a list.
is_single_learner <- function(learners) {
  is.list(learners) && !is.list(learners[[1]])
}#IS_SINGLE_LEARNER

# Ensure all learner specs have $what set.
# Accepts $fun as deprecated alias. Removes $fun after resolving.
normalize_learners <- function(learners) {
  resolve <- function(l) {
    if (!is.null(l$what)) return(l$what)
    if (!is.null(l$fun)) {
      if (is.null(getOption("ddml.fun_deprecated_warned"))) {
        message("Note: 'fun' in learner specifications ",
                "is deprecated. Use 'what' instead.")
        options(ddml.fun_deprecated_warned = TRUE)
      }#IF
      return(l$fun)
    }#IF
    stop("Learner must have a 'what' or 'fun' element.",
         call. = FALSE)
  }#RESOLVE

  if (is_single_learner(learners)) {
    learners$what <- resolve(learners)
    learners$fun <- NULL
    return(learners)
  }#IF
  for (i in seq_along(learners)) {
    learners[[i]]$what <- resolve(learners[[i]])
    learners[[i]]$fun <- NULL
  }#FOR
  learners
}#NORMALIZE_LEARNERS

# Build a CEF-like result from pre-computed per-learner predictions.
# When crossval_resid + subsamples are available, recomputes
# per-fold weights from inner-CV residuals (exact). Otherwise
# uses sample-fold residuals (approximate, exact for shortstacking).
build_CEF_from_crossfit <- function(y, crossfit_fitted_eq,
                                    ensemble_type,
                                    custom_ensemble_weights,
                                    crossval_resid = NULL,
                                    subsamples = NULL,
                                    auxiliary_fitted_bylearner = NULL) {
  nobs <- length(y)
  cf <- as.matrix(crossfit_fitted_eq)
  nlearners <- ncol(cf)
  dummy_learners <- lapply(seq_len(nlearners),
    function(i) list(what = identity))

  # Per-learner OOS residuals, MSPE, and R-squared
  oos_resid <- drop(y) - cf
  mspe <- colMeans(oos_resid^2)

  if (!is.null(crossval_resid) && !is.null(crossval_resid[[1]]) &&
      !is.null(subsamples)) {
    # Per-fold weight recomputation from inner-CV residuals
    K <- length(subsamples)
    nensb <- NULL
    oos_fitted <- matrix(0, nobs, 1)
    all_weights <- vector("list", K)
    for (k in seq_len(K)) {
      train_idx <- setdiff(seq_len(nobs), subsamples[[k]])
      fakecv_k <- list(
        oos_resid = crossval_resid[[k]],
        mspe = colMeans(crossval_resid[[k]]^2))
      ew_k <- ensemble_weights(
        y[train_idx], cf[train_idx, ],
        type = ensemble_type,
        learners = dummy_learners,
        cv_results = fakecv_k,
        custom_weights = custom_ensemble_weights,
        silent = TRUE)
      all_weights[[k]] <- ew_k$weights
      if (is.null(nensb)) nensb <- ncol(ew_k$weights)
      if (ncol(oos_fitted) < nensb) {
        oos_fitted <- matrix(0, nobs, nensb)
      }#IF
      oos_fitted[subsamples[[k]], ] <- cf[subsamples[[k]], ] %*%
        ew_k$weights
    }#FOR
    weights <- array(0, dim = c(nlearners, nensb, K))
    for (k in seq_len(K)) weights[, , k] <- all_weights[[k]]
    dimnames(weights) <- list(NULL, colnames(all_weights[[1]]),
                              paste("sample fold ", seq_len(K)))
  } else {
    # Global weights from sample-fold residuals
    fakecv <- list(oos_resid = oos_resid,
                   mspe = mspe)
    ew <- ensemble_weights(
      y, cf, type = ensemble_type,
      learners = dummy_learners,
      cv_results = fakecv,
      custom_weights = custom_ensemble_weights,
      silent = TRUE)
    weights <- ew$weights
    oos_fitted <- cf %*% weights
  }#IFELSE

  # Propagate ensemble type names to oos_fitted columns
  ens_names <- if (length(dim(weights)) == 3) {
    colnames(weights[, , 1])
  } else {
    colnames(weights)
  }
  if (!is.null(ens_names)) colnames(oos_fitted) <- ens_names

  # Fold-level auxiliary predictions (ATE/ATT/LATE extrapolation)
  auxiliary_fitted <- NULL
  if (!is.null(auxiliary_fitted_bylearner)) {
    K <- length(auxiliary_fitted_bylearner)
    auxiliary_fitted <- vector("list", K)
    if (length(dim(weights)) == 3) {
      for (k in seq_len(K)) {
        auxiliary_fitted[[k]] <-
          as.matrix(auxiliary_fitted_bylearner[[k]]) %*%
          weights[, , k]
      }#FOR
    } else {
      for (k in seq_len(K)) {
        auxiliary_fitted[[k]] <-
          as.matrix(auxiliary_fitted_bylearner[[k]]) %*%
          weights
      }#FOR
    }#IFELSE
  }#IF
  y_var <- as.numeric(stats::var(y))
  r2 <- if (y_var > 0) 1 - mspe / y_var else
    rep(NA_real_, length(mspe))

  list(oos_fitted = oos_fitted,
       weights = weights,
       mspe = mspe,
       r2 = r2,
       auxiliary_fitted = auxiliary_fitted,
       crossfit_fitted = cf,
       crossfit_resid = oos_resid,
       crossval_resid = crossval_resid)
}#BUILD_CEF_FROM_CROSSFIT

validate_fitted_splits_pair <- function(fitted, splits,
                                        w_cv = FALSE) {
  if (is.null(fitted)) return(invisible(NULL))
  if (is.null(splits)) {
    stop("Argument 'splits' must be supplied when 'fitted' is supplied.")
  }#IF
  if (is.null(splits$subsamples)) {
    stop("splits must contain 'subsamples' when 'fitted' is supplied.")
  }#IF
  if (w_cv && is.null(splits$cv_subsamples)) {
    stop(paste("splits must contain 'cv_subsamples' for data-driven",
               "stacking when 'fitted' is supplied."))
  }#IF
}#VALIDATE_FITTED_SPLITS_PAIR

build_fitted_entry <- function(res, save_crossval,
                               include_auxiliary = FALSE) {
  entry <- list(ensemble_fitted = res$oos_fitted,
                crossfit_fitted = res$crossfit_fitted,
                crossfit_resid = res$crossfit_resid)
  if (save_crossval) {
    entry$crossval_resid <- res$crossval_resid
  }#IF
  if (include_auxiliary) {
    entry$auxiliary_fitted <- res$auxiliary_fitted
    entry$auxiliary_fitted_bylearner <-
      res$auxiliary_fitted_bylearner
  }#IF
  entry
}#BUILD_FITTED_ENTRY

build_fitted_from_list <- function(res_list, save_crossval) {
  lapply(res_list, build_fitted_entry,
         save_crossval = save_crossval)
}#BUILD_FITTED_FROM_LIST

get_crossfit_resid_for_eq <- function(fitted, eq) {
  entry <- fitted[[eq]]
  if (!is.null(entry) && !is.null(entry$crossfit_resid)) {
    return(entry$crossfit_resid)
  }#IF
  m <- regmatches(eq, regexec("^([A-Za-z]+)(\\d+)(_\\w+)$",
                              eq))[[1]]
  if (length(m) == 4) {
    group <- paste0(m[2], m[4])
    idx <- as.integer(m[3])
    entry <- fitted[[group]]
    if (is.list(entry) && length(entry) >= idx) {
      return(entry[[idx]]$crossfit_resid)
    }#IF
  }#IF
  NULL
}#GET_CROSSFIT_RESID_FOR_EQ

normalize_splits <- function(splits = NULL,
                             by_label = NULL, ...) {
  dots <- list(...)
  subsamples <- dots[["subsamples"]]
  cv_subsamples <- dots[["cv_subsamples"]]
  # Handle cv_subsamples_list (older deprecated name)
  if (!is.null(dots[["cv_subsamples_list"]])) {
    if (!is.null(cv_subsamples))
      stop("Specify cv_subsamples or cv_subsamples_list, ",
           "not both.")
    message("Note: cv_subsamples_list has been renamed to ",
            "cv_subsamples.")
    cv_subsamples <- dots[["cv_subsamples_list"]]
  }#IF
  # Handle grouped split args (ATE/ATT/LATE)
  subsamples_by <- NULL
  cv_subsamples_by <- NULL
  if (!is.null(by_label)) {
    sub_by <- paste0("subsamples_by", by_label)
    cv_sub_by <- paste0("cv_subsamples_by", by_label)
    subsamples_by <- dots[[sub_by]]
    cv_subsamples_by <- dots[[cv_sub_by]]
  }#IF
  legacy_used <- !is.null(subsamples) || !is.null(cv_subsamples) ||
    !is.null(subsamples_by) || !is.null(cv_subsamples_by)
  if (legacy_used) {
    warning("Deprecated split arguments detected. ",
            "Use 'splits' instead.", call. = FALSE)
  }#IF
  if (is.null(splits) && !legacy_used) return(NULL)
  if (is.null(splits)) splits <- list()
  if (!is.list(splits)) stop("'splits' must be a list.")
  if (is.null(splits$subsamples) && !is.null(subsamples))
    splits$subsamples <- subsamples
  if (is.null(splits$cv_subsamples) && !is.null(cv_subsamples))
    splits$cv_subsamples <- cv_subsamples
  if (!is.null(by_label)) {
    if (is.null(splits[[sub_by]]) && !is.null(subsamples_by))
      splits[[sub_by]] <- subsamples_by
    if (is.null(splits[[cv_sub_by]]) &&
        !is.null(cv_subsamples_by))
      splits[[cv_sub_by]] <- cv_subsamples_by
  }#IF
  splits
}#NORMALIZE_SPLITS

# Input validation checks for DDML estimators
validate_inputs <- function(y = NULL, D = NULL, X = NULL, Z = NULL,
                            learners = NULL,
                            sample_folds = NULL, cv_folds = NULL,
                            ensemble_type = NULL, trim = NULL,
                            weights = NULL,
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
      is_single <- is_single_learner(learners)
      if (!is_single) {
        for (l in learners) {
          if (!is.list(l) ||
              (is.null(l$what) && is.null(l$fun))) {
            stop("Each stacking learner must have a ",
                 "'what' or 'fun' element.",
                 call. = FALSE)
          }#IF
        }#FOR
      }#IF
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

  if (!is.null(weights) && !is.null(y)) {
    if (!is.numeric(weights) || anyNA(weights)) {
      stop("weights must be a numeric vector with no NAs.")
    }
    if (length(weights) != nobs) {
      stop("weights must have the same length as y.")
    }
  }
}#VALIDATE_INPUTS

validate_custom_weights <- function(custom_weights, learners) {
  if (is.null(custom_weights)) return(invisible(NULL))
  if (!is.numeric(custom_weights)) {
    stop("custom_ensemble_weights must be numeric.")
  }
  custom_weights <- as.matrix(custom_weights)
  n_learners <- if (is_single_learner(learners)) 1 else length(learners)
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

# Update ensemble info from CEF result weights or oos_fitted.
update_ensemble_info <- function(res_weights = NULL,
                                 oos_fitted = NULL) {
  if (!is.null(res_weights)) {
    ensemble_type <- dimnames(res_weights)[[2]]
  } else if (!is.null(oos_fitted)) {
    ensemble_type <- colnames(as.matrix(oos_fitted))
  } else {
    ensemble_type <- NULL
  }#IFELSE
  nensb <- if (is.null(ensemble_type)) 1L
    else length(ensemble_type)
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
                             parallel = NULL,
                             fitted = NULL) {
  nM <- ncol(M)
  res_list <- vector("list", nM)
  for (k in seq_len(nM)) {
    fit_k <- if (!is.null(fitted)) {
      fitted[[k]]
    }#IF
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
      parallel = parallel,
      fitted = fit_k)
  }#FOR
  res_list
}#COMPUTE_CEF_LIST
