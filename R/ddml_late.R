#' Estimator for the Local Average Treatment Effect
#'
#' @family ddml estimators
#'
#' @description Estimator for the local average treatment effect.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_late} provides a
#'     Double/Debiased Machine Learning estimator for the local average
#'     treatment effect. Under the standard instrumental variable assumptions
#'     (conditional independence, exclusion restriction, relevance, and 
#'     monotonicity) with a binary instrument \eqn{Z} and a binary treatment 
#'     \eqn{D}, the parameter is identified by the following reduced form 
#'     parameter:
#'
#' \deqn{\theta_0^{\textrm{LATE}} = \frac{E[E[Y|Z=1, X] - E[Y|Z=0, X]]}{E[E[D|Z=1, X] - E[D|Z=0, X]]}}
#'
#' where \eqn{W \equiv (Y, D, X, Z)} is the observed random vector.
#'
#' \strong{Nuisance Parameters:} The nuisance parameters are
#'     \eqn{\eta = (\ell_0, \ell_1, r_0, r_1, p)} taking true values
#'     \eqn{\ell_{z,0}(X) = E[Y|Z=z, X]}, \eqn{r_{z,0}(X) = E[D|Z=z, X]},
#'     and \eqn{p_0(X) = \Pr(Z=1|X)}.
#'
#' \strong{Neyman Orthogonal Score / Moment Equation:} The Neyman orthogonal score is:
#'
#' \deqn{m(W; \theta, \eta) = \frac{Z(Y - \ell_1(X))}{p(X)} - \frac{(1-Z)(Y-\ell_0(X))}{1-p(X)} + \ell_1(X) - \ell_0(X) - \theta\left(\frac{Z(D - r_1(X))}{p(X)} - \frac{(1-Z)(D-r_0(X))}{1-p(X)} + r_1(X) - r_0(X)\right)}
#'
#' \strong{Jacobian:}
#'
#' \deqn{J = -E[r_1(X) - r_0(X)]}
#'
#' See \code{\link{ddml-intro}} for how the influence function
#' and inference are derived from these components.
#'
#' @inheritParams ddml-intro
#' @inheritParams ddml_apo
#' @param Z Binary instrumental variable.
#' @param learners_DXZ,learners_ZX Optional arguments to allow for different
#'     base learners for estimation of \eqn{E[D \vert X, Z]}, \eqn{E[Z \vert X]}. Setup is
#'     identical to \code{learners}.
#' @param custom_ensemble_weights_DXZ,custom_ensemble_weights_ZX Optional
#'     arguments to allow for different
#'     custom ensemble weights for \code{learners_DXZ},\code{learners_ZX}. Setup
#'     is identical to \code{custom_ensemble_weights}. Note:
#'     \code{custom_ensemble_weights} and
#'     \code{custom_ensemble_weights_DXZ},\code{custom_ensemble_weights_ZX} must
#'     have the same number of columns.
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_late}, recommended keys are \code{subsamples},
#'     \code{subsamples_byZ}, \code{cv_subsamples}, and
#'     \code{cv_subsamples_byZ}.
#' @param ... Additional arguments passed to internal methods.
#'
#' @return \code{ddml_late} returns an object of S3 class
#'     \code{ddml_late} and \code{ddml}. See
#'     \code{\link{ddml-intro}} for the common output structure.
#'     Additional pass-through fields: \code{learners},
#'     \code{learners_DXZ}, \code{learners_ZX}.
#' @export
#'
#' @references
#' Imbens G, Angrist J (1994). "Identification and Estimation of Local Average
#'     Treatment Effects." Econometrica, 62(2), 467-475.
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' Z = AE98[, "samesex"]
#' X = AE98[, c("age","agefst","black","hisp","othrace","educ")]
#'
#' # Estimate the local average treatment effect using a single base learner,
#' #     ridge.
#' late_fit <- ddml_late(y, D, Z, X,
#'                       learners = list(what = mdl_glmnet,
#'                                       args = list(alpha = 0)),
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(late_fit)
#'
#' \donttest{
#' # Estimate the local average treatment effect using short-stacking with base
#' #     learners ols, lasso, and ridge. We can also use custom_ensemble_weights
#' #     to estimate the LATE using every individual base learner.
#' weights_everylearner <- diag(1, 3)
#' colnames(weights_everylearner) <- c("mdl:ols", "mdl:lasso", "mdl:ridge")
#' late_fit <- ddml_late(y, D, Z, X,
#'                       learners = list(list(what = ols),
#'                                       list(what = mdl_glmnet),
#'                                       list(what = mdl_glmnet,
#'                                            args = list(alpha = 0))),
#'                       ensemble_type = 'nnls',
#'                       custom_ensemble_weights = weights_everylearner,
#'                       shortstack = TRUE,
#'                       sample_folds = 2,
#'                       silent = TRUE)
#' summary(late_fit)
#' }
ddml_late <- function(y, D, Z, X,
                      learners,
                      learners_DXZ = learners,
                      learners_ZX = learners,
                      sample_folds = 10,
                      ensemble_type = "nnls",
                      shortstack = FALSE,
                      cv_folds = 10,
                      custom_ensemble_weights = NULL,
                      custom_ensemble_weights_DXZ = custom_ensemble_weights,
                      custom_ensemble_weights_ZX = custom_ensemble_weights,
                      cluster_variable = seq_along(y),
                      stratify = TRUE,
                      trim = 0.01,
                      silent = FALSE,
                      parallel = NULL,
                      fitted = NULL,
                      splits = NULL,
                      save_crossval = TRUE,
                      ...) {
  cl <- match.call()

  # Preliminaries --------------------------------------------------------------

  dots <- list(...)
  messages <- resolve_messages(dots, "ddml_late", list(
    y_Z0 = "E[Y|Z=0,X]",
    y_Z1 = "E[Y|Z=1,X]",
    D_Z0 = "E[D|Z=0,X]",
    D_Z1 = "E[D|Z=1,X]",
    Z_X = "E[Z|X]"))

  validate_inputs(y = y, D = D, X = X, Z = Z,
                  learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  cluster_variable = cluster_variable,
                  require_binary_D = FALSE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DXZ,
                          learners_DXZ)
  validate_custom_weights(custom_ensemble_weights_ZX,
                          learners_ZX)

  nobs <- length(y)

  validate_fitted_splits_pair(fitted, splits)

  t0 <- proc.time()[3]
  announce_start(messages, parallel, silent)

  # Reduced-form estimation ----------------------------------------------------

  # Map LATE splits/fitted to ATE's expected format.
  ate_splits <- if (!is.null(splits)) list(
    y_X_D0 = splits$y_X_Z0,
    y_X_D1 = splits$y_X_Z1,
    D_X = splits$Z_X)
  ate_fitted_rf <- if (!is.null(fitted)) list(
    y_X_D0 = fitted$y_X_Z0,
    y_X_D1 = fitted$y_X_Z1,
    D_X = fitted$Z_X)

  # Reduced form: Y ~ Z via ddml_ate
  ate_rf <- ddml_ate(
    y = y, D = Z, X = X,
    learners = learners, learners_DX = learners_ZX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_ZX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    stratify = stratify,
    trim = trim, parallel = parallel, silent = silent,
    splits = ate_splits,
    fitted = ate_fitted_rf,
    save_crossval = save_crossval,
    messages = list(start = "", finish = "",
                    y_D1 = messages$y_Z1,
                    y_D0 = messages$y_Z0,
                    D_X = messages$Z_X))

  # First stage: D ~ Z via ddml_ate (reuses splits + propensity)
  ate_fs <- ddml_ate(
    y = D, D = Z, X = X,
    learners = learners_DXZ, learners_DX = learners_ZX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights_DXZ,
    custom_ensemble_weights_DX = custom_ensemble_weights_ZX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    stratify = stratify,
    trim = trim, parallel = parallel, silent = silent,
    splits = ate_rf$splits,
    fitted = list(y_X_D0 = fitted$D_X_Z0,
                  y_X_D1 = fitted$D_X_Z1,
                  D_X = ate_rf$fitted$D_X),
    save_crossval = save_crossval,
    messages = list(start = "", finish = "",
                    y_D1 = messages$D_Z1,
                    y_D0 = messages$D_Z0, D_X = ""))

  ensemble_type <- ate_rf$ensemble_type
  nensb <- ncol(ate_rf$coefficients)

  # Target parameter & influence function --------------------------------------

  late <- as.vector(ate_rf$coefficients) / as.vector(ate_fs$coefficients)

  scores <- array(NA_real_, dim = c(nobs, 1, nensb))
  J <- array(NA_real_, dim = c(1, 1, nensb))
  inf_func <- array(NA_real_, dim = c(nobs, 1, nensb))
  dinf_dtheta <- array(NA_real_, dim = c(nobs, 1, 1, nensb))
  
  for (j in seq_len(nensb)) {
    theta_fs <- ate_fs$coefficients[, j]
    phi_rf <- ate_rf$inf_func[, 1, j]
    phi_fs <- ate_fs$inf_func[, 1, j]
    
    scores[, 1, j] <- phi_rf - phi_fs * late[j]
    J[1, 1, j] <- -theta_fs
    
    J_inv <- csolve(matrix(J[, , j], 1, 1))
    inf_func[, 1, j] <- -matrix(scores[, 1, j], nobs, 1) %*% t(J_inv)
    
    dinf_dtheta[, 1, 1, j] <- phi_fs * J_inv[1, 1]
  }#FOR

  coef_names <- "LATE"
  coef <- matrix(late, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # Output ---------------------------------------------------------------------

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef,
    ensemble_weights = list(
      y_X_Z0 = ate_rf$ensemble_weights$y_X_D0,
      y_X_Z1 = ate_rf$ensemble_weights$y_X_D1,
      D_X_Z0 = ate_fs$ensemble_weights$y_X_D0,
      D_X_Z1 = ate_fs$ensemble_weights$y_X_D1,
      Z_X = ate_rf$ensemble_weights$D_X),
    mspe = list(
      y_X_Z0 = ate_rf$mspe$y_X_D0,
      y_X_Z1 = ate_rf$mspe$y_X_D1,
      D_X_Z0 = ate_fs$mspe$y_X_D0,
      D_X_Z1 = ate_fs$mspe$y_X_D1,
      Z_X = ate_rf$mspe$D_X),
    r2 = list(
      y_X_Z0 = ate_rf$r2$y_X_D0,
      y_X_Z1 = ate_rf$r2$y_X_D1,
      D_X_Z0 = ate_fs$r2$y_X_D0,
      D_X_Z1 = ate_fs$r2$y_X_D1,
      Z_X = ate_rf$r2$D_X),
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Local Average Treatment Effect",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_Z0 = ate_rf$fitted$y_X_D0,
      y_X_Z1 = ate_rf$fitted$y_X_D1,
      D_X_Z0 = ate_fs$fitted$y_X_D0,
      D_X_Z1 = ate_fs$fitted$y_X_D1,
      Z_X = ate_rf$fitted$D_X),
    splits = list(
      y_X_Z0 = ate_rf$splits$y_X_D0,
      y_X_Z1 = ate_rf$splits$y_X_D1,
      D_X_Z0 = ate_fs$splits$y_X_D0,
      D_X_Z1 = ate_fs$splits$y_X_D1,
      Z_X = ate_rf$splits$D_X),
    call = cl,
    subclass = "ddml_late",
    # ddml_late-specific fields
    learners = learners,
    learners_DXZ = learners_DXZ,
    learners_ZX = learners_ZX
  )
}#DDML_LATE
