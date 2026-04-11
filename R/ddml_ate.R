#' Estimator for the Average Treatment Effect
#'
#' @family ddml estimators
#'
#' @description Estimator for the average treatment effect and the average
#'     treatment effect on the treated.
#'
#' @details
#' \strong{Parameter of Interest:} \code{ddml_ate} and \code{ddml_att} provide
#'     Double/Debiased Machine Learning estimators for the average treatment
#'     effect and the average treatment effect on the treated, respectively.
#'     Under conditional unconfoundedness and overlap, the parameters
#'     are identified by the following reduced form parameters:
#'
#' \deqn{\theta_0^{\textrm{ATE}} = E[E[Y|D=1, X] - E[Y|D=0, X]]}
#'
#' and the average treatment effect on the treated (ATT) is defined as
#'
#' \deqn{\theta_0^{\textrm{ATT}} = E[Y|D=1] - E[E[Y|D=0, X]|D = 1].}
#'
#' where \eqn{W \equiv (Y, D, X)} is the observed random vector.
#'
#' \strong{Neyman Orthogonal Score:} The Neyman orthogonal scores are:
#'
#' \deqn{m^{\textrm{ATE}}(W; \theta, \eta) = \frac{D(Y - \ell_1(X))}{r(X)} - \frac{(1-D)(Y-\ell_0(X))}{1-r(X)} + \ell_1(X) - \ell_0(X) - \theta}
#'
#' \deqn{m^{\textrm{ATT}}(W; \theta, \eta) = \frac{D(Y - \ell_0(X))}{\pi} - \frac{r(X)(1-D)(Y-\ell_0(X))}{\pi(1-r(X))} - \frac{D}{\pi}\theta}
#'
#' where the nuisance parameters are \eqn{\eta = (\ell_0, \ell_1, r, \pi)} taking true values
#'     \eqn{\ell_{d,0}(X) = E[Y|D=d, X]}, \eqn{r_0(X) = \Pr(D=1|X)}, and \eqn{\pi_0 = \Pr(D=1)}.
#'
#' \strong{Jacobian:}
#'
#' \deqn{J^{\textrm{ATE}} = -1}
#'
#' \deqn{J^{\textrm{ATT}} = -1}
#'
#' See \code{\link{ddml-intro}} for how the influence function
#' and inference are derived from these components.
#'
#' @inheritParams ddml-intro
#' @inheritParams ddml_apo
#' @param D The binary endogenous variable of interest.
#' @param splits An optional list of sample split objects. For
#'     \code{ddml_ate}/\code{ddml_att}, recommended keys are
#'     \code{subsamples}, \code{subsamples_byD}, \code{cv_subsamples},
#'     and \code{cv_subsamples_byD}.
#' @param ... Additional arguments passed to internal methods.
#'
#' @return \code{ddml_ate} and \code{ddml_att} return objects of S3
#'     class \code{ddml_ate}/\code{ddml_att} and \code{ddml}. See
#'     \code{\link{ddml-intro}} for the common output structure.
#'     Additional pass-through fields: \code{learners},
#'     \code{learners_DX}.
#' @export
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
#'                     learners = list(list(what = ols),
#'                                     list(what = mdl_glmnet),
#'                                     list(what = mdl_glmnet,
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
  messages <- resolve_messages(dots, "ddml_ate", list(
    y_D0 = "E[Y|D=0,X]",
    y_D1 = "E[Y|D=1,X]",
    D_X = "E[D|X]"))

  validate_inputs(y = y, D = D, X = X, learners = learners,
                  sample_folds = sample_folds,
                  cv_folds = cv_folds,
                  ensemble_type = ensemble_type, trim = trim,
                  cluster_variable = cluster_variable,
                  require_binary_D = TRUE)
  validate_custom_weights(custom_ensemble_weights, learners)
  validate_custom_weights(custom_ensemble_weights_DX, learners_DX)

  nobs <- length(y)

  validate_fitted_splits_pair(fitted, splits)

  indxs <- get_sample_splits(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    D = D, stratify = stratify,
    subsamples = splits$D_X$subsamples,
    subsamples_byD = list(
      splits$y_X_D0$subsamples,
      splits$y_X_D1$subsamples),
    cv_subsamples = splits$D_X$cv_subsamples,
    cv_subsamples_byD = list(
      splits$y_X_D0$cv_subsamples,
      splits$y_X_D1$cv_subsamples))
  check_subsamples(indxs$subsamples, indxs$subsamples_byD,
                   stratify, D)

  t0 <- proc.time()[3]
  announce_start(messages, parallel, silent)

  # Reduced-form estimation ----------------------------------------------------

  shared_splits <- list(subsamples = indxs$subsamples,
                        cv_subsamples = indxs$cv_subsamples)
  apo_splits_1 <- list(
    y_X = list(subsamples = indxs$subsamples_byD[[2]],
               cv_subsamples = indxs$cv_subsamples_byD[[2]]),
    D_X = shared_splits)
  apo_splits_0 <- list(
    y_X = list(subsamples = indxs$subsamples_byD[[1]],
               cv_subsamples = indxs$cv_subsamples_byD[[1]]),
    D_X = shared_splits)

  # E[g(1,X)] via ddml_apo (also estimates P(D=1|X) internally)
  apo_1 <- ddml_apo(
    y = y, D = D, X = X, d = 1, weights = NULL,
    learners = learners, learners_DX = learners_DX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_DX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = apo_splits_1,
    fitted = list(y_X = fitted$y_X_D1,
                  D_X = fitted$D_X),
    messages = list(start = "", finish = "",
                    y_X = messages$y_D1,
                    D_X = messages$D_X))

  # Flip propensity: P(D=0|X) = 1 - P(D=1|X).
  fitted_D_X_0 <- list(cf_fitted = 1 - apo_1$fitted$D_X$cf_fitted)

  # E[g(0,X)] via ddml_apo
  apo_0 <- ddml_apo(
    y = y, D = D, X = X, d = 0, weights = NULL,
    learners = learners, learners_DX = learners_DX,
    sample_folds = sample_folds, cv_folds = cv_folds,
    custom_ensemble_weights = custom_ensemble_weights,
    custom_ensemble_weights_DX = custom_ensemble_weights_DX,
    cluster_variable = cluster_variable,
    ensemble_type = ensemble_type, shortstack = shortstack,
    trim = trim, parallel = parallel, silent = silent,
    splits = apo_splits_0,
    fitted = list(y_X = fitted$y_X_D0,
                  D_X = fitted_D_X_0),
    messages = list(start = "", finish = "",
                    y_X = messages$y_D0, D_X = ""))

  ensemble_type <- apo_1$ensemble_type
  nensb <- ncol(apo_1$coefficients)

  # Target parameter & influence function --------------------------------------

  ate <- as.vector(apo_1$coefficients) - as.vector(apo_0$coefficients)

  scores <- array(NA_real_, dim = c(nobs, 1, nensb))
  J <- array(NA_real_, dim = c(1, 1, nensb))
  inf_func <- array(NA_real_, dim = c(nobs, 1, nensb))
  dinf_dtheta <- array(NA_real_, dim = c(nobs, 1, 1, nensb))
  for (j in seq_len(nensb)) {
    inf_func[, 1, j] <- apo_1$inf_func[, 1, j] - apo_0$inf_func[, 1, j]
    scores[, 1, j] <- inf_func[, 1, j]
    J[1, 1, j] <- -1
    dinf_dtheta[, 1, 1, j] <- -1
  }#FOR
  coef_names <- "ATE"
  coef <- matrix(ate, nrow = 1, ncol = nensb)
  rownames(coef) <- coef_names
  colnames(coef) <- ensemble_type

  # Output ---------------------------------------------------------------------

  announce_finish(t0, messages, silent)

  ddml(
    coefficients = coef,
    ensemble_weights = list(
      y_X_D0 = apo_0$ensemble_weights$y_X,
      y_X_D1 = apo_1$ensemble_weights$y_X,
      D_X = apo_1$ensemble_weights$D_X),
    mspe = list(y_X_D0 = apo_0$mspe$y_X,
                y_X_D1 = apo_1$mspe$y_X,
                D_X = apo_1$mspe$D_X),
    r2 = list(y_X_D0 = apo_0$r2$y_X,
              y_X_D1 = apo_1$r2$y_X,
              D_X = apo_1$r2$D_X),
    inf_func = inf_func, dinf_dtheta = dinf_dtheta,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = "Average Treatment Effect",
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = if (shortstack) NULL else cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = list(
      y_X_D0 = apo_0$fitted$y_X,
      y_X_D1 = apo_1$fitted$y_X,
      D_X = apo_1$fitted$D_X),
    splits = list(
      y_X_D0 = list(
        subsamples = indxs$subsamples_byD[[1]],
        cv_subsamples = indxs$cv_subsamples_byD[[1]]),
      y_X_D1 = list(
        subsamples = indxs$subsamples_byD[[2]],
        cv_subsamples = indxs$cv_subsamples_byD[[2]]),
      D_X = list(
        subsamples = indxs$subsamples,
        cv_subsamples = indxs$cv_subsamples)),
    call = cl,
    subclass = "ddml_ate",
    # ddml_ate-specific fields
    learners = learners,
    learners_DX = learners_DX
  )
}#DDML_ATE
