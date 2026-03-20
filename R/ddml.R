#' Intro to Double/Debiased Machine Learning
#'
#' @name ddml-intro
#'
#' @description All \code{ddml_*} estimators (\code{\link{ddml_plm}},
#' \code{\link{ddml_pliv}}, \code{\link{ddml_fpliv}},
#' \code{\link{ddml_ate}}, \code{\link{ddml_att}},
#' \code{\link{ddml_late}}, \code{\link{ddml_apo}}) return
#' objects that inherit from S3 class \code{"ddml"}.
#'
#' Each object is a list containing the components described below.
#' Estimator-specific fields (e.g., pass-through learner
#' arguments) are documented on the individual estimator pages.
#'
#' The \code{ddml()} constructor can also be used directly to build
#' a \code{"ddml"} object from user-supplied score components,
#' enabling implementation of custom DML estimators that inherit
#' all S3 methods.
#'
#' @details All \code{ddml_*} estimators target a low-dimensional
#' parameter \eqn{\theta_0} identified by a moment condition
#'
#' \deqn{E[m(W; \theta_0, \eta_0)] = 0,}
#'
#' where \eqn{W} denotes observed random variables and
#' \eqn{\eta_0} is a (potentially high-dimensional) nuisance
#' parameter. Throughout, the score \eqn{m} is assumed to be
#' \emph{Neyman orthogonal}.
#'
#' Estimation proceeds via cross-fitting: the sample is randomly
#' partitioned into \eqn{K} folds \eqn{\{I_k\}_{k=1}^K}. For
#' each fold \eqn{k}, nuisance parameters are estimated on the
#' complementary folds (\eqn{\hat\eta_{-k}}) and the scores are
#' evaluated on fold \eqn{k}. The DML estimator
#' \eqn{\hat\theta} solves
#'
#' \deqn{\frac{1}{n} \sum_{k=1}^{K} \sum_{i \in I_k}
#' m(W_i; \hat\theta, \hat\eta_{-k}) = 0.}
#'
#' Inference is based on the influence function. Define the Jacobian
#'
#' \deqn{J(\theta, \eta) = E\!\left[
#'   \frac{\partial m(W; \theta, \eta)}
#'   {\partial \theta'}\right]}
#'
#' and the influence function
#'
#' \deqn{\phi_\theta(W_i; \theta, \eta, J)
#'   = -J^{-1}\,m(W_i; \theta, \eta).}
#'
#' The variance of \eqn{\hat\theta} is then estimated by
#'
#' \deqn{\hat{V} = \frac{1}{n} \sum_i
#'   \phi_\theta(W_i; \hat\theta, \hat\eta_{-k(i)},
#'   \hat{J})\,\phi_\theta(W_i; \hat\theta,
#'   \hat\eta_{-k(i)}, \hat{J})'},
#'
#' where \eqn{\hat{J}} is the sample analog of the Jacobian:
#'
#' \deqn{\hat{J} = \frac{1}{n} \sum_i
#'   \frac{\partial m(W_i; \hat\theta, \hat\eta_{-k(i)})}
#'   {\partial \theta'}.}
#'
#' HC1 and HC3 variance estimators are described in
#' \code{\link{vcov.ral}}. The leverage
#' (see \code{\link{hatvalues.ral}}) for the DML estimator
#' is
#'
#' \deqn{h_\theta(W_i; \theta, \eta, J)
#'   = \mathrm{tr}\!\left(
#'   -J^{-1} \frac{1}{n}
#'   \frac{\partial m(W_i; \theta, \eta)}
#'   {\partial \theta'}\right),}
#'
#' and its sample analog is
#' \eqn{\hat{h}_{\theta,i} = h_\theta(W_i; \hat\theta,
#' \hat\eta_{-k(i)}, \hat{J})}, stored in
#' \code{dinf_dtheta}.
#'
#' Under regularity conditions and sufficient convergence of
#' \eqn{\hat\eta}, the DML estimator is asymptotically normal:
#'
#' \deqn{\sqrt{n}\,\hat{V}^{-1/2}(\hat\theta - \theta_0)
#' \overset{d}{\to} N(0, I).}
#'
#' Further details and regularity conditions are given in
#' Chernozhukov et al. (2018). The specific forms of the
#' score \eqn{m} and Jacobian \eqn{J} for each estimator
#' are documented on their respective help pages (e.g.,
#' \code{\link{ddml_plm}}, \code{\link{ddml_ate}}).
#'
#' @references
#' Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E,
#' Wiemann T (2026). "An Introduction to Double/Debiased Machine
#' Learning." Journal of Economic Literature, forthcoming.
#'
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B,
#' Newey W, Robins J (2018). "Double/debiased machine learning
#' for treatment and structural parameters." The Econometrics
#' Journal, 21(1), C1-C68.
#'
#' @section Common output components:
#' \describe{
#' \item{\code{coefficients}}{A matrix of estimated target
#'     parameters: rows correspond to components of
#'     \eqn{\theta}, columns to ensemble types.}
#' \item{\code{ensemble_weights}}{A named list. Each element
#'     is a weight matrix (or 3D array when
#'     \code{shortstack = TRUE}) showing the weight assigned
#'     to each base learner by the ensemble procedure for
#'     the corresponding nuisance equation.}
#' \item{\code{mspe}}{A named list of numeric vectors
#'     containing per-learner out-of-sample MSPEs, computed
#'     from cross-fitted residuals.}
#' \item{\code{r2}}{A named list of numeric vectors
#'     containing per-learner out-of-sample R-squared values.}
#' \item{\code{inf_func}}{A 3D array of evaluated influence
#'     functions (\code{n x p x nensb}).}
#' \item{\code{dinf_dtheta}}{An optional 4D array of dimension
#'     \code{(n x p x p x nensb)} containing the derivatives of the
#'     influence functions with respect to \eqn{\theta}. Used
#'     internally by \code{\link{hatvalues.ral}} for HC3
#'     inference.}
#' \item{\code{scores}}{A 3D array of evaluated Neyman
#'     orthogonal scores (\code{n x p x nensb}).}
#' \item{\code{J}}{A 3D array of evaluated Jacobians
#'     (\code{p x p x nensb}).}
#' \item{\code{fitted}}{A named list of per-equation
#'     cross-fitted prediction objects. Can be passed back
#'     via the \code{fitted} argument together with
#'     \code{splits} to skip cross-fitting on
#'     re-estimation.}
#' \item{\code{splits}}{The data splitting structure
#'     (subsamples, CV subsamples, and any stratification
#'     indices).}
#' \item{\code{ensemble_type}}{Character vector of ensemble
#'     types used.}
#' \item{\code{cluster_variable}}{The cluster variable
#'     vector used for sample splitting and inference.}
#' \item{\code{nobs}}{Number of observations.}
#' \item{\code{sample_folds}}{Number of cross-fitting folds.}
#' \item{\code{shortstack}}{Logical indicating whether
#'     short-stacking was used.}
#' \item{\code{call}}{The matched call.}
#' \item{\code{coef_names}}{Character vector of coefficient
#'     names.}
#' \item{\code{estimator_name}}{Character string identifying
#'     the estimator (e.g., \code{"Partially Linear Model"}).}
#' }
#'
#' @section S3 methods:
#' The following generic methods are available for all
#' \code{ddml} objects: \code{\link{summary.ddml}},
#' \code{\link{coef.ral}}, \code{\link{vcov.ral}},
#' \code{\link{confint.ral}}, \code{\link{hatvalues.ral}},
#' \code{\link{nobs.ral}}, \code{\link{tidy.ddml}},
#' \code{\link{glance.ddml}}, and
#' \code{\link{diagnostics}}.
#'
#' @param y The outcome variable.
#' @param D A matrix of endogenous variables.
#' @param X A (sparse) matrix of control variables.
#' @param learners May take one of two forms, depending on whether a
#' single learner or stacking with multiple learners is used for
#' estimation of the conditional expectation functions.
#' If a single learner is used, \code{learners} is a list with
#' two named elements:
#' \itemize{
#'     \item{\code{what} The base learner function. The function
#'         must be such that it predicts a named input \code{y}
#'         using a named input \code{X}.}
#'     \item{\code{args} Optional arguments to be passed to
#'         \code{what}.}
#' }
#' If stacking with multiple learners is used, \code{learners} is
#' a list of lists, each containing three named elements:
#' \itemize{
#'     \item{\code{what} The base learner function. The function
#'         must be such that it predicts a named input \code{y}
#'         using a named input \code{X}.}
#'     \item{\code{args} Optional arguments to be passed to
#'         \code{what}.}
#'     \item{\code{assign_X} An optional vector of column indices
#'         corresponding to control variables in \code{X} that
#'         are passed to the base learner.}
#' }
#' Omission of the \code{args} element results in default
#' arguments being used in \code{what}. Omission of
#' \code{assign_X} results in inclusion of all variables in
#' \code{X}.
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into
#' final estimate of the conditional expectation functions.
#' Possible values are:
#' \itemize{
#'     \item{\code{"nnls"} Non-negative least squares.}
#'     \item{\code{"nnls1"} Non-negative least squares with the
#'         constraint that all weights sum to one.}
#'     \item{\code{"singlebest"} Select base learner with minimum
#'         MSPE.}
#'     \item{\code{"ols"} Ordinary least squares.}
#'     \item{\code{"average"} Simple average over base learners.}
#' }
#' Multiple ensemble types may be passed as a vector of strings.
#' @param shortstack Boolean to use short-stacking.
#' @param cv_folds Number of folds used for cross-validation in
#' ensemble construction.
#' @param custom_ensemble_weights A numerical matrix with
#' user-specified ensemble weights. Each column corresponds to a
#' custom ensemble specification, each row corresponds to a base
#' learner in \code{learners} (in chronological order). Optional
#' column names are used to name the estimation results
#' corresponding the custom ensemble specification.
#' @param cluster_variable A vector of cluster indices.
#' @param silent Boolean to silence estimation updates.
#' @param parallel An optional named list with parallel processing
#' options. When \code{NULL} (the default), computation is
#' sequential. Supported fields:
#' \describe{
#'     \item{\code{cores}}{Number of cores to use.}
#'     \item{\code{export}}{Character vector of object names to
#'         export to parallel workers (for custom learners that
#'         reference global objects).}
#'     \item{\code{packages}}{Character vector of additional
#'         package names to load on workers (for custom learners
#'         that use packages not imported by \code{ddml}).}
#' }
#' @param fitted An optional named list of per-equation cross-fitted
#' predictions, typically obtained from a previous fit via
#' \code{fit$fitted}. When supplied (together with \code{splits}),
#' base learners are not re-fitted; only ensemble weights are
#' recomputed. This allows fast re-estimation with a different
#' \code{ensemble_type}. See \code{\link{ddml_plm}} for
#' an example.
#' @param splits An optional list of sample split objects, typically
#' obtained from a previous fit via \code{fit$splits}. Must be
#' supplied when \code{fitted} is provided. Can also be used
#' standalone to provide pre-computed sample folds.
#' @param save_crossval Logical indicating whether to store the inner
#' cross-validation residuals used for ensemble weight
#' computation. Default \code{TRUE}. When \code{TRUE}, subsequent
#' pass-through calls with data-driven ensembles (e.g.,
#' \code{"nnls"}) reproduce per-fold weights exactly. Set to
#' \code{FALSE} to reduce object size at the cost of approximate
#' weight recomputation.
#' @param ... Additional arguments passed to internal methods.
#'
#' @family ddml estimators
NULL

#' Construct a \code{ddml} Object.
#'
#' @family utilities
#'
#' @description Build a \code{"ddml"} object from user-supplied score
#' components. The resulting object inherits all S3 methods
#' available for \code{ddml} objects, including
#' \code{\link{summary.ddml}}, \code{\link{confint.ral}},
#' \code{\link{vcov.ral}}, and \code{\link{tidy.ddml}}.
#'
#' @param coefficients A \code{(p x nensb)} matrix of estimated
#' target parameters. Rows correspond to components of
#' \eqn{\theta}, columns to ensemble types.
#' @param scores A 3D array of evaluated Neyman orthogonal scores
#' with dimensions \code{(n x p x nensb)}.
#' @param J A 3D array of evaluated Jacobians with dimensions
#' \code{(p x p x nensb)}.
#' @param inf_func A 3D array of evaluated influence functions
#' with dimensions \code{(n x p x nensb)}.
#' @param nobs Number of observations.
#' @param coef_names Character vector of coefficient names
#' (length \code{p}).
#' @param estimator_name Character string identifying the estimator
#' (e.g., \code{"My Custom Estimator"}).
#' @param ensemble_type Character vector of ensemble types. Defaults
#' to \code{colnames(coefficients)}.
#' @param cluster_variable A vector of cluster indices. Defaults to
#' \code{seq_len(nobs)}.
#' @param sample_folds Number of cross-fitting folds used. Optional.
#' @param cv_folds Number of cross-validation folds used. Optional.
#' @param shortstack Logical indicating whether short-stacking was
#' used. Default \code{FALSE}.
#' @param ensemble_weights A named list of ensemble weight matrices.
#' Optional.
#' @param mspe A named list of per-learner MSPEs. Optional.
#' @param r2 A named list of per-learner R-squared values. Optional.
#' @param fitted A named list of per-equation cross-fitted prediction
#' objects. Optional.
#' @param splits A list of sample split objects. Optional.
#' @param call The matched call. Defaults to \code{match.call()}.
#' @param subclass Optional character string for a subclass name. If
#' provided, the object will have class
#' \code{c(subclass, "ddml")}.
#' @param dinf_dtheta An optional 4D array of dimensions \code{(nobs x p x p x nensb)}
#'     containing the derivatives of the influence functions.
#' @param ... Additional named components to include in the object.
#'
#' @return An object of S3 class \code{"ddml"} (or
#' \code{c(subclass, "ddml")} if \code{subclass} is specified).
#' See \code{\link{ddml-intro}} for the output structure.
#'
#' @export
#'
#' @examples
#' # A minimal example: construct a ddml object from pre-computed
#' #     score components for a simple mean estimator.
#' n <- 100
#' y <- rnorm(n)
#' theta <- mean(y)
#'
#' scores <- array(y - theta, dim = c(n, 1, 1))
#' J <- array(-1, dim = c(1, 1, 1))
#' psi_b <- list(matrix(y, ncol = 1))
#' psi_a <- list(array(-1, dim = c(n, 1, 1)))
#' inf_func <- array(y - theta, dim = c(n, 1, 1))
#' dinf_dtheta <- array(1, dim = c(n, 1, 1, 1))
#' coef <- matrix(theta, 1, 1, dimnames = list("mean", "custom"))
#'
#' fit <- ddml(coefficients = coef, scores = scores, J = J,
#'         inf_func = inf_func, nobs = n, coef_names = "mean",
#'         dinf_dtheta = dinf_dtheta,
#'         estimator_name = "Sample Mean")
#' summary(fit)
ddml <- function(coefficients, scores, J, inf_func,
                 nobs, coef_names, estimator_name,
                 ensemble_type = colnames(coefficients),
                 cluster_variable = seq_len(nobs),
                 sample_folds = NULL,
                 cv_folds = NULL,
                 shortstack = FALSE,
                 ensemble_weights = NULL,
                 mspe = NULL, r2 = NULL,
                 fitted = NULL, splits = NULL,
                 call = match.call(),
                 subclass = NULL,
                 dinf_dtheta = NULL,
                 ...) {
  # Validate required fields
  if (!is.matrix(coefficients)) {
    stop("'coefficients' must be a matrix.", call. = FALSE)
  }#IF
  if (length(dim(scores)) != 3) {
    stop("'scores' must be a 3D array.", call. = FALSE)
  }#IF
  if (length(dim(J)) != 3) {
    stop("'J' must be a 3D array.", call. = FALSE)
  }#IF
  if (!is.numeric(inf_func) || length(dim(inf_func)) != 3) {
    stop("'inf_func' must be a 3D numeric array.", call. = FALSE)
  }#IF
  if (!is.null(dinf_dtheta) && !is.array(dinf_dtheta)) {
    stop("'dinf_dtheta' must be a 4D array or NULL.", call. = FALSE)
  }#IF

  # Dimension consistency
  p <- nrow(coefficients)
  nensb <- ncol(coefficients)
  if (dim(scores)[1] != nobs || dim(scores)[2] != p ||
      dim(scores)[3] != nensb) {
    stop("'scores' dimensions must be (nobs x p x nensb).",
         call. = FALSE)
  }#IF
  if (dim(J)[1] != p || dim(J)[2] != p ||
      dim(J)[3] != nensb) {
    stop("'J' dimensions must be (p x p x nensb).", call. = FALSE)
  }#IF
  if (dim(inf_func)[1] != nobs || dim(inf_func)[2] != p ||
      dim(inf_func)[3] != nensb) {
    stop("'inf_func' dimensions must be (nobs x p x nensb).", call. = FALSE)
  }#IF
  if (!is.null(dinf_dtheta)) {
    if (length(dim(dinf_dtheta)) != 4 ||
        dim(dinf_dtheta)[1] != nobs || dim(dinf_dtheta)[2] != p ||
        dim(dinf_dtheta)[3] != p || dim(dinf_dtheta)[4] != nensb) {
      stop("'dinf_dtheta' dimensions must be (nobs x p x p x nensb).", call. = FALSE)
    }#IF
  }#IF

  # Assemble the object
  obj <- c(list(
    coefficients = coefficients,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    inf_func = inf_func,
    dinf_dtheta = dinf_dtheta,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = estimator_name,
    ensemble_type = ensemble_type,
    nobs = nobs,
    nfit = ncol(coefficients),
    fit_labels = colnames(coefficients),
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = fitted,
    splits = splits,
    call = call), list(...))

  cls <- if (!is.null(subclass)) c(subclass, "ddml", "ral") else
    c("ddml", "ral")
  class(obj) <- cls
  obj
}#DDML

# S3 methods ================================================================
# Inference methods (coef, nobs, vcov, confint, hatvalues, tidy, glance)
# are inherited from the "ral" superclass in ral.R.

#' Subscript a summary.ddml object (deprecated).
#' @param x An object of class \code{summary.ddml}.
#' @param ... Indices passed to \code{[}.
#' @keywords internal
#' @export
`[.summary.ddml` <- function(x, ...) {
  message("Note: subscripting a summary.ddml object with ",
          "'[' is deprecated. Use x$coefficients[...] instead.")
  x$coefficients[...]
}#`[.SUMMARY.DDML`

#' Summary for DDML Estimators
#'
#' @description Computes a coefficient table with estimates,
#' standard errors, z-values, and p-values for all
#' ensemble types. Standard errors are based on a
#' heteroskedasticity-robust sandwich variance; see
#' \code{\link{vcov.ral}} for the HC0/HC1/HC3 formulas.
#'
#' @param object An object of class \code{ddml}.
#' @param type Character. HC type (\code{"HC0"},
#'     \code{"HC1"}, or \code{"HC3"}). Default
#'     \code{"HC1"}.
#' @param ... Currently unused.
#'
#' @return An object of class \code{summary.ddml} with:
#' \describe{
#'     \item{\code{coefficients}}{A 3-dimensional array
#'         (\eqn{p \times 4 \times}{p x 4 x} nensb) of estimates, standard
#'         errors, z-values, and p-values.}
#'     \item{\code{type}}{The HC type used.}
#'     \item{\code{nobs}}{Number of observations.}
#'     \item{\code{sample_folds}}{Number of cross-fitting
#'         folds.}
#'     \item{\code{ensemble_type}}{Ensemble type labels.}
#' }
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' summary(plm_fit)
#' summary(plm_fit, type = "HC3")
#' }
#'
#' @seealso \code{\link{vcov.ral}}
#'
#' @family ddml inference
#' @method summary ddml
#' @export
summary.ddml <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))

  # DML-specific fit labels
  single_learner <- is_single_learner(object$learners)
  ens_type <- if (single_learner) "single base learner" else
    object$ensemble_type
  object$fit_labels <- ens_type

  # Delegate table computation to ral
  result <- summary.ral(object, type = type, ...)

  # Attach DML-specific fields
  result$model_type    <- class(object)[1]
  result$sample_folds  <- object$sample_folds
  result$shortstack    <- object$shortstack
  result$ensemble_type <- ens_type

  class(result) <- c(
    paste0("summary.", class(object)[1]),
    "summary.ddml")
  result
}#SUMMARY.DDML

#' @rdname summary.ddml
#'
#' @param x An object of class \code{summary.ddml}.
#' @param digits Number of significant digits. Default 3.
#'
#' @method print summary.ddml
#' @export
print.summary.ddml <- function(x, digits = 3, ...) {
  model_name <- x$estimator_name
  if (is.null(model_name)) model_name <- x$model_type

  cat("DDML estimation:", model_name, "\n")
  cat("Obs:", x$nobs,
      "  Folds:", x$sample_folds)
  if (!is.null(x$shortstack) && x$shortstack) {
    cat("  Stacking: short-stack")
  }#IF
  if (!is.null(x$type) && x$type != "HC1") {
    cat("  SE:", x$type)
  }#IF
  cat("\n\n")

  print_coef_tables(x$coefficients, fit_label = "Ensemble type",
                    digits = digits)

  invisible(x)
}#PRINT.SUMMARY.DDML

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a DDML Object
#'
#' DML-specific tidy method. Adds \code{ensemble_type}
#' column labeling based on the estimator's learner/ensemble
#' configuration. Delegates to \code{tidy.ral} for the base
#' table computation.
#'
#' @param x A \code{ddml} object.
#' @param ensemble_idx Integer index of the ensemble type to
#'     report. Defaults to 1. Set to \code{NULL} for all.
#' @param conf.int Logical. Include confidence intervals?
#'     Default \code{FALSE}.
#' @param conf.level Confidence level. Default 0.95.
#' @param type Character. HC type. Default \code{"HC1"}.
#' @param uniform Logical. Uniform CIs? Default \code{FALSE}.
#' @param bootstraps Integer. Bootstrap draws. Default 999.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, and \code{ensemble_type}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' tidy(plm_fit)
#' tidy(plm_fit, conf.int = TRUE)
#' }
#'
#' @export
#' @method tidy ddml
tidy.ddml <- function(x, ensemble_idx = 1, conf.int = FALSE,
                      conf.level = 0.95,
                      type = "HC1",
                      uniform = FALSE,
                      bootstraps = 999L, ...) {
  res <- tidy.ral(x, fit_idx = ensemble_idx,
                  conf.int = conf.int,
                  conf.level = conf.level,
                  type = type, uniform = uniform,
                  bootstraps = bootstraps)
  # Rename fit_label -> ensemble_type for DML compatibility
  names(res)[names(res) == "fit_label"] <- "ensemble_type"
  res
}#TIDY.DDML

#' Glance at a DDML Object
#'
#' DML-specific glance method. Includes DML fields like
#' \code{sample_folds}, \code{shortstack}, and
#' \code{model_type}.
#'
#' @param x A \code{ddml} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' glance(plm_fit)
#' }
#'
#' @export
#' @method glance ddml
glance.ddml <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    sample_folds = x$sample_folds,
    shortstack = if (is.null(x$shortstack)) FALSE else
      x$shortstack,
    ensemble_type = paste(x$ensemble_type, collapse = ", "),
    model_type = class(x)[1],
    estimator_name = if (is.null(x$estimator_name)) {
      class(x)[1]
    } else {
      x$estimator_name
    },
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML
