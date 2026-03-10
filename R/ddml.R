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
#' parameter, typically a vector of conditional expectation
#' functions. The score \eqn{m} is \emph{Neyman orthogonal},
#' meaning the moment condition is locally insensitive to
#' perturbations of \eqn{\eta} around \eqn{\eta_0}.
#'
#' Each estimator in \code{ddml} uses a score that decomposes
#' linearly in \eqn{\theta}:
#'
#' \deqn{m(W_i; \theta, \eta) = \psi_b(W_i; \eta) + \psi_a(W_i; \eta)\,\theta.}
#'
#' See the individual estimator pages for the specific forms
#' of \eqn{\psi_a} and \eqn{\psi_b}.
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
#' For the linear scores used in \code{ddml}, this yields the
#' closed-form solution
#'
#' \deqn{\hat\theta = -\hat{J}^{-1}\,\hat\psi_b,}
#'
#' where \eqn{\hat\psi_b = n^{-1}\sum_{k}\sum_{i \in I_k}
#' \psi_b(W_i; \hat\eta_{-k})} and
#' \eqn{\hat{J} = n^{-1}\sum_k \sum_{i \in I_k}
#' \psi_a(W_i; \hat\eta_{-k})} is the sample Jacobian.
#'
#' Inference is based on the sandwich variance estimator
#'
#' \deqn{\hat\Sigma = \hat{J}^{-1} \left(\frac{1}{n}
#' \sum_i m_i m_i^\top\right)
#' \hat{J}^{-\top} / n}
#'
#' where \eqn{m_i = m(W_i; \hat\theta, \hat\eta_{-k(i)})};
#' see \code{\link{vcov.ddml}} for HC0/HC1/HC3 variants.
#' Under regularity conditions and sufficient convergence of \eqn{\hat\eta}, 
#' the DML estimator is asymptotically normal:
#'
#' \deqn{\sqrt{n}\,\hat\Sigma^{-1/2}(\hat\theta - \theta_0)
#' \overset{d}{\to} N(0, I).}
#'
#' Further details and regularity conditions are given in
#' Chernozhukov et al. (2018).
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
#' \item{\code{psi_a}, \code{psi_b}}{Score component lists
#'     (length \code{nensb}). \code{psi_a[[j]]} is an
#'     \code{(n x p x p)} array and \code{psi_b[[j]]} is an
#'     \code{(n x p)} matrix. Used internally by
#'     \code{\link{hatvalues.ddml}}.}
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
#' \code{\link{coef.ddml}}, \code{\link{vcov.ddml}},
#' \code{\link{confint.ddml}}, \code{\link{hatvalues.ddml}},
#' \code{\link{nobs.ddml}}, \code{\link{tidy.ddml}},
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
#' \code{\link{summary.ddml}}, \code{\link{confint.ddml}},
#' \code{\link{vcov.ddml}}, and \code{\link{tidy.ddml}}.
#'
#' @param coefficients A \code{(p x nensb)} matrix of estimated
#' target parameters. Rows correspond to components of
#' \eqn{\theta}, columns to ensemble types.
#' @param scores A 3D array of evaluated Neyman orthogonal scores
#' with dimensions \code{(n x p x nensb)}.
#' @param J A 3D array of evaluated Jacobians with dimensions
#' \code{(p x p x nensb)}.
#' @param psi_a A list of length \code{nensb}. Each element is an
#' \code{(n x p x p)} array of individual-level Jacobian
#' contributions.
#' @param psi_b A list of length \code{nensb}. Each element is an
#' \code{(n x p)} matrix of individual-level score
#' contributions.
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
#' coef <- matrix(theta, 1, 1,
#'            dimnames = list("mean", "custom"))
#'
#' fit <- ddml(coefficients = coef, scores = scores, J = J,
#'         psi_a = psi_a, psi_b = psi_b,
#'         nobs = n, coef_names = "mean",
#'         estimator_name = "Sample Mean")
#' summary(fit)
ddml <- function(coefficients, scores, J,
                 psi_a, psi_b,
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
  if (!is.list(psi_a) || !is.list(psi_b)) {
    stop("'psi_a' and 'psi_b' must be lists.", call. = FALSE)
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
    stop("'J' dimensions must be (p x p x nensb).",
         call. = FALSE)
  }#IF
  if (length(psi_a) != nensb || length(psi_b) != nensb) {
    stop("'psi_a' and 'psi_b' must have length nensb.",
         call. = FALSE)
  }#IF

  # Assemble the object
  obj <- c(list(
    coefficients = coefficients,
    ensemble_weights = ensemble_weights,
    mspe = mspe,
    r2 = r2,
    psi_a = psi_a, psi_b = psi_b,
    scores = scores, J = J,
    coef_names = coef_names,
    estimator_name = estimator_name,
    ensemble_type = ensemble_type,
    nobs = nobs,
    sample_folds = sample_folds,
    cv_folds = cv_folds,
    shortstack = shortstack,
    cluster_variable = cluster_variable,
    fitted = fitted,
    splits = splits,
    call = call), list(...))

  cls <- if (!is.null(subclass)) c(subclass, "ddml") else "ddml"
  class(obj) <- cls
  obj
}#DDML

# == S3 Methods =============================================================

#' @rdname ddml
#' @param x An object of class \code{ddml}.
#' @param name Name of the element to extract.
#' @export
`$.ddml` <- function(x, name) {
  if (name == "weights") {
    message("Note: '$weights' is deprecated for ddml ",
            "objects. Use '$ensemble_weights' instead.")
    return(.subset2(x, "ensemble_weights"))
  }#IF
  .subset2(x, name)
}#`$.DDML`

#' Extract Model Coefficients
#'
#' @description Extracts the estimated coefficients
#' from a DDML model for the specified or default
#' ensemble type.
#'
#' @param object An object of class \code{ddml}.
#' @param ... Currently unused.
#'
#' @return Named vector (single ensemble) or matrix
#' (multiple ensembles).
#'
#' @examples
#' \donttest{
#' # Fit a PLM and extract coefficients
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' coef(plm_fit)
#' }
#'
#' @family ddml inference
#' @method coef ddml
#' @export
coef.ddml <- function(object, ...) {
  cf <- object$coefficients
  if (is.matrix(cf) && ncol(cf) == 1) {
    nm <- rownames(cf)
    cf <- as.vector(cf)
    names(cf) <- nm
  }#IF
  cf
}#COEF.DDML

#' Extract Number of Observations
#'
#' @description Returns the number of observations
#' used to fit the DDML model.
#'
#' @param object An object of class \code{ddml}.
#' @param ... Currently unused.
#'
#' @return An integer specifying the number of observations.
#'
#' @method nobs ddml
#' @importFrom stats nobs
#' @export
nobs.ddml <- function(object, ...) {
  object$nobs
}#NOBS.DDML

#' Extract Generalized Leverage (Hat Values)
#'
#' @description Computes the generalized leverage (hat values) for a DDML 
#' estimator. These values are used internally to compute
#' heteroskedasticity-robust HC3 standard errors.
#'
#' @details See \code{\link{ddml-intro}} for the DML framework and
#' the definition of \eqn{\psi_a} and \eqn{\hat{J}}.
#' For the linear scores used in \code{ddml},
#' the generalized leverage simplifies to
#'
#' \eqn{h_{ii} = \mathrm{tr}(\psi_a(W_i;
#' \hat\eta)\,(n\hat{J})^{-1}).}
#'
#' @param model An object of class \code{ddml}.
#' @param ensemble_idx Integer index of the ensemble type to extract leverage
#' values for. Defaults to 1.
#' @param ... Currently unused.
#'
#' @return A numeric vector of generalized leverage values.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' h = hatvalues(plm_fit)
#' head(h)
#' }
#' 
#' @seealso \code{\link{vcov.ddml}} for the use of leverage in HC3 standard errors.
#'
#' @family ddml inference
#' @importFrom stats hatvalues
#' @export
hatvalues.ddml <- function(model, ensemble_idx = 1, ...) {
  validate_method_args(model, ensemble_idx = ensemble_idx)
  psi_a_j <- model$psi_a[[ensemble_idx]]
  J_j <- model$J[, , ensemble_idx, drop = FALSE]
  dim(J_j) <- dim(J_j)[1:2]
  J_inv <- csolve(J_j)
  p <- ncol(J_j)
  n <- model$nobs
  dim(psi_a_j) <- c(n, p * p)
  as.vector((psi_a_j %*% as.vector(t(J_inv))) / n)
}#HATVALUES.DDML

#' Variance-Covariance Matrix for DDML Estimators
#'
#' @description Computes a heteroskedasticity-robust
#' variance-covariance matrix for the DDML estimator
#' \eqn{\hat\theta}.
#'
#' @details See \code{\link{ddml-intro}} for the DML framework,
#' including the definitions of the score \eqn{m_i},
#' the Jacobian \eqn{\hat{J}}, and the base sandwich
#' estimator \eqn{\hat\Sigma}. This function provides
#' three variants:
#'
#' \strong{HC0}:
#' \deqn{V_{\mathrm{HC0}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i m_i m_i'\right)
#'   \hat{J}^{-\top} / n}
#'
#' \strong{HC1} (default):
#' \deqn{V_{\mathrm{HC1}} = V_{\mathrm{HC0}}
#'   \times \frac{n}{n - p}}
#'
#' where \eqn{p} is the dimension of \eqn{\theta}.
#'
#' \strong{HC3}:
#' \deqn{V_{\mathrm{HC3}} = \hat{J}^{-1}
#'   \left(\frac{1}{n}\sum_i
#'   \frac{m_i m_i'}{(1 - h_{ii})^2}\right)
#'   \hat{J}^{-\top} / n}
#'
#' where \eqn{h_{ii}} is the generalized leverage;
#' see \code{\link{hatvalues.ddml}}.
#'
#' @param object An object of class \code{ddml}.
#' @param ensemble_idx Integer index of the ensemble type to
#' use. Defaults to 1 (first ensemble type).
#' @param type Character string specifying the
#' variance-covariance estimator. One of \code{"HC1"}
#' (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A \eqn{p \times p}{p x p} variance-covariance matrix.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' vcov(plm_fit)
#' vcov(plm_fit, type = "HC3")
#' }
#'
#' @seealso \code{\link{hatvalues.ddml}},
#' \code{\link{confint.ddml}},
#' \code{\link{summary.ddml}}
#'
#' @family ddml inference
#' @method vcov ddml
#' @export
vcov.ddml <- function(object, ensemble_idx = 1,
                      type = "HC1", ...) {
  type <- validate_method_args(object,
    ensemble_idx = ensemble_idx, type = type)

  h <- if (type == "HC3") {
    stats::hatvalues(object, ensemble_idx = ensemble_idx)
  }#IF

  sc <- object$scores[, , ensemble_idx, drop = FALSE]
  dim(sc) <- dim(sc)[1:2]
  J_j <- object$J[, , ensemble_idx, drop = FALSE]
  dim(J_j) <- dim(J_j)[1:2]
  p <- ncol(sc)

  # Cluster aggregation: rowsum scores to cluster level
  clustered <- !is.null(object$cluster_variable) &&
    length(unique(object$cluster_variable)) < nrow(sc)
  if (clustered) {
    sc <- rowsum(sc, object$cluster_variable)
  }#IF

  n_eff <- nrow(sc)

  if (type == "HC3") {
    if (clustered) {
      h <- as.vector(tapply(h, object$cluster_variable, sum))
    }#IF
    sc <- sc / (1 - h)
  }#IF

  meat <- crossprod(sc) / n_eff
  J_inv <- csolve(J_j)
  V <- J_inv %*% meat %*% t(J_inv) / n_eff

  # HC1 degrees-of-freedom correction
  if (type == "HC1") {
    V <- V * n_eff / (n_eff - p)
  }#IF
  
  rownames(V) <- colnames(V) <- object$coef_names
  V
}#VCOV.DDML

#' Confidence Intervals for DDML Estimators
#'
#' @description Computes confidence intervals for one or more 
#' parameters in a fitted DDML model.
#'
#' @param object An object of class \code{ddml}.
#' @param parm A specification of which parameters are to be
#' given confidence intervals, either a vector of numbers
#' or a vector of names. If missing, all parameters are
#' considered.
#' @param level Confidence level. Default 0.95.
#' @inheritParams vcov.ddml
#' @param ... Currently unused.
#'
#' @return A matrix with columns for lower and upper bounds.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                 learners = list(what = ols),
#'                 sample_folds = 2, silent = TRUE)
#' confint(plm_fit)
#' confint(plm_fit, parm = "D1")
#' confint(plm_fit, level = 0.90)
#' }
#' 
#' @seealso \code{\link{vcov.ddml}}
#'
#' @family ddml inference
#' @method confint ddml
#' @importFrom stats vcov
#' @export
confint.ddml <- function(object, parm, level = 0.95,
                         ensemble_idx = 1,
                         type = "HC1", ...) {
  cf <- object$coefficients[, ensemble_idx]
  cf_names <- object$coef_names

  if (missing(parm)) {
    parm <- cf_names
  } else if (is.numeric(parm)) {
    parm <- cf_names[parm]
  } else {
    parm <- intersect(parm, cf_names)
    if (length(parm) == 0) {
      stop("None of the specified 'parm' were found in ",
           "the model coefficients.", call. = FALSE)
    }#IF
  }#IFELSE
  
  cf <- cf[parm]

  V <- vcov(object, ensemble_idx = ensemble_idx,
            type = type)
  se <- sqrt(diag(V))[parm]
  
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- c((1 - level) / 2, (1 + level) / 2) * 100
  colnames(ci) <- paste0(format(pct, digits = 3), " %")
  rownames(ci) <- parm
  ci
}#CONFINT.DDML

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
#' \code{\link{vcov.ddml}} for the HC0/HC1/HC3 formulas.
#'
#' @param object An object of class \code{ddml}.
#' @inheritParams vcov.ddml
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
#' @seealso \code{\link{vcov.ddml}}
#'
#' @family ddml inference
#' @method summary ddml
#' @export
summary.ddml <- function(object, type = "HC1", ...) {
  type <- match.arg(type, c("HC0", "HC1", "HC3"))
  
  single_learner <- is_single_learner(object$learners)
  ens_type <- if (single_learner) {
    "single base learner"
  } else {
    object$ensemble_type
  }#IFELSE
  
  nensb <- length(ens_type)
  p <- nrow(object$coefficients)

  inf <- array(0, dim = c(p, 4, nensb))

  for (j in seq_len(nensb)) {
    theta_j <- object$coefficients[, j]

    V <- stats::vcov(object, ensemble_idx = j, type = type)
    se <- sqrt(diag(V))
    z_val <- theta_j / se
    p_val <- 2 * stats::pnorm(abs(z_val),
                               lower.tail = FALSE)

    inf[, 1, j] <- theta_j
    inf[, 2, j] <- se
    inf[, 3, j] <- z_val
    inf[, 4, j] <- p_val
  }#FOR

  dimnames(inf) <- list(
    object$coef_names,
    c("Estimate", "Std. Error", "z value", "Pr(>|z|)"),
    ens_type
  )

  result <- list(
    coefficients = inf,
    type = type,
    model_type = class(object)[1],
    estimator_name = object$estimator_name,
    nobs = object$nobs,
    sample_folds = object$sample_folds,
    shortstack = object$shortstack,
    ensemble_type = ens_type)
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

  nensb <- dim(x$coefficients)[3]
  for (j in seq_len(nensb)) {
    if (nensb > 1) {
      cat("Ensemble type:",
          dimnames(x$coefficients)[[3]][j], "\n")
    }#IF
    tbl <- x$coefficients[, , j]
    if (!is.matrix(tbl)) {
      tbl <- matrix(tbl, nrow = 1,
                    dimnames = list(
                      dimnames(x$coefficients)[[1]],
                      dimnames(x$coefficients)[[2]]))
    }#IF
    stats::printCoefmat(tbl, digits = digits,
                        has.Pvalue = TRUE,
                        signif.stars = TRUE)
    if (j < nensb) cat("\n")
  }#FOR

  invisible(x)
}#PRINT.SUMMARY.DDML

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a ddml Object
#'
#' Extracts coefficient estimates, standard errors, test
#' statistics, and p-values from a \code{ddml} estimator in a
#' format compatible with \pkg{modelsummary} and the
#' \pkg{broom} ecosystem.
#'
#' @param x A \code{ddml} object.
#' @param ensemble_idx Integer index of the ensemble type to
#' report. Defaults to 1 (first ensemble type). Set to
#' \code{NULL} to return results for all ensemble types.
#' @param conf.int Logical. Include confidence interval
#' columns? Default \code{FALSE}.
#' @param conf.level Confidence level for intervals.
#' Default 0.95.
#' @param type Character string specifying the
#' variance-covariance estimator. One of \code{"HC1"}
#' (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#' \code{estimate}, \code{std.error}, \code{statistic},
#' \code{p.value}, and \code{ensemble_type}. If
#' \code{conf.int = TRUE}, also \code{conf.low} and
#' \code{conf.high}.
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
                      type = "HC1", ...) {
  type <- match.arg(type, c("HC1", "HC0", "HC3"))

  s <- summary(x, type = type)
  inf <- s$coefficients
  nensb <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(ensemble_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    if (any(ensemble_idx < 1) || any(ensemble_idx > nensb)) {
      stop(sprintf("ensemble_idx must be between 1 and %d", nensb),
           call. = FALSE)
    }#IF
    j_seq <- ensemble_idx
  }#IFELSE

  # Pre-allocate rows
  n_rows <- length(j_seq) * p
  term <- rep(dimnames(inf)[[1]], length(j_seq))
  ensemble_type <- rep(dimnames(inf)[[3]][j_seq], each = p)
  
  estimate <- numeric(n_rows)
  std.error <- numeric(n_rows)
  statistic <- numeric(n_rows)
  p.value <- numeric(n_rows)

  idx <- 1
  for (j in j_seq) {
    for (k in seq_len(p)) {
      estimate[idx] <- inf[k, 1, j]
      std.error[idx] <- inf[k, 2, j]
      statistic[idx] <- inf[k, 3, j]
      p.value[idx] <- inf[k, 4, j]
      idx <- idx + 1
    }#FOR
  }#FOR

  res <- data.frame(
    term = term,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    ensemble_type = ensemble_type,
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    z <- stats::qnorm((1 + conf.level) / 2)
    res$conf.low <- res$estimate - z * res$std.error
    res$conf.high <- res$estimate + z * res$std.error
  }#IF

  res
}#TIDY.DDML

#' Glance at a ddml Object
#'
#' Returns a one-row summary of model-level statistics,
#' compatible with \pkg{modelsummary} and the \pkg{broom}
#' ecosystem.
#'
#' @param x A \code{ddml} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with columns
#' \code{nobs}, \code{sample_folds}, \code{shortstack},
#' \code{ensemble_type}, and \code{model_type}.
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
#' @seealso \code{\link{tidy.ddml}},
#' \code{\link{summary.ddml}}
#'
#' @export
#' @method glance ddml
glance.ddml <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    sample_folds = x$sample_folds,
    shortstack = if (is.null(x$shortstack)) FALSE else x$shortstack,
    ensemble_type = paste(x$ensemble_type, collapse = ", "),
    model_type = class(x)[1],
    estimator_name = if (is.null(x$estimator_name)) class(x)[1] else x$estimator_name,
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML
