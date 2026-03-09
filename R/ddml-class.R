#' Common Output Structure for DDML Estimators
#'
#' @name ddml-class
#'
#' @description All \code{ddml_*} estimators (\code{\link{ddml_plm}},
#'     \code{\link{ddml_pliv}}, \code{\link{ddml_fpliv}},
#'     \code{\link{ddml_ate}}, \code{\link{ddml_att}},
#'     \code{\link{ddml_late}}, \code{\link{ddml_apo}}) return
#'     objects that inherit from S3 class \code{"ddml"}.
#'
#' Each object is a list containing the components described below.
#'     Estimator-specific fields (e.g., pass-through learner
#'     arguments) are documented on the individual estimator pages.
#'
#' @details All \code{ddml_*} estimators target a low-dimensional
#'     parameter \eqn{\theta_0} identified by a moment condition
#'
#' \deqn{E[m(W; \theta_0, \eta_0)] = 0,}
#'
#'     where \eqn{W} denotes observed random variables and
#'     \eqn{\eta_0} is a (potentially high-dimensional) nuisance
#'     parameter, typically a vector of conditional expectation
#'     functions. The score \eqn{m} is \emph{Neyman orthogonal},
#'     meaning the moment condition is locally insensitive to
#'     perturbations of \eqn{\eta} around \eqn{\eta_0}.
#'
#'     Each estimator in \code{ddml} uses a score that decomposes
#'     linearly in \eqn{\theta}:
#'
#' \deqn{m(W_i; \theta, \eta) = \psi_b(W_i; \eta) +
#'     \psi_a(W_i; \eta)\,\theta.}
#'
#'     See the individual estimator pages for the specific forms
#'     of \eqn{\psi_a} and \eqn{\psi_b}.
#'
#'     Estimation proceeds via cross-fitting: the sample is randomly
#'     partitioned into \eqn{K} folds \eqn{\{I_k\}_{k=1}^K}. For
#'     each fold \eqn{k}, nuisance parameters are estimated on the
#'     complementary folds (\eqn{\hat\eta_{-k}}) and the scores are
#'     evaluated on fold \eqn{k}. The DML estimator
#'     \eqn{\hat\theta} solves
#'
#' \deqn{\frac{1}{n} \sum_{k=1}^{K} \sum_{i \in I_k}
#'     m(W_i; \hat\theta, \hat\eta_{-k}) = 0.}
#'
#'     For the linear scores used in \code{ddml}, this yields the
#'     closed-form solution
#'
#' \deqn{\hat\theta = -\hat{J}^{-1}\,\hat\psi_b,}
#'
#'     where \eqn{\hat\psi_b = n^{-1}\sum_{k}\sum_{i \in I_k}
#'     \psi_b(W_i; \hat\eta_{-k})} and
#'     \eqn{\hat{J} = n^{-1}\sum_k \sum_{i \in I_k}
#'     \psi_a(W_i; \hat\eta_{-k})} is the sample Jacobian.
#'
#'     Inference is based on the sandwich variance estimator
#'
#' \deqn{\hat\Sigma = \hat{J}^{-1} \left(\frac{1}{n}
#'     \sum_i m_i m_i^\top\right)
#'     \hat{J}^{-\top} / n}
#'
#'     where \eqn{m_i = m(W_i; \hat\theta, \hat\eta_{-k(i)})};
#'     see \code{\link{vcov.ddml}} for HC0/HC1/HC3 variants.
#'     Under regularity conditions and sufficient convergence of \eqn{\hat\eta}, 
#'     the DML estimator is asymptotically normal:
#'
#' \deqn{\sqrt{n}\,\hat\Sigma^{-1/2}(\hat\theta - \theta_0)
#'     \overset{d}{\to} N(0, I).}
#'
#'     Further details and regularity conditions are given in
#'     Chernozhukov et al. (2018).
#'
#' @references
#' Ahrens A, Chernozhukov V, Hansen C B, Kozbur D, Schaffer M E,
#'     Wiemann T (2026). "An Introduction to Double/Debiased Machine
#'     Learning." Journal of Economic Literature, forthcoming.
#'
#' Chernozhukov V, Chetverikov D, Demirer M, Duflo E, Hansen C B,
#'     Newey W, Robins J (2018). "Double/debiased machine learning
#'     for treatment and structural parameters." The Econometrics
#'     Journal, 21(1), C1-C68.
#'
#' @section Common output components:
#' \describe{
#'     \item{\code{coefficients}}{A matrix of estimated target
#'         parameters: rows correspond to components of
#'         \eqn{\theta}, columns to ensemble types.}
#'     \item{\code{ensemble_weights}}{A named list. Each element
#'         is a weight matrix (or 3D array when
#'         \code{shortstack = TRUE}) showing the weight assigned
#'         to each base learner by the ensemble procedure for
#'         the corresponding nuisance equation.}
#'     \item{\code{mspe}}{A named list of per-learner MSPEs from
#'         the cross-validation step in ensemble construction.}
#'     \item{\code{r2}}{A named list of per-learner out-of-sample
#'         R-squared values.}
#'     \item{\code{psi_a}, \code{psi_b}}{Score component lists
#'         (length \code{nensb}). \code{psi_a[[j]]} is an
#'         \code{(n x p x p)} array and \code{psi_b[[j]]} is an
#'         \code{(n x p)} matrix. Used internally by
#'         \code{\link{hatvalues.ddml}}.}
#'     \item{\code{scores}}{A 3D array of evaluated Neyman
#'         orthogonal scores (\code{n x p x nensb}).}
#'     \item{\code{J}}{A 3D array of evaluated Jacobians
#'         (\code{p x p x nensb}).}
#'     \item{\code{fitted}}{A named list of per-equation
#'         cross-fitted prediction objects. Can be passed back
#'         via the \code{fitted} argument together with
#'         \code{splits} to skip cross-fitting on
#'         re-estimation.}
#'     \item{\code{splits}}{The data splitting structure
#'         (subsamples, CV subsamples, and any stratification
#'         indices).}
#'     \item{\code{ensemble_type}}{Character vector of ensemble
#'         types used.}
#'     \item{\code{cluster_variable}}{The cluster variable
#'         vector used for sample splitting and inference.}
#'     \item{\code{nobs}}{Number of observations.}
#'     \item{\code{sample_folds}}{Number of cross-fitting folds.}
#'     \item{\code{shortstack}}{Logical indicating whether
#'         short-stacking was used.}
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{coef_names}}{Character vector of coefficient
#'         names.}
#'     \item{\code{estimator_name}}{Character string identifying
#'         the estimator (e.g., \code{"Partially Linear Model"}).}
#' }
#'
#' @section S3 methods:
#' The following generic methods are available for all
#'     \code{ddml} objects: \code{\link{summary.ddml}},
#'     \code{\link{coef.ddml}}, \code{\link{vcov.ddml}},
#'     \code{\link{confint.ddml}}, \code{\link{hatvalues.ddml}},
#'     \code{\link{nobs.ddml}}, \code{\link{tidy.ddml}},
#'     \code{\link{glance.ddml}}, and
#'     \code{\link{diagnostics}}.
#'
#' @param y The outcome variable.
#' @param D A matrix of endogenous variables.
#' @param X A (sparse) matrix of control variables.
#' @param learners May take one of two forms, depending on whether a
#'     single learner or stacking with multiple learners is used for
#'     estimation of the conditional expectation functions.
#'     If a single learner is used, \code{learners} is a list with
#'     two named elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function
#'             must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to
#'             \code{what}.}
#'     }
#'     If stacking with multiple learners is used, \code{learners} is
#'     a list of lists, each containing three named elements:
#'     \itemize{
#'         \item{\code{what} The base learner function. The function
#'             must be such that it predicts a named input \code{y}
#'             using a named input \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to
#'             \code{what}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to control variables in \code{X} that
#'             are passed to the base learner.}
#'     }
#'     Omission of the \code{args} element results in default
#'     arguments being used in \code{what}. Omission of
#'     \code{assign_X} results in inclusion of all variables in
#'     \code{X}.
#' @param sample_folds Number of cross-fitting folds.
#' @param ensemble_type Ensemble method to combine base learners into
#'     final estimate of the conditional expectation functions.
#'     Possible values are:
#'     \itemize{
#'         \item{\code{"nnls"} Non-negative least squares.}
#'         \item{\code{"nnls1"} Non-negative least squares with the
#'             constraint that all weights sum to one.}
#'         \item{\code{"singlebest"} Select base learner with minimum
#'             MSPE.}
#'         \item{\code{"ols"} Ordinary least squares.}
#'         \item{\code{"average"} Simple average over base learners.}
#'     }
#'     Multiple ensemble types may be passed as a vector of strings.
#' @param shortstack Boolean to use short-stacking.
#' @param cv_folds Number of folds used for cross-validation in
#'     ensemble construction.
#' @param custom_ensemble_weights A numerical matrix with
#'     user-specified ensemble weights. Each column corresponds to a
#'     custom ensemble specification, each row corresponds to a base
#'     learner in \code{learners} (in chronological order). Optional
#'     column names are used to name the estimation results
#'     corresponding the custom ensemble specification.
#' @param cluster_variable A vector of cluster indices.
#' @param silent Boolean to silence estimation updates.
#' @param parallel An optional named list with parallel processing
#'     options. When \code{NULL} (the default), computation is
#'     sequential. Supported fields:
#'     \describe{
#'         \item{\code{cores}}{Number of cores to use.}
#'         \item{\code{export}}{Character vector of object names to
#'             export to parallel workers (for custom learners that
#'             reference global objects).}
#'         \item{\code{packages}}{Character vector of additional
#'             package names to load on workers (for custom learners
#'             that use packages not imported by \code{ddml}).}
#'     }
#' @param fitted An optional named list of per-equation cross-fitted
#'     predictions, typically obtained from a previous fit via
#'     \code{fit$fitted}. When supplied (together with \code{splits}),
#'     base learners are not re-fitted; only ensemble weights are
#'     recomputed. This allows fast re-estimation with a different
#'     \code{ensemble_type}. See \code{\link{ddml_plm}} for
#'     an example.
#' @param splits An optional list of sample split objects, typically
#'     obtained from a previous fit via \code{fit$splits}. Must be
#'     supplied when \code{fitted} is provided. Can also be used
#'     standalone to provide pre-computed sample folds.
#' @param save_crossval Logical indicating whether to store the inner
#'     cross-validation residuals used for ensemble weight
#'     computation. Default \code{TRUE}. When \code{TRUE}, subsequent
#'     pass-through calls with data-driven ensembles (e.g.,
#'     \code{"nnls"}) reproduce per-fold weights exactly. Set to
#'     \code{FALSE} to reduce object size at the cost of approximate
#'     weight recomputation.
#' @param ... Deprecated arguments are still accepted for backward
#'     compatibility but should be replaced with the \code{splits}
#'     argument. See individual estimator pages for details.
#'
#' @family ddml estimators
NULL
