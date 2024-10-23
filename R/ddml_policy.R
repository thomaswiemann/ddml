#' Estimators of a Multi-Action Policy's Expected Value.
#'
#' @family ddml
#'
#' @seealso [ddml::summary.ddml_policy()]
#'
#' @description Estimators of the expected value of a multi-action policy.
#'
#' @details \code{ddml_policy} provides a double/debiased machine learning
#'     estimator for the expected value of a multi-action policy
#'     \eqn{\pi:\operatorname{supp} X \to \operatorname{supp} D = \{d_1, \ldots, d_K\}}
#'     given by
#'
#' \eqn{E[Y(\pi(X))],}
#'
#' where \eqn{Y(\cdot)} are the potential payoffs at different values of the
#'     policy.
#'
#' @inheritParams ddml_ate
#' @param D The observed discrete (potentially multi-valued) treatment variable.
#' @param policy The policy-assigned treatment variable.
#' @param policy_levels A vector of the unique values the policy. These are the
#'     treatment levels that the reduced form functions will be estimated for.
#' @param omega A vector of sampling weights \eqn{\omega_i} used for
#'     aggregating scores.
#' @param subsamples_byD List of lists corresponding to the unique treatment
#'     levels. Each list contains vectors with sample indices for
#'     cross-fitting.
#' @param cv_subsamples_byD List of lists, each corresponding to one of the
#'     unique treatment levels. Each of the lists contains lists, each
#'     corresponding to a subsample and contains vectors with subsample indices
#'     for cross-validation.
#' @param oos_pred A list containing arrays for the reduced form estimates. Can
#'     Can be passed from the output of previously-fitted \code{ddml_policy}
#'     objects. This allows rapid evaluation of multiple policies without
#'     re-estimation of the nuisance functions.
#'
#' @return \code{ddml_policy} returns an object of S3 class \code{ddml_policy}.
#'     An object of class \code{ddml_policy} is a list containing the following
#'     components:
#'     \describe{
#'         \item{\code{policy_value}}{A vector with the expected policy value
#'             estimates.}
#'         \item{\code{weights}}{A list of matrices, providing the weight
#'             assigned to each base learner (in chronological order) by the
#'             ensemble procedure.}
#'         \item{\code{mspe}}{A list of matrices, providing the MSPE of each
#'             base learner (in chronological order) computed by the
#'             cross-validation step in the ensemble construction.}
#'         \item{\code{psi_a}, \code{psi_b}}{Matrices needed for the computation
#'             of scores. Used in [ddml::summary.ddml_policy()].}
#'         \item{\code{oos_pred}}{List of arrays, providing the reduced form
#'             predicted values.}
#'         \item{\code{learners},\code{learners_DX},\code{cluster_variable},
#'             \code{subsamples},\code{subsamples_byD},
#'             \code{cv_subsamples},\code{cv_subsamples_byD},
#'             \code{ensemble_type},\code{omega}}{Pass-through of
#'             selected user-provided arguments. See above.}
#'     }
#' @export
#'
#' @references
#' TO DO!
#'
#' @examples
#' # TO DO!
ddml_policy <- function(y, D, X,
                        policy,
                        policy_levels = sort(unique(policy)),
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
                        balance_on = "clusters",
                        subsamples = NULL,
                        subsamples_byD = NULL,
                        cv_subsamples = NULL,
                        cv_subsamples_byD = NULL,
                        oos_pred = NULL,
                        trim = 0.01,
                        omega = rep(1, length(y)),
                        silent = FALSE) {

  # Data parameters
  nobs <- length(y)
  n_policy_levels <- length(policy_levels)
  is_D <- rep(list(NULL), n_policy_levels)
  for (d in 1:n_policy_levels) {
    is_D[[d]] <- which(D == policy_levels[d])
  }#FOR

  # Was oos_pred passed. If so, jump straight to policy evaluation.
  if (is.null(oos_pred)) {

    # Check whether ddml uses conventional stacking w/ data driven weights
    w_cv <- !shortstack &
      any(ensemble_type %in% c("nnls", "nnls1", "singlebest", "ols")) &
      (class(learners[[1]]) != "function" |
         class(learners_DX[[1]]) != "function")

    # Create crossfitting and cv tuples
    indxs <- get_all_indx(cluster_variable = cluster_variable,
                          sample_folds = sample_folds, cv_folds = cv_folds,
                          D = D, by_D = TRUE, stratify = stratify,
                          balance_on = balance_on,
                          subsamples = subsamples,
                          subsamples_byD = subsamples_byD,
                          cv_subsamples = cv_subsamples,
                          cv_subsamples_byD = cv_subsamples_byD,
                          compute_cv_indices = w_cv,
                          compute_aux_X_indices = TRUE)
    subsamples = indxs$subsamples
    subsamples_byD = indxs$subsamples_byD
    cv_subsamples  = indxs$cv_subsamples
    cv_subsamples_byD = indxs$cv_subsamples_byD

    # Print to progress to console
    if (!silent) cat("DDML estimation in progress. \n")

    # Compute est. of E[y|D=d,X] for all policy values d
    y_Dd_X_res_list <- rep(list(NULL), n_policy_levels)
    for (d in 1:n_policy_levels) {
      y_Dd_X_res_list[[d]]$fit <- get_CEF(y[is_D[[d]]],
                                          X[is_D[[d]], , drop = F],
                                          learners = learners,
                                          ensemble_type = ensemble_type,
                                          shortstack = shortstack,
                                          custom_ensemble_weights =
                                            custom_ensemble_weights,
                                          subsamples =
                                            indxs$subsamples_byD[[d]],
                                          cv_subsamples =
                                            indxs$cv_subsamples_byD[[d]],
                                          auxiliary_X =
                                            get_auxiliary_X(indxs$aux_indx[[d]],
                                                            X),
                                          silent = silent,
                                          progress = paste0("E[Y|D=",
                                                            policy_levels[d],
                                                            ",X]: "))
      y_Dd_X_res_list[[d]]$d <- policy_levels[d]
    }#FOR

    # Compute est. of Pr(D=d|X) for all policy values d
    D_X_res_list <- rep(list(NULL), n_policy_levels - 1)
    for (d in 1:(n_policy_levels - 1)) {
      # Pr(D=d|X)
      D_X_res_list[[d]] <- get_CEF(1 * (D == policy_levels[d]), X,
                                          learners = learners_DX,
                                          ensemble_type = ensemble_type,
                                          shortstack = shortstack,
                                          custom_ensemble_weights =
                                            custom_ensemble_weights_DX,
                                          subsamples = indxs$subsamples,
                                          cv_subsamples =
                                            indxs$cv_subsamples,
                                          silent = silent,
                                          progress = paste0("Pr[D=",
                                                            policy_levels[d],
                                                            ",X]: "))
    }#FOR

    # Update ensemble type to account for (optional) custom weights
    ensemble_type <- dimnames(y_Dd_X_res_list[[1]][[1]]$weights)[[2]]
    nensb <- ifelse(is.null(ensemble_type), 1, length(ensemble_type))

    # Check whether multiple ensembles are computed simultaneously
    multiple_ensembles <- nensb > 1

    # Construct reduced form variables
    g_X_Dd <- extrapolate_CEF(D = D,
                              CEF_res_byD = y_Dd_X_res_list,
                              aux_indx = indxs$aux_indx)
    m_X <- array(0, dim = c(nobs, nensb, n_policy_levels))
    for (d in 1:(n_policy_levels-1)) {
      m_X[ , , d] <- D_X_res_list[[d]]$oos_fitted
    }#FOR
    m_X[, , n_policy_levels] <-  1 - apply(m_X, c(1, 2), sum)
  } else {
    # If oos_pred was passed, simply collect the reduced form variables
    g_X_Dd <- oos_pred$EY_Dd_X
    m_X <- oos_pred$EDd_X
    nensb <- dim(m_X)[2]
  }#IFELSE

  # Trim propensity scores, return warnings
  m_X_tr <- m_X
  for (d in 1:n_policy_levels) {
    m_X_tr[ , , d] <-
      trim_propensity_scores(as.matrix(m_X[, , d, drop = FALSE]),
                             trim, ensemble_type)
  }#FOR

  # Compute the score using the constructed reduced form variables
  y_copy <- matrix(rep(y, nensb), nobs, nensb)
  psi_b <- matrix(rep(0, nobs * nensb), nobs, nensb)
  for (d in 1:n_policy_levels) {
    policy_copy <-
      matrix(rep(policy == policy_levels[d], nensb), nobs, nensb) * 1
    Dd_copy <- matrix(rep(D == policy_levels[d], nensb), nobs, nensb) * 1
    psi_b <- psi_b + policy_copy *
      (Dd_copy * (y_copy - g_X_Dd[, , d]) / m_X_tr[, , d] + g_X_Dd[, , d])
  }#FOR
  psi_a <- matrix(-1, nobs, nensb)

  # Scale the scores by omega (optional)
  if (!identical(omega, rep(1, nobs))) {
    psi_b <- sweep(psi_b, 1, omega, "*")
    psi_a <- sweep(psi_a, 1, omega, "*")
  }#IF

  # Policy value
  policy_value <- -colMeans(psi_b) / colMeans(psi_a)
  names(policy_value) <- ensemble_type

  # Organize complementary ensemble output
  weights <- mspe <- list(NULL)
  if (exists("y_Dd_X_res_list")) {
    for (d in 1:n_policy_levels) {
      weights[[d]] <- y_Dd_X_res_list[[d]][[1]]$weights
      mspe[[d]] <- y_Dd_X_res_list[[d]][[1]]$mspe
    }#FOR

    for (d in 1:(n_policy_levels-1)) {
      weights[[n_policy_levels + d]] <- D_X_res_list[[d]]$weights
      mspe[[n_policy_levels + d]] <- D_X_res_list[[d]]$mspe
    }#FOR
    names(weights) <- names(mspe) <-
      c(sapply(policy_levels, function (d) paste0("y_X_p", d)),
        sapply(policy_levels[-n_policy_levels],
               function (d) paste0("D", d, "_X")))
  }#IF

  # Organize reduced form predicted values
  oos_pred <- list(EY_Dd_X = g_X_Dd, EDd_X = m_X)

  # Organize output
  ddml_fit <- list(policy_value = policy_value,
                   weights = weights, mspe = mspe,
                   psi_a = psi_a, psi_b = psi_b,
                   oos_pred = oos_pred,
                   learners = learners,
                   learners_DX = learners_DX,
                   cluster_variable = cluster_variable,
                   subsamples = subsamples,
                   subsamples_byD = subsamples_byD,
                   cv_subsamples  = cv_subsamples,
                   cv_subsamples_byD = cv_subsamples_byD,
                   ensemble_type = ensemble_type, omega = omega)

  # Print estimation progress
  if (!silent) cat("DDML estimation completed. \n")

  # Amend class and return
  class(ddml_fit) <- "ddml_policy"
  return(ddml_fit)
}#DDML_Policy

#' Inference Method for the Multi-Action Policy Value Estimator.
#'
#' @description Inference method for multi-action policy value estimator. By
#'     default, standard errors are heteroskedasiticty-robust. If the \code{ddml}
#'     estimator was computed using a \code{cluster_variable}, the standard
#'     errors are also cluster-robust by default.
#'
#' @param object An object of class \code{ddml_policy} as fitted by
#'     [ddml::ddml_policy()].
#' @param ... Currently unused.
#'
#' @return A matrix with inference results.
#'
#' @export
#'
#' @examples
#' # TO DO!
summary.ddml_policy <- function(object, ...) {
  # Check whether stacking was used, replace ensemble type if TRUE
  single_learner <- ("what" %in% names(object$learners))
  if (single_learner) object$ensemble_type <- " "
  # Compute and return inference results
  coefficients <- organize_interactive_inf_results(coef = object$policy_value,
                                                   psi_a = object$psi_a,
                                                   psi_b = object$psi_b,
                                                   ensemble_type =
                                                     object$ensemble_type,
                                                   cluster_variable =
                                                     object$cluster_variable)
  class(coefficients) <- c("summary.ddml_policy", class(coefficients))
  coefficients
}#SUMMARY.DDML_POLICY

#' Print Method for the Multi-Action Policy Value Estimator.
#'
#' @description Print method for the multi-action policy value estimator.
#'
#' @param x An object of class \code{summary.ddml_policy} as returned by
#'     [ddml::summary.ddml_policy()].
#' @param digits The number of significant digits used for printing.
#' @param ... Currently unused.
#'
#' @return NULL.
#'
#' @export
#'
#' @examples
#' # TO DO!
print.summary.ddml_policy <- function(x, digits = 3, ...) {
  cat("Policy value estimation results: \n \n")
  class(x) <- class(x)[-1]
  print(x, digits = digits)
}#PRINT.SUMMARY.DDML_POLICY
