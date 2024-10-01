# Function to organize inference results for plm, pliv, and fpliv.
organize_inf_results <- function(fit_obj_list, ensemble_type, cluster_variable,
                                 ...) {
  # Check whether a single or multiple ensemble_types were computed
  if (!methods::is(fit_obj_list, "list")) {
    # Compute standard errors, t-vals, and p-vals
    inf_results <- compute_inf_results_by_ensemble(fit_obj_list,
                                                   cluster_variable, ...)
    coef_names <- names(fit_obj_list$coefficients)
  } else {
    # Compute output for each ensemble_type
    ncoef <- length(fit_obj_list[[1]]$coefficients)
    nens <- length(ensemble_type)
    inf_results <- array(0, dim = c(ncoef, 4, nens))
    for (ens in 1:nens) {
      inf_results[, , ens] <-
        compute_inf_results_by_ensemble(fit_obj_list[[ens]], cluster_variable,
                                        ...)
    }#FOR
    coef_names <- names(fit_obj_list[[1]]$coefficients)
  }#IFELSE
  # Name array dimensions
  dimnames(inf_results) <- list(coef_names,
                                c("Estimate", "Std. Error",
                                  "t value", "Pr(>|t|)"),
                                ensemble_type)
  # Return inference results
  inf_results
}#COMPUTE_INF_RESULTS

# Function to compute std. errors, t-vals, and p-vals for plm, pliv, and fpliv.
compute_inf_results_by_ensemble <- function(fit_obj, cluster_variable, ...) {
  # Data parameters
  ncoef <- length(fit_obj$coefficients)
  cluster <- !identical(seq_along(fit_obj$residuals), cluster_variable)
  # Compute standard error, t-values, and p-vales
  if (cluster) {
    Sigma <- sandwich::vcovCL(fit_obj, cluster = cluster_variable, ...)
  } else {
    Sigma <- sandwich::vcovHC(fit_obj, ...)
  }#IFELSE
  std_errors <- sqrt(diag(Sigma))
  t_values <- fit_obj$coefficients / std_errors
  p_values <-  2 * vapply(abs(t_values), stats::pnorm,
                          FUN.VALUE = 1, lower.tail = F)

  # Store results in a matrix
  inf_results <- array(0, dim = c(ncoef, 4, 1))
  inf_results[, 1, 1] <- fit_obj$coefficients
  inf_results[, 2, 1] <- std_errors
  inf_results[, 3, 1] <- t_values
  inf_results[, 4, 1] <- p_values
  # Return results
  inf_results
}#COMPUTE_INF_RESULTS_BY_ENSEMBLE

# Function to organize inference results for ate and late.
organize_interactive_inf_results <- function(coef, psi_a, psi_b,
                                             ensemble_type, cluster_variable) {
  # Data parameters
  nens <- length(ensemble_type)
  # Compute inference results by ensemble type
  inf_results <- matrix(0, nens, 4)
  for (ens in 1:nens) {
    inf_results[ens, ] <-
      compute_interactive_inf_results_by_ensemble(coef = coef[ens],
                                                  psi_a = psi_a[, ens],
                                                  psi_b = psi_b[, ens],
                                                  cluster_variable)
  }#FOR
  # Name and return output
  dimnames(inf_results) <- list(ensemble_type,
                                c("Estimate", "Std. Error",
                                  "t value", "Pr(>|t|)"))
  inf_results
}#ORGANIZE_INTERACTIVE_INF_RESULTS

# Function to compute std. errors, t-vals, and p-vals for ate and late.
compute_interactive_inf_results_by_ensemble <- function(coef, psi_a, psi_b,
                                                        cluster_variable) {
  # Data parameters
  nobs <- length(psi_a)
  cluster <- !identical(1:nobs, cluster_variable)
  # Compute scores
  scores <- psi_a * coef + psi_b
  # Compute standard error, t-values, and p-vales
  if (cluster) {
    # Aggregate scores to cluster level for clustered standard errors
    scores <- stats::aggregate(scores, by = list(cluster_variable),
                               FUN = sum)[, 2]
    psi_a <- stats::aggregate(psi_a, by = list(cluster_variable),
                              FUN = sum)[, 2]
    nobs <- length(scores)
  }#IF
  std_error <- sqrt(mean(scores^2) / nobs) / abs(mean(psi_a))
  t_values <- coef / std_error
  p_value <-  2 * vapply(abs(t_values), stats::pnorm,
                         FUN.VALUE = 1, lower.tail = F)
  # Store results in a matrix
  inf_results <- c(coef, std_error, t_values, p_value)
  # Return results
  inf_results
}#COMPUTE_INTERACTIVE_INF_RESULTS_BY_ENSEMBLE
