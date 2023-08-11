# Simple function to generate subsamples.
generate_subsamples <- function(nobs, sample_folds) {
  sampleframe <- rep(1:sample_folds, ceiling(nobs/sample_folds))
  sample_groups <- sample(sampleframe, size=nobs, replace=F)
  subsamples <- sapply(c(1:sample_folds),
                       function(x) {which(sample_groups == x)},
                       simplify = F)
  return(subsamples)
}#GENERATE_SUBSAMPLES

# Simple generalized inverse wrapper.
csolve <- function(X) {
  # Attempt inversion
  X_inv <- tryCatch(solve(X), error = function(e) NA)
  # If inversion failed, calculate generalized inverse
  if (any(is.na(X_inv))) {
    X_inv <- MASS::ginv(X)
  }#IF
  # Return (generalized) inverse
  return(X_inv)
}#CSOLVE

# Function to organize inference results for plm, pliv, and fpliv.
organize_inf_results <- function(fit_obj_list, ensemble_type, ...) {
  # Check whether a single or multiple ensemble_types were computed
  if (!methods::is(fit_obj_list, "list")) {
    # Compute standard errors, t-vals, and p-vals
    inf_results <- compute_inf_results_by_ensemble(fit_obj_list, ...)
    coef_names <- names(fit_obj_list$coefficients)
  } else {
    # Compute output for each ensemble_type
    ncoef <- length(fit_obj_list[[1]]$coefficients)
    nens <- length(ensemble_type)
    inf_results <- array(0, dim = c(ncoef, 4, nens))
    for (ens in 1:nens) {
      inf_results[, , ens] <-
        compute_inf_results_by_ensemble(fit_obj_list[[ens]], ...)
    }#FOR
    coef_names <- names(fit_obj_list[[1]]$coefficients)
  }#IFELSE
  # Name array dimensions
  dimnames(inf_results) <- list(coef_names,
                                c("Estimate", "Std. Error",
                                  "t value", "Pr(>|t|)"),
                                ensemble_type)
  # Print inference results
  return(inf_results)
}#COMPUTE_INF_RESULTS

# Function to compute std. errors, t-vals, and p-vals for plm, pliv, and fpliv.
compute_inf_results_by_ensemble <- function(fit_obj, ...) {
  # Data parameters
  ncoef <- length(fit_obj$coefficients)
  # Compute standard error, t-values, and p-vales
  Sigma <- sandwich::vcovHC(fit_obj, ...)
  std_errors <- sqrt(diag(Sigma))
  t_values <- fit_obj$coefficients / std_errors
  p_values <-  2 * sapply(abs(t_values), stats::pnorm, lower.tail = F)

  # Store results in a matrix
  inf_results <- array(0, dim = c(ncoef, 4, 1))
  inf_results[, 1, 1] <- fit_obj$coefficients
  inf_results[, 2, 1] <- std_errors
  inf_results[, 3, 1] <- t_values
  inf_results[, 4, 1] <- p_values
  # Return results
  return(inf_results)
}#COMPUTE_INF_RESULTS_BY_ENSEMBLE

# Function to organize inference results for ate and late.
organize_interactive_inf_results <- function(coef, psi_a, psi_b,
                                             ensemble_type) {
  # Data parameters
  nens <- length(ensemble_type)
  # Compute inference results by ensemble type
  inf_results <- matrix(0, nens, 4)
  for (ens in 1:nens) {
    inf_results[ens, ] <-
      compute_interactive_inf_results_by_ensemble(coef = coef[ens],
                                                  psi_a = psi_a[, ens],
                                                  psi_b = psi_b[, ens])
  }#FOR
  # Name and print output
  dimnames(inf_results) <- list(ensemble_type,
                                c("Estimate", "Std. Error",
                                  "t value", "Pr(>|t|)"))
  print(inf_results)
}#ORGANIZE_INTERACTIVE_INF_RESULTS

# Function to compute std. errors, t-vals, and p-vals for ate and late.
compute_interactive_inf_results_by_ensemble <- function(coef, psi_a, psi_b) {
  # Data parameters
  nobs <- length(psi_a)
  # Compute scores
  scores <- psi_a * coef + psi_b
  # Compute standard error, t-values, and p-vales
  std_error <- sqrt(mean(scores^2) / nobs) / abs(mean(psi_a))
  t_value <- coef / std_error
  p_value <-  2 * sapply(abs(t_value), stats::pnorm, lower.tail = F)
  # Store results in a matrix
  inf_results <- c(coef, std_error, t_value, p_value)
  # Return results
  return(inf_results)
}#COMPUTE_INTERACTIVE_INF_RESULTS_BY_ENSEMBLE
