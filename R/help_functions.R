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
  for (j in length(nensb)) {
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
