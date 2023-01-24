#' Title
#'
#' @param nobs abc
#' @param sample_folds abc
#'
#' @return object
#' @export
#'
#' @examples
#' generate_subsamples(1000, 3)
generate_subsamples <- function(nobs, sample_folds) {
  sampleframe <- rep(1:sample_folds, ceiling(nobs/sample_folds))
  sample_groups <- sample(sampleframe, size=nobs, replace=F)
  subsamples <- sapply(c(1:sample_folds),
                       function(x) {which(sample_groups == x)},
                       simplify = F)
  return(subsamples)
}#GENERATE_SUBSAMPLES

# Utility to print progress ot console
update_progress <- function(silent) {
  if (!silent) cat(" -- Done! \n")
}#UPDATE_PROGRESS


# Simple generalized inverse wrapper.
csolve <- function(X) {
  # Attmept inversion
  X_inv <- tryCatch(solve(X), error = function(e) NA)
  # If inversion failed, calculate generalized inverse
  if (any(is.na(X_inv))) {
    X_inv <- MASS::ginv(X)
    #warning("Inversion failed. The generalized inverse was calulated instead.")
  }#IF
  #Return (generalized) inverse
  return(X_inv)
}#CSOLVE
