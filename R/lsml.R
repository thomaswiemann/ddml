#' Compute an LS + ML estimator.
#'
#' Compute an LS + ML estimator.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param model A constructor function for a nonlinear model.
#' @param args Arguments to be passed to \code{model}.
#' @param split_X A list of two vectors of indices. The first contains the
#'     indices of variables in \code{X} that are used in the linear model. The
#'     second contains the indices of variables in \code{X} that are used in the
#'     nonlinear model.
#' @param max_iter An upperbound on the number of iterations until convergence.
#' @param eps The convergence criterion.
#' @param rnd_seed A random seed number to fix inherit randomness in some ML
#'     procedures (e.g., varibale selection random forests) across iterations.
#' @param silent A boolean indicating whether the iteration progress should be
#'     printed to the console.
#'
#' @return \code{lsml} returns an object of S3 class "\code{lsml}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class.
#'
#' An object of class "\code{lsml}" is a list containig the following components:
#' \describe{
#' \item{\code{ls_fit}}{An object of class \code{ols}.}
#' \item{\code{ml_fit}}{An object corresponding to the constructor function
#'     passed via \code{model}.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X}}{A feature matrix.}
#' \item{\code{split_X}}{Passthrough of the input under the same name.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export lsml
lsml <- function(y, X,
                 model, args = list(),
                 split_X = replicate(2, c(1:ncol(X)), F),
                 max_iter = 10, eps = 1e-6,
                 rnd_seed = sample(c(1:999999), 1),
                 silent = TRUE){
  # Inital OLS fit with all features
  ls_fit  <- ols(y, X[, split_X[[1]]])
  resid_ols <- y - predict(ls_fit)
  # Loop until convergence
  args$X <- X[, split_X[[2]]]
  for (k in 1:max_iter) {
    # Repeat same rnd seed
    if(!is.null(rnd_seed)) set.seed(rnd_seed)
    # Fit ML model on ols residual
    args$y <- resid_ols
    ml_fit <- do.call(what = model, args = args)
    fitted_ml <- predict(ml_fit, X[, split_X[[2]]])
    # Fit LS with ML fitted values
    ls_fit <- ols(y, cbind(X[, split_X[[1]]], fitted_ml))
    resid_ols <- y  - predict(ls_fit, cbind(X[, split_X[[1]]],
                                            mean(fitted_ml)))
    # Check convergence condition
    mspe_k <- mean((y  - predict(ls_fit))^2)
    if(k > 1 && abs(mspe_k-mspe_nk) < eps) break
    mspe_nk <- mspe_k
    if(!silent) print(paste0("Progress: ", k / max_iter))
  }#FOR
  # Check for convergence
  if(k == max_iter) warning('Coefficients did not converge.')
  # organize output
  output <- list(ls_fit = ls_fit, ml_fit = ml_fit,
                 y = y, X = X,
                 split_X = split_X)
  class(output) <- 'lsml'
  return(output)
}#LSML

# Complementary methods ========================================================
#' Predict method for lsml fits.
#'
#' Predict method for lsml fits.
#'
#' @export predict.lsml
predict.lsml <- function(obj, newdata = NULL){
  # Check for new data
  if(is.null(newdata)) newdata <- obj$X
  # Pedict ML model
  fitted_ml <- predict(obj$ml_fit, newdata[, obj$split_X[[2]]])
  # Predict LS + ML model
  fitted <- predict(obj$ls_fit, cbind(newdata[, obj$split_X[[1]]], fitted_ml))
  # Return fitted values
  return(fitted)
}#PREDICT.LSML

#' Instrument selection for lsml fits.
#'
#' Instrument selection for lsml fits. Always returns \code{TRUE}.
#'
#' @export any_iv.lsml
any_iv.lsml <- function(obj, ...){
  TRUE
}#ANY_IV.LSML
