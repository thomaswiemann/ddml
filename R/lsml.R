#' Compute (weighted) least squares estimator.
#'
#' Compute (weighted) least squares estimator.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param const A boolean indicating inclusion of a constant.
#' @param w An optional weights vector.
#'
#' @return \code{ols} returns an object of S3 class "\code{ols}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class. The function \code{summary} computes the corresponding
#'     standard errors, t-statistics, and p-values.
#'
#' An object of class "\code{ols}" is a list containig the following components:
#' \describe{
#' \item{\code{coef}}{A vector of least squares coefficients.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X}}{A feature matrix.}
#' \item{\code{const}}{A boolean indicating inclusion of a constant.}
#' \item{\code{w}}{An optional weights vector.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(100*3), 100, 3) # Simulate features
#' y <- 1 + X %*% c(-1, 1, 0) + rnorm(100) # Simulate linear model
#' ols(y, X, const = T) # Compute least squares fit
#'
#' @export lsml
lsml <- function(response, features,
                   model, args,
                   f_split = replicate(2, c(1:ncol(features)), F),
                   max_iter = 10, eps = 1e-6,
                   rnd_seed = sample(c(1:999999), 1),
                   silent = TRUE){
  # inital OLS fit with all features
  fit.LS <- OLS(response, features[,f_split[[1]]])
  resid.LS <- response - predict(fit.LS)

  # loop
  for(k in 1:max_iter){
    # Repeat same rnd seed
    if(!is.null(rnd_seed)) set.seed(rnd_seed)

    # Fit ML model on LS residual
    args$response <- resid.LS; args$features <- features[,f_split[[2]]]
    fit.ML <- do.call(what = model, args = args)
    fitted.ML <- predict(fit.ML, features[,f_split[[2]]])

    # Fit LS with ML fitted values
    fit.LS <- OLS(response, cbind(features[,f_split[[1]]], fitted.ML))
    resid.LS <- response  - predict(fit.LS, cbind(features[,f_split[[1]]], mean(fitted.ML))) # OLS fitted without h(X)

    # Check convergence condition
    MSPE.k <- mean((response  - predict(fit.LS))^2)
    if(k > 1 && abs(MSPE.k-MSPE.nk) < eps) break
    MSPE.nk <- MSPE.k
    print(MSPE.k)

    if(!silent) print(paste0("Progress: ", k/max_iter))
  }#FOR
  if(k == max_iter) warning('Coefficients did not converge.')

  # organize output
  fit.myLSML <- list(fit.LS = fit.LS, fit.ML = fit.ML,
                     response = response, features = features,
                     f_split = f_split)
  class(fit.myLSML) <- 'MYLSML'
  return(fit.myLSML)
}#LSML

# Complementary methods ========================================================
#' Predict method for lsml fits.
#'
#' Predict method for lsml fits.
#'
#' @export predict.ols
predict.MYLSML <- function(fit.myLSML,
                           new_features = NULL){
  # Check for new data
  if(is.null(new_features)) new_features <- fit.myLSML$features

  # Pedict ML model
  fitted.ML <- predict(fit.myLSML$fit.ML, new_features[,fit.myLSML$f_split[[2]]])
  # Predict LS + ML model
  fitted <- predict(fit.myLSML$fit.LS, cbind(new_features[,fit.myLSML$f_split[[1]]],fitted.ML))

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
