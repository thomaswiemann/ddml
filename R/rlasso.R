#' Compute Lasso with plug-in penalty level and loadings.
#'
#' Compute Lasso with plug-in penalty level and loadings.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param include An optional vector of feature indices with zero penalty
#'     loadings.
#' @param c A number strictly greater than 1 for the computation of the penalty
#'     level.
#' @param gamma A positive number for the computation of the penalty level.
#' @param post A boolean indicating whether a post-Lasso ols model should be
#'     trained.
#' @param partial A booled indicating whether the variables with zerp penalty
#'     loadings should be partialled out prior to the Lasso computation.
#' @param iter_resid A positive integer indicating the number of penalty loading
#'     iterations.
#' @param d A positive integer indicating the number of variables considered for
#'     the inital ols fit for computation of the penalty loadings.
#' @param HC_robust A boolean indicating whether penalty loadings should be
#'     computed under heteroskedasticity.
#'
#' @return \code{rlasso} returns an object of S3 class "\code{rlasso}".
#'
#' The function \code{predict} computes fitted values for a trained model of
#'     this class.
#'
#' An object of class "\code{rlasso}" is a list containig the following
#'     components:
#' \describe{
#' \item{\code{coef_lasso}}{A vector of Lasso coefficients.}
#' \item{\code{coef_postlasso}}{A vector of post Lasso ols coefficients.}
#' \item{\code{retained_X}}{A vector with indices of retained variables.}
#' \item{\code{nretained}}{The number of retained variables.}
#' \item{\code{y}}{A response vector.}
#' \item{\code{X}}{A feature matrix.}
#' \item{\code{mean_y}}{The mean of the response.}
#' \item{\code{mdl_glmnet}}{An S3 object of type "glmnet" with the trained
#'     Lasso model.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(100*100), 100, 100) # Simulate features
#' y <- 1 + X %*% (10*runif(100) * (runif(100) < 0.05)) + rnorm(100)
#' rlasso(y, X) # Compute rlasso
#'
#' @export rlasso
rlasso <- function(y, X,
                   include = NULL,
                   c = 1.1, gamma = 0.1 / log(nobs),
                   post = TRUE,
                   partial = FALSE,
                   iter_resid = 15, d = 5,
                   HC_robust = FALSE) {
  # Partial out controls or center y
  if (partial) {
    y <- y - predict(ols(y, X[, include], const = T))
    X <- cbind(1, apply(X[, -include], 2, function(x, X_) {
        x - predict(ols(x, X_))
      }, X_ = X[, include]))
    include <- 1 # always include a constant
  }#IF
  mean_y <- mean(y)
  y <- y - mean_y
  # Data parameters
  nobs <- length(y)
  ncol_X <- ncol(X)
  ninclude <- length(include)
  # Obtain penalty level and loadings via plug-in procedure
  penalty <- rlasso_penalty(y, X, include = include,
                            c = c, gamma = gamma,
                            iter_resid = iter_resid, d = d,
                            HC_robust = HC_robust)
  lambda <- penalty$lambda # penalty level
  psi <- penalty$psi # penalty loadings
  # Calculate glmnet rescaling of penalty loadings
  scale_psi <-  1 / mean(psi)
  # Compute Lasso on the first stage and retain non-zero features
  mdl_glmnet <- glmnet::glmnet(X, y,
                               family = "gaussian",
                               lambda = lambda / (2 * nobs * scale_psi),
                               penalty.factor = psi,
                               standardize = FALSE,
                               intercept = TRUE)
  coef_lasso <- mdl_glmnet$beta
  # Check how many variables were retained
  nonzero_coef <- c(1, which((mdl_glmnet$beta != 0)[, 1]))
  retained_X <- setdiff(nonzero_coef, 1)
  nretained <- length(retained_X)
  # Calculate post Lasso coefficients
  coef_postlasso <- Matrix::Matrix(0, ncol_X, 1)
  if(post){
    coef_postlasso[nonzero_coef] <- ols(y, X[, nonzero_coef])$coef
  }#IF
  # Return post Lasso coefficient estimate and number of retained instruments
  output <- list(coef_lasso = coef_lasso, coef_postlasso = coef_postlasso,
                 retained_X = retained_X, nretained = nretained,
                 y = y, X = X,
                 mean_y = mean_y,
                 mdl_glmnet = mdl_glmnet)
  class(output) <- "rlasso"
  return(output)
}#RLASSO

# Complementary methods ========================================================
#' Predict method for rlasso fits.
#'
#' Predict method for rlasso fits.
#'
#' @export predict.rlasso
predict.rlasso <- function(obj, newdata = NULL,
                          post = TRUE){
  # Get coefficients and data from lasso.fit
  if(post){
    coef <- obj$coef_postlasso
  }#IF
  if(is.null(newdata)){
    newdata <- obj$X
  }#IF
  # Calculate and return fitted values
  if(post){
    fitted <- newdata %*% coef
  } else {
    fitted <- predict(obj$mdl_glmnet, newdata)
  }#IFELSE
  # Re-center fitted values
  fitted <- as.matrix(fitted + obj$mean_y)
  return(fitted)
}#PREDICT.RLASSO

#' Instrument selection for rlasso fits.
#'
#' Instrument selection for rlasso fits. Returns \code{TRUE} if at least one
#'     instrument is retained..
#'
#' @export any_iv.rlasso
any_iv.rlasso <- function(obj, index_iv, ...){
  # Check whether any instruments are retained
  length(intersect(obj$retained, index_iv)) > 0
}#ANY_IV.RLASSO

# Complementary functions ======================================================
#' Compute plug-in penalty level and loadings for rlasso.
#'
#' Compute plug-in penalty level and loadings for rlasso.
#'
#' @param y A response vector.
#' @param X A feature matrix.
#' @param include An optional vector of feature indices with zero penalty
#'     loadings.
#' @param c A number strictly greater than 1 for the computation of the penalty
#'     level.
#' @param gamma A positive number for the computation of the penalty level.
#' @param iter_resid A positive integer indicating the number of penalty loading
#'     iterations.
#' @param d A positive integer indicating the number of variables considered for
#'     the inital ols fit for computation of the penalty loadings.
#' @param HC_robust A boolean indicating whether penalty loadings should be
#'     computed under heteroskedasticity.
#'
#' @return \code{rlasso_penalty} returns a list with the following components:
#' \describe{
#' \item{\code{lambda}}{The plug-in penalty level.}
#' \item{\code{psi}}{A vector of the plug-in penalty loadings.}
#' \item{\code{sigma}}{The estimated post-lasso residual variance.}
#' }
#'
#' @examples
#' # Add example here.
#'
#' @export rlasso_penalty
rlasso_penalty <- function(y, X,
                           include = NULL,
                           c = 1.1, gamma = 0.1/log(nobs),
                           iter_resid = 1, d = 5,
                           HC_robust = FALSE){
  # Data parameters
  nobs <- length(y)
  ncol_X <- ncol(X)
  ninclude <- length(include)
  if(ninclude == ncol_X) stop("Calculation of the plug-in peanlty terms is not
                              possible when no features can be excluded.")

  # Mean-center
  y <- y - mean(y)

  # Penalty level
  lambda <- 2 * c * sqrt(nobs) *
    qnorm(1 - (gamma / (2 * (ncol_X - ninclude))), 0, 1)

  # Iterative residual-estimation
  for (k in 0:iter_resid) {
    # Obtain kth-step residuals
    if (k == 0) {
      # For initial residuals, select the d features for which the absolute
      #     correlation with the response is the greatest. Only consider
      #     correlation for features that are not necessarily incldued in the
      #     model.
      if (d == 0) {
        resid.k <- response
      } else {
        minus_include <- setdiff(c(1:ncol_X), include)
        cor_y_x <- ccor(X[, minus_include, drop = F], y)
        sorder_cor_y_x <- order(abs(cor_y_x), decreasing = T)[1:d]
        # Select which features have greatest correlation coefficient and
        #     features that are necessarily included in the model.
        selected_X <- c(include, c(1:ncol_X)[minus_include][sorder_cor_y_x])
        # Calculate ols residuals using subset of selected features.
        resid_k <- y - predict(ols(y, X[, selected_X]))
      }#IFELSE
    } else {
      # Calculate glmnet rescaling of penalty factors
      scale_psi_k <- 1 / mean(psi_k)
      # Obtain Lasso residuals using the penalty level and loadings computed
      #     with the (k-1)th residuals
      mdl_glmnet_k <- glmnet::glmnet(X, y,
                                family = "gaussian",
                                lambda = lambda / (2 * nobs * scale_psi_k),
                                penalty.factor = psi_k,
                                standardize = FALSE,
                                intercept = TRUE)
      nonzero_coef <- c(1, which((mdl_glmnet_k$beta != 0)[, 1]))
      # Run post-Lasso ols and calculate residuals
      resid_k <- y - predict(ols(y, X[, nonzero_coef]))
    }#IFELSE

    # Calculate the kth-step penalty loadings for all features
    sigma_k <- sqrt(Matrix::mean(resid_k^2))
    if (HC_robust) {
      # W <- Diagonal(x=as.numeric(resid.k))
      # psi.k <- sqrt(diag(ccov(W%*%features, 0)))
      # psi.k <-  sqrt(as.matrix((1/N)*t(crossprod(resid.k^2, features^2)) -
      #     t((1/N)*crossprod(resid.k, features))^2)) # no DOF adjustment
      mean_X <- Matrix::colMeans(X)
      psi_k <- t(Matrix::crossprod(resid_k^2, X^2) / nobs +
                   (mean_X^2)*(sigma_k^2) -
                   2 * mean_X * Matrix::crossprod(resid_k^2, X) / nobs)
      if (all(psi_k>0)) {
        psi_k <- sqrt(psi_k)
      } else {
        warning('Possible issues with numerical zeros having negative signs.')
        psi_k <- sqrt(abs(psi_k))
      }#IFELSE

      # psi.k <- sqrt(crossprod(resid.k^2, scale(features, colMeans(features),
      #     FALSE)^2)/N)
    } else {
      # Penalty loadings under homoskedasticity (no dof adjustment)
      psi_k <- sigma_k * sqrt(Matrix::colMeans(X^2) - Matrix::colMeans(X)^2)
    } #IFELSE

    # Set penalty loadings to zero for those that are necessarily included
    if (!("matrix" %in% class(psi_k))) psi_k <- as.matrix(psi_k)
    psi_k[include] <- 0
  }#FOR

  # Organize and return output
  output <- list(lambda = lambda,
                 psi = psi_k,
                 sigma = sigma_k)
  return(output)
}#RLASSO_PENALTY
