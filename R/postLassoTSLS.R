#' An implementation of post-Lasso TSLS.
#'
#' An implementation of post-Lasso TSLS.
#'
#' @export postLassoTSLS
postLassoTSLS <- function(response, control, D, instrument,
                          splitSample = FALSE,
                          splitSample.IVChoice = FALSE,
                          X.select = FALSE,
                          penaltyMethod = 'r',
                          K.resid = 15, d = 5,
                          c = 1.1, gamma = 0.1/log(N),
                          heteroskedasticity = FALSE,
                          lambda = NULL, psi = NULL,
                          lambda.list = NULL, K.CV = 10,
                          cores = 1,
                          nullEst = FALSE){
  # Data parameters
  N <- length(response)
  ncontrol <- ncol(control)
  ninstrument <- ncol(instrument)
  lambda.fixed <- lambda

  # Adjust default gamma in case of split sample TSLS
  if((splitSample | splitSample.IVChoice) & gamma == 0.1/log(N))
    gamma <- 0.1/log(N/(1+splitSample+splitSample.IVChoice))

  # Draw sub-samples for split sample and IV choice, if necessary
  sample <- sample.IVChoice <- replicate(4, c(1:N), simplify=F) # full samples
  if(splitSample){
    sample[[1]] <- sample(c(1:N), floor(N/2), replace = F)
    sample[[2]] <- c(1:N)[-sample[[1]]]
  }#IF
  for(s in 1:(splitSample+1)){
    split.N <-length(sample[[s]]) # half-sample size adjustment when necessary
    if(splitSample.IVChoice){
      sample.IVChoice[[1+2*(s-1)]] <- sample(c(1:split.N), floor(split.N/2), replace = F)
      sample.IVChoice[[2+2*(s-1)]] <- c(1:split.N)[-sample.IVChoice[[1+2*(s-1)]]]
    } else {
      sample.IVChoice[[1+2*(s-1)]] <- sample.IVChoice[[2+2*(s-1)]] <- c(1:split.N)
    } #IFELSE
  }#FOR

  # Loop over sub samples for split sample estimates
  coef.pLTSLS <- nRetained.Z <- nRetained.X <- matrix(0, 2^(splitSample+splitSample.IVChoice), 1)
  for(s in 1:(1+splitSample)){
    # Loop over sub samples for split sample IV choice estimates
    for(s.IV in 1:(1+splitSample.IVChoice)){
      # Center control, and instrument using only sub-sample under consideration
      features.s1 <- cbind(control, instrument)[sample[[s]], ][sample.IVChoice[[s.IV + 2*(s-1)]], ]

      # Center D (the endogeneous variable of interest)
      D.s1 <- D[sample[[s]]][sample.IVChoice[[s.IV + 2*(s-1)]]]
      mean.D.s1 <- mean(D.s1); D.s1 <- D.s1 -  mean.D.s1
      N.s1 <- length(D.s1)

      # Obtain penalty level and loadings
      if(penaltyMethod=='r'){
        # Plug-in penalty always recalculates penalty level and loadings
        penalty <- calcPenalty(D.s1, features.s1,
                               include = c(1:ncontrol),
                               heteroskedasticity,
                               K.resid, d,
                               c, gamma,
                               lambda = lambda.fixed)
        if(is.null(lambda.fixed)) lambda <- penalty$lambda
        psi <- penalty$psi
      } else if(penaltyMethod=='CV'){
        # Cross-validated penalty level may be given or may be re-computed
        if(is.null(lambda.fixed)){
          if(is.null(lambda.list)) stop("Please supply a grid of penalty levels when recomputing the CV-penalty.")
          lambda <- penalty.CV(response = D.s1, features = features.s1,
                               include = c(1:ncontrol),
                               lambda.list = lambda.list,
                               heteroskedasticity = heteroskedasticity,
                               K.resid = K.resid, d = d,
                               K.CV = K.CV,
                               test.run = TRUE, cores = cores)$lambda
        }#IF
        # Calculate penalty loadings on the full sample given a CV-penalty level
        if(is.null(psi)){
          psi <- calcPenalty(D.s1, features.s1,
                             include = c(1:ncontrol),
                             heteroskedasticity,
                             K.resid, d,
                             c, gamma,
                             lambda)$psi
        }#IF
      }#ELSEIF

      # Calculate glmnet rescaling of penalty factors
      scale.psi <-  1/mean(psi)
      # Compute Lasso on the first stage and retain non-zero features
      lasso.firstStage <- glmnet::glmnet(features.s1, D.s1,
                                         family = "gaussian",
                                         lambda = lambda/(2*N.s1*scale.psi),
                                         penalty.factor = psi,
                                         standardize = FALSE,
                                         intercept = TRUE)
      retain <- c(1, which(as.matrix(!(lasso.firstStage$beta == 0))))

      #print(colnames(features.s1)[setdiff(retain, c(1:ncontrol))])
      nRetained.Z[(s.IV + 2*(s-1))*(splitSample.IVChoice) +
                    s*(!splitSample.IVChoice)] <- length(setdiff(retain, c(1:ncontrol))) # how many of the retained variables are not controls
      if(X.select){
        nRetained.X[(s.IV + 2*(s-1))*(splitSample.IVChoice) +
                      s*(!splitSample.IVChoice)] <- length(setdiff(retain, c((ncontrol+1):(ncontrol+ninstrument))))
      }#IF

      # Center control, instrument, and D in second sub-sample for split sample IV choice
      if(splitSample.IVChoice){
        features.s2 <- cbind(control, instrument)[sample[[s]], ][sample.IVChoice[[2^(s.IV==1) + 2*(s-1)]], ]
        features.s2 <- features.s2[, retain] # remove omitted and non-retained features
        D.s2 <- D[sample[[s]]][sample.IVChoice[[2^(s.IV==1) + 2*(s-1)]]] -  mean.D.s1
        control.s2 <- Matrix::Matrix(control[sample[[s]], ])[sample.IVChoice[[2^(s.IV==1) + 2*(s-1)]], ]
      } else {
        features.s2 <- features.s1[, retain]
        D.s2 <- D.s1
        control.s2 <-  Matrix::Matrix(control[sample[[s]], ])[sample.IVChoice[[s.IV + 2*(s-1)]], ]
      }#IFELSE

      # Remove redundant variables from memory
      rm(list=c("features.s1", "D.s1", "lasso.firstStage"))

      # Finish iteration if no instruments are retained (the associated coefficient is 0)
      if(nRetained.Z[(s.IV + 2*(s-1))*(splitSample.IVChoice) +
                     s*(!splitSample.IVChoice)]==0 & !nullEst) next

      # Compute post Lasso first stage
      ZZ1 <- as.matrix(crossprod(features.s2))
      DZ1 <- Matrix::crossprod(cbind(D.s2, control.s2), features.s2)
      FS1 <- Matrix::tcrossprod(csolve(ZZ1), DZ1)

      # Adjust control, instrument, and response in the second half-sample when necessary
      if(splitSample){
        features.s3 <- cbind(control, instrument)[sample[[2^(s==1)]], ]
        features.s3 <- features.s3[, retain] # remove omitted and non-retained features
        D.s3 <- D[sample[[2^(s==1)]]] -  mean.D.s1
        control.s3 <- control[sample[[2^(s==1)]], ]
        response.s3 <- response[sample[[2^(s==1)]]]
      } else if (splitSample.IVChoice){
        features.s3 <- features.s2
        D.s3 <-  D.s2
        control.s3 <- control.s2
        response.s3 <- response[sample[[1]]][sample.IVChoice[[2^(s.IV==1)]]]
      } else {
        features.s3 <- features.s2
        D.s3 <-  D.s2
        control.s3 <- control.s2
        response.s3 <- response[sample[[1]]][sample.IVChoice[[1]]]
      }#IFELSE

      # Calculate SSIV coefficient
      DZ2 <- Matrix::crossprod(cbind(D.s3, control.s3), features.s3)
      ZY2 <- Matrix::crossprod(features.s3, response.s3)
      TSLS.coef <- Matrix::crossprod(csolve(as.matrix(DZ2%*%FS1)),
                                     Matrix::crossprod(FS1, ZY2))
      coef.pLTSLS[(s.IV + 2*(s-1))*(splitSample.IVChoice) +
                    s*(!splitSample.IVChoice)] <- TSLS.coef[1]
    }#FOR
  }#FOR

  # Combine half-sample TSLS estimates
  coef.pLTSLS <- mean(coef.pLTSLS) # A&F (2020) simply average the estimates

  # Return post Lasso coefficient estimate and number of retained instruments
  output <- list(coef = coef.pLTSLS,
                 nRetained.Z = mean(nRetained.Z),
                 anyRetained = all(nRetained.Z>0),
                 nRetained.X = mean(nRetained.X))
  return(output)
}#POSTLASSOTSLS

calcPenalty <- function(response, features,
                        include = NULL, # variables that are always included in the model
                        heteroskedasticity = FALSE,
                        K.resid = 15, d = 5,
                        c = 1.1, gamma = 0.1/log(N),
                        lambda = NULL){
  # Data parameters
  N <- length(response)
  nfeatures <- ncol(features)
  ninclude <- length(include)
  if(ninclude == nfeatures) stop("Calculation of peanlty level and loadings not possible when no features can be excluded.")

  # Mean-center
  response <- response - mean(response)

  # Penalty level
  if(is.null(lambda)) lambda <- 2*c*sqrt(N)*qnorm(1-(gamma/(2*(nfeatures-ninclude))), 0, 1)

  # Iterative residual-estimation
  for(k in 0:K.resid){
    # Obtain kth-step residuals
    if(k == 0){
      # For initial residuals, select the d features for which the absolute correlation with the response is the greatest
      # Only consider correlation for features that are not necessarily incldued in the model
      # Calucluate correlation coefficient between response and (included) features
      if(d==0){
        resid.k <- response
      } else {
        minus.include <- setdiff(c(1:nfeatures), include)
        cor.y.x <- ccor(features[, minus.include], response)
        # Select which features have greatest correlation coefficient + features that are necessarily included in the model
        selected.features <- c(include, c(1:nfeatures)[minus.include][order(abs(cor.y.x), decreasing = TRUE)[1:d]])
        # Calculate OLS residuals using subset of selected features
        resid.k <- response - predict(ols(response, features[, selected.features]))
      }#IFELSE
    } else {
      # Calculate glmnet rescaling of penalty factors
      scale.psi.k <- 1/mean(psi.k)
      # Obtain LASSO residuals using the penalty level and loadings computed with the (k-1)th residuals
      lasso.k <- glmnet::glmnet(features, response,
                                family = "gaussian",
                                lambda = lambda/(2*N*scale.psi.k),
                                penalty.factor = psi.k,
                                standardize = FALSE,
                                intercept = TRUE)
      retain <- c(1, which(as.matrix(!(lasso.k$beta == 0))))
      # Run post Lasso OLS and calculate residuals
      resid.k <- response - predict(ols(response, features[, retain]))
    }#IFELSE

    # Calculate the kth-step penalty loadings for all features
    sigma.k <- sqrt(mean(resid.k^2))
    if(heteroskedasticity == TRUE){
      # W <- Diagonal(x=as.numeric(resid.k))
      # psi.k <- sqrt(diag(ccov(W%*%features, 0)))
      #psi.k <-  sqrt(as.matrix((1/N)*t(crossprod(resid.k^2, features^2)) - t((1/N)*crossprod(resid.k, features))^2)) # no DOF adjustment
      features.mean <- colMeans(features)
      psi.k <-  t(crossprod(resid.k^2, features^2)/N + (features.mean^2)*(sigma.k^2) - 2*features.mean*crossprod(resid.k^2, features)/N)
      if(all(psi.k>0)){
        psi.k <- sqrt(psi.k)
      } else {
        warning('Possible issues with numerical zeros having negative signs. Double check.')
        psi.k <- sqrt(abs(psi.k))
      }#IFELSE

      #psi.k <- sqrt(crossprod(resid.k^2, scale(features, colMeans(features), FALSE)^2)/N)
    } else {
      # Penalty loadings under homoskedasticity
      #psi.k <- sigma.k*sqrt(diag(ccov(features, 0)))
      psi.k <- sigma.k * sqrt(colMeans(features^2) - colMeans(features)^2) # no DOF adjustment
      #psi.k <- sigma.k * sqrt(colMeans(scale(features, colMeans(features), FALSE)^2))
    } #IFELSE

    # Set penalty loadings to zero for those that are necessarily included
    if(class(psi.k)!='matrix') psi.k <- as.matrix(psi.k)
    psi.k[include] <- 0

  }#FOR

  # Calculate residual variance if not yet defined
  if(!exists("sigma.k")) sigma.k <- sqrt(mean(resid.k^2))

  # Return final penalty level and loadings
  output <- list(lambda = lambda,
                 psi = psi.k,
                 sigma = sigma.k)
  return(output)
}#CALCPENALTY
