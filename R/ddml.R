#' ddml: An implementation of LIE-conform Double/Debiased Machine Learning
#'
#' @docType package
#' @name ddml
#'
#' @examples
#'
#' # Demo using the Berry, Levinson, Pakes (1995) data, inspired by the empirical
#' # example of Chernozhukov, Hansen, and Spindler (2015) (CHS2015, hereafter).
#'
#' # Reproduce ols and tsls estimates from CHS2015
#' n <- length(BLP_1995$model.id)
#' y <- log(BLP_1995$share) -  log(BLP_1995$outshr)
#' x1 <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air", "mpd", "space", "price")]))
#'
#' # Construct demand instruments
#' X_ <- x1[, 1:5]; ncol_X <- ncol(X_) # exclude price
#' sum_other <- sum_rival <- matrix(0, n, 5)
#' for (i in 1:n) {
#'   other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
#'     (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
#'   rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
#'     (BLP_1995$cdid == BLP_1995$cdid[i])
#'   sum_other[i, ] <- colSums(X_[other_ind, , drop = FALSE])
#'   sum_rival[i, ] <- colSums(X_[rival_ind, , drop = FALSE])
#' }#FOR
#' Z_ <- cbind(sum_other, sum_rival); ncol_Z <- ncol(Z_)
#'
#' # ols and tsls. Note tsls is slightly different from the results reported in
#' # CHS2015 due to a slight instrument- construction error in the corresponding
#' # matlab code
#' ols(y, x1)$coef[6] # Gives essentially identical OLS results to BLP 1995 Table III
#' tsls(y, D = BLP_1995$price, Z = Z_, X = X_)$coef[1]
#'
#' # Construct additional instruments from interaction terms
#' tu = BLP_1995$trend/19;
#' mpdu = BLP_1995$mpd/7;
#' spaceu = BLP_1995$space/2;
#' XL_ <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air")], mpdu, spaceu, tu,
#'              BLP_1995$hpwt^2, BLP_1995$hpwt^3, mpdu^2, mpdu^3,
#'              spaceu^2, spaceu^3, tu^2, tu^3,
#'              BLP_1995$hpwt * BLP_1995$air,  mpdu * BLP_1995$air,
#'              spaceu * BLP_1995$air, tu * BLP_1995$air, BLP_1995$hpwt * mpdu,
#'              BLP_1995$hpwt * spaceu, BLP_1995$hpwt * tu, mpdu * spaceu,
#'              mpdu * tu, spaceu * tu))
#' ncol_XL <- ncol(XL_)
#' sum_otherL <- sum_rivalL <- matrix(0, n, 24)
#' for (i in 1:n) {
#'   other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
#'     (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
#'   rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
#'     (BLP_1995$cdid == BLP_1995$cdid[i])
#'   sum_otherL[i, ] <- colSums(XL_[other_ind, , drop = FALSE])
#'   sum_rivalL[i, ] <- colSums(XL_[rival_ind, , drop = FALSE])
#' }#FOR
#' ZL_ <- cbind(sum_otherL,sum_rivalL)
#' ncol_ZL <- ncol(ZL_)
#'
#' # tsls with full instrument set
#' tsls(y, D = BLP_1995$price, Z = ZL_, X = XL_)$coef[1]
#'
#' # It's possible to revisit the estimation with DDMML IV procedures. Some
#' # illustrative examples follow.
#
#' # Lasso DDML IV
#' # If a single model is considered for estimation, we may pass it as a single
#' # argument named "models". A Lasso estimator is considered here.
#' models <- list(what = mdl_glmnet,
#'                args = list(alpha = 1))
#' ddml_iv_fit <- ddml_iv(y, D = BLP_1995$price,
#'                        Z = ZL_, X = XL_,
#'                        models,
#'                        sample_folds = 2,
#'                        silent = T)
#' ddml_iv_fit$coef[1] # coefficient corresponding to D
#'
#' # Ensemble DDML IV
#' # When further felxibility regarding the functional form is desired, one may
#' # pass a multitude of models for residualization and the IV first stage.
#' # Depending on the type of ensemble, the procedure then combines the passed
#' # models such that the crossvalidation MSPE is minimized.
#'
#' # Definition of the considered models can be flexible. The below considers
#' # two sets of boosted trees, one Lasso regression, one ridge regression,
#' # one elastic net, and one ordinary least squares regression.
#' models_1 <- list(list(fun = mdl_xgboost,
#'                       args = list(num_parallel_tree = 1,
#'                                   nrounds = 6,
#'                                   max.depth = 1),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)),
#'                  list(fun = mdl_xgboost,
#'                       args = list(num_parallel_tree = 1,
#'                                   nrounds = 6,
#'                                   max.depth = 3),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)),
#'                  list(fun = mdl_glmnet,
#'                       args = list(alpha = 1),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)),
#'                  list(fun = mdl_glmnet,
#'                       args = list(alpha = 0),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)),
#'                  list(fun = mdl_glmnet,
#'                       args = list(alpha = 0.5),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)),
#'                  list(fun = ols,
#'                       args = list(),
#'                       assign_X = c(1:ncol_XL),
#'                       assign_Z = c(1:ncol_ZL)))
#
#' # One may also pass different sets of data to the procedures. For example, one
#' # may want to consider both the original set of instruments as well as the
#' # extended set of instruments. This is done below.
#' X_c <- cbind(X_, XL_); ncol_Xc <- ncol(X_c)
#' Z_c <- cbind(Z_, ZL_); ncol_Zc <- ncol(Z_c)
#' models_2 <- list(list(fun = mdl_xgboost,
#'                       args = list(num_parallel_tree = 1,
#'                                   nrounds = 6,
#'                                   max.depth = 1),
#'                       assign_X = c(1:ncol_X),
#'                       assign_Z = c(1:ncol_Z)), # initial set of instruments
#'                  list(fun = mdl_xgboost,
#'                       args = list(num_parallel_tree = 1,
#'                                   nrounds = 6,
#'                                   max.depth = 3),
#'                       assign_X = c(1:ncol_X),
#'                       assign_Z = c(1:ncol_Z)), # initial set of instruments
#'                  list(fun = mdl_glmnet,
#'                       args = list(alpha = 1),
#'                       assign_X = c((ncol_X + 1):ncol_Xc),
#'                       assign_Z = c((ncol_Z + 1):ncol_Zc)), # extended set
#'                  list(fun = ols,
#'                       args = list(),
#'                       assign_X = c(1:ncol_X),
#'                       assign_Z = c(1:ncol_Z)), # initial set
#'                  list(fun = ols,
#'                       args = list(),
#'                       assign_X = c((ncol_X + 1):ncol_Xc),
#'                       assign_Z = c((ncol_Z + 1):ncol_Zc))) # extended set
#
#' # Having specified the model arguments, an ensemble procedure may combine them
#' # through informed or uninformed (i.e., without crossvaldiation) averages, or
#' # simply select a single crossvalidation MSPE-minimizing model.
#'
#' # Consider both model specifications under alternative ensembles. Computation
#' # takes a few seconds.
#' ddml_iv_fit_1 <- ddml_iv(y, D = BLP_1995$price,
#'                        Z = Z_c, X = X_c,
#'                        models = models_1,
#'                        ens_type = c("cv"),
#'                        cv_folds = 5,
#'                        sample_folds = 2,
#'                        silent = T)
#' ddml_iv_fit_1$coef[1]
#
#' ddml_iv_fit_2 <- ddml_iv(y, D = BLP_1995$price,
#'                          Z = Z_c, X = X_c,
#'                          models = models_2,
#'                          ens_type = c("stacking"),
#'                          cv_folds = 5,
#'                          sample_folds = 2,
#'                          silent = T)
#' ddml_iv_fit_2$coef[1]
#'
#' # Equipped with the multiple coefficient estimates, we may recalculate the
#' # Number of products with inelastic demand as in CHS2015. BLP1995's estimates
#' suggest about 746. DDML_IV estimates suggest lower values. (Note the
#' estimates differ from those in CHS2015 quite substantially. This is likely
#' due to a construction error in the instruments.)
#' sum(ols(y, x1)$coef[6] * (BLP_1995$price) * (1 - BLP_1995$share) > -1)
#' sum(tsls(y, D = BLP_1995$price, Z = Z_, X = X_)$coef[1] * (BLP_1995$price) *
#'       (1 - BLP_1995$share) > -1)
#' sum(ddml_iv_fit_1$coef[1] * (BLP_1995$price) * (1 - BLP_1995$share) > -1)
#' sum(ddml_iv_fit_2$coef[1] * (BLP_1995$price) * (1 - BLP_1995$share) > -1)
#'
#' # Sidenote: Variation from random sample splits may be reduced when the arguments
#' # cv_folds and sample_folds are increased. Due to the computational burden, it
#' # may be advisable to utilize multiple cores. See a toy example below. (Should
#' # take roughly 20-40 seconds on a modern laptop.)
#'
#' # Computation with multiple cores
#' system.time({ddml_iv_fit_3 <- ddml_iv(y, D = BLP_1995$price,
#'                                       Z = Z_c, X = X_c,
#'                                       models = models_1,
#'                                       ens_type = c("stacking"),
#'                                       cv_folds = 10,
#'                                       sample_folds = 5,
#'                                       setup_parallel = list(type = "static",
#'                                                             cores = 2),
#'                                       silent = T)})
#
#' ddml_iv_fit_3$coef[1]
#'
#' # Computation one a single core
#' system.time({ddml_iv_fit_4 <- ddml_iv(y, D = BLP_1995$price,
#'                                       Z = Z_c, X = X_c,
#'                                       models = models_1,
#'                                       ens_type = c("stacking"),
#'                                       cv_folds = 10,
#'                                       sample_folds = 5,
#'                                       silent = T)})
#' ddml_iv_fit_4$coef[1]
#'
NULL
