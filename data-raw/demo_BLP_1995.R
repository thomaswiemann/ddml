#' ddml Demo using the Berry, Levinson, Pakes (1995) data  =====================
#'
#' This R script is a demonstration of the ddml R package. We revisit the
#'     empirical example of Chernozhukov, Hansen, and Spindler (2015)
#'     (CHS2015, hereafter), which extends the instruments of Berry, Levinson,
#'     and Pakes (1995) (BLP1995, hereafter) and applies an instrument selection
#'     procedure based on the Lasso. We consider the same instrument extension
#'     and apply a selection of ensemble procedures that combines conventional
#'     linear estimators with computational alternatives including Lasso-based
#'     approaches and random forests.
#'
#' The BLP1995 data is included as an exemplary dataset with the ddml package
#'     under the name `BLP_1995'.

# Reproduction of LS and TSLS Estimates from CHS2015 ===========================

library(ddml)

# Data prepraration
n <- length(BLP_1995$model.id)
y <- as.matrix(log(BLP_1995$share) -  log(BLP_1995$outshr))
D <- as.matrix(BLP_1995$price)
x1 <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air", "mpd", "space", "price")]))
colnames(y) <- "share"
colnames(D) <- "price"

#' Construct the BLP1995 instruments. Instruments are sums of product
#'     characteristics (excluding price and other potentially endogenous
#'     variables) of other products offered by the firm as well as over
#'     competing firms. The exact instrument specification here follows the
#'     approach of CHS2015.
X_ <- x1[, 1:5]; ncol_X <- ncol(X_) # exclude price
sum_other <- sum_rival <- matrix(0, n, 5)
for (i in 1:n) {
  other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
  rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i])
  sum_other[i, ] <- colSums(X_[other_ind, , drop = FALSE])
  sum_rival[i, ] <- colSums(X_[rival_ind, , drop = FALSE])
}#FOR
Z_ <- cbind(sum_other, sum_rival); ncol_Z <- ncol(Z_)

#' Calculate LS and TSLS estimates corresponding to price. Note that the TSLS
#'     estimates differ from those in CHS2015. This is due to a slight
#'     instrument-construction error in the code of the CHS2015.
summary(ols(y, x1), type = "HC1")$res[6, ]
summary(tsls(y, D = BLP_1995$price, Z = Z_, X = X_), type = "HC1")$res[1, ]

# Extending the BLP1995 Instruments ============================================

#' Here, we follow the extension of the instruments as in CHS2015. In
#'     particular, various interaction terms are considered.
tu = BLP_1995$trend/19;
mpdu = BLP_1995$mpd/7;
spaceu = BLP_1995$space/2;
XL_ <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air")], mpdu, spaceu, tu,
                       BLP_1995$hpwt^2, BLP_1995$hpwt^3, mpdu^2, mpdu^3,
                       spaceu^2, spaceu^3, tu^2, tu^3, BLP_1995$hpwt *
                         BLP_1995$air,  mpdu * BLP_1995$air, spaceu *
                         BLP_1995$air, tu * BLP_1995$air, BLP_1995$hpwt *
                         mpdu, BLP_1995$hpwt * spaceu, BLP_1995$hpwt * tu,
                       mpdu * spaceu,  mpdu * tu, spaceu * tu))
ncol_XL <- ncol(XL_)
sum_otherL <- sum_rivalL <- matrix(0, n, 24)
for (i in 1:n) {
  other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
  rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i])
  sum_otherL[i, ] <- colSums(XL_[other_ind, , drop = FALSE])
  sum_rivalL[i, ] <- colSums(XL_[rival_ind, , drop = FALSE])
}#FOR
ZL_ <- cbind(sum_otherL,sum_rivalL); ncol_ZL <- ncol(ZL_)

# TSLS estimates on the extended set of instruments
summary(tsls(y, D = BLP_1995$price, Z = ZL_, X = XL_), type = "HC1")$res[1, ]

# Lasso DDML IV ================================================================

#' As a first pass, we revisit the estimation of the TSLS estimates using a
#'     simple DDML IV procedure based on a single Lasso estimator. In
#'     particular, we are interested in selection among instruments but not
#'     among controls. For this purpose, we specify two sets of models for the
#'     first and second stage computations, respectively.

# First stage model
models_FS <- list(what = rlasso,
                  args = list(include = c(1:ncol_XL)))

# Second stage model
models_SS <- list(what = ols,
                  args = list())

# DDML IV. We consider cross-residiualization across 10 sample folds here.
ddml_rlasso_fit <- ddml_iv(y, D = D,
                           Z = ZL_, X = XL_,
                           models = models_SS,
                           models_FS = models_FS,
                           sample_folds = 5,
                           silent = T)
ddml_rlasso_fit$coef[1]

# Export orthogonalized variables (optional)
if (FALSE) {
  export_ddml(ddml_rlasso_fit, "rlasso")
}#IF

# Neurel Net DDML IV ===========================================================

#' As an alternative, we may also consider neural networks. A wrapper for keras
#'     is already implemented. All that is left is to specify the architecture
#'     and pass it as an option to the wrapper. For specifying the architecture,
#'     keras high-level syntax is particularly intuitive. To use it, we
#'     explicitly load the keras library.

library(keras)

nnet <- keras::keras_model_sequential() %>%
  keras::layer_dense(units = 30, activation = "relu") %>%
  keras::layer_dense(units = 30, activation = "relu") %>%
  keras::layer_dense(units = 20, activation = "relu") %>%
  keras::layer_dense(units = 20, activation = "relu") %>%
  keras::layer_dense(units = 10, activation = "relu") %>%
  keras::layer_dense(units = 1)

# First stage model
model <- list(what = mdl_keras,
              args = list(model = nnet,
                          epochs = 200))

# DDML IV. We consider cross-residiualization across 10 sample folds here.
ddml_keras_fit <- ddml_iv(y, D = D,
                          Z = ZL_, X = XL_,
                          models = model,
                          sample_folds = 5,
                          silent = T)
ddml_keras_fit$coef[1]

#' Note that the architecture of a neural network is an intricate
#'     hyperparameter. specifying an ensemble of neural networks, where the
#'     choosen architecture performs well in cross-validation, may be helpful
#'     when prior knowledge is not available.

# Ensemble DDML IV =============================================================

#' A more ambitious approach may consider a variety of computational models for
#'     the estimation of the first and second stage. Here, we consider a
#'     combination of linear models, including unpenalized, Lasso and Ridge
#'     regression, random forests, generalized random forests
#'     (Athey et al., 2019), and boosted trees, using either the inital set of
#'     instruments as in BLP1995 or the extended set as in CHS2015. As before,
#'     we specify  seperate first and second stage model sets (although this is
#'     not strictly necessary).

#' Combine the data. It's convinient to specify sets of indices for easier
#'     column selection in the model specification. When using random forests,
#'     it's important that variables have columnames.
X_c <- cbind(X_, XL_); colnames(X_c) <- c(1:ncol(X_c))
Z_c <- cbind(Z_, ZL_); colnames(Z_c) <- c(1:ncol(Z_c))
set_X <- 1:ncol(X_); set_XL <- setdiff(c(1:ncol(X_c)), set_X)
set_Z <- 1:ncol(Z_); set_ZL <- setdiff(c(1:ncol(Z_c)), set_Z)

#' First stage models. When computing ensemble procedures, the control and
#'     instruments must be explicitely assigned via `assign_X' and `assign_Z'.
models_FS <- list(list(fun = ols,
                       args = list(),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = ols,
                       args = list(),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = rlasso,
                       args = list(include = 1:ncol(XL_)),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_glmnet,
                       args = list(alpha = 0.5),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_glmnet,
                       args = list(alpha = 0),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_xgboost,
                       args = list(num_parallel_tree = 1,
                                   nrounds = 6,
                                   max.depth = 3),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = mdl_randomForest,
                       args = list(ntree  = 100),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = mdl_grf,
                       args = list(num.trees  = 100),
                       assign_X = set_X,
                       assign_Z = set_Z))

# Second stage models. Here, we only omit the rlasso procedure, which selects
#'     only among instruments (and is thus identical to ols in the second
#'     stage).
models_SS <- list(list(fun = ols,
                       args = list(),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = ols,
                       args = list(),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_glmnet,
                       args = list(alpha = 0.5),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_glmnet,
                       args = list(alpha = 0),
                       assign_X = set_XL,
                       assign_Z = set_ZL),
                  list(fun = mdl_xgboost,
                       args = list(num_parallel_tree = 1,
                                   nrounds = 6,
                                   max.depth = 3),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = mdl_randomForest,
                       args = list(ntree  = 100),
                       assign_X = set_X,
                       assign_Z = set_Z),
                  list(fun = mdl_grf,
                       args = list(num.trees  = 100),
                       assign_X = set_X,
                       assign_Z = set_Z))

#' Models may be aggregated to ensembles in various ways. Here, we consider
#'     both uninformed averages across models as well as combinations chosen
#'     such that the crossvalidation MSPE is minimized.
ens_type <- c("stacking_01", "stacking_nn", "cv", "stacking", "average")

#' We may now estimate the DDML IV coefficient. Computation may take up to 1min.
ensemble_fit <- ddml_iv(y, D = D,
                        Z = Z_c, X = X_c,
                        models = models_SS,
                        models_FS = models_FS,
                        ens_type = ens_type,
                        cv_folds = 5,
                        sample_folds = 2,
                        silent = T)
ensemble_fit$coef

# Export orthogonalized variables (optional)
if (FALSE) {
  export_ddml(ensemble_fit, "ensemble")
}#IF

#' To better understand the composition of the ensembles, it is often useful to
#'     inspect the ensemble coefficients. These may readily be retrieved from
#'     the fitted object. Here, we see that the generalized random forest has
#'     been assigned the most weight in the crossvalidation informed ensemble
#'     procedures (all excpet `average').
ensemble_fit$weights


# Elasticities =================================================================

#' Equipped with the multiple coefficient estimates, we may recalculate the
#'     number of products with inelastic demand as in CHS2015. BLP1995's
#'     estimates suggest about 746. CHS2015's results suggests that there is
#'     inelastic demand for only 305 products when variables are selected from
#'     the extended set. (Note: these numbers differ from those reported in
#'     CHS2015 due to the abovementioned construction error.) The ensemble DDML
#'     IV estimates suggest that there is inelastic demand for more products
#'     than suggested by the BLP1995 estimates.

# BLP1995 OLS implied number of products with inelastic demand
sum(ols(y, x1)$coef[6] * (BLP_1995$price) *
      (1 - BLP_1995$share) > -1)

# BLP1995 TSLS implied number of products with inelastic demand
sum(tsls(y, D = BLP_1995$price, Z = Z_, X = X_)$coef[1] * (BLP_1995$price) *
      (1 - BLP_1995$share) > -1)

# Ensemble DDML IV implied number of products with inelastic demand
sum(ensemble_fit$coef[4] * (BLP_1995$price) * (1 - BLP_1995$share) > -1)

# Lasso DDML IV implied number of products with inelastic demand
sum(ddml_rlasso_fit$coef * (BLP_1995$price) * (1 - BLP_1995$share) > -1)
