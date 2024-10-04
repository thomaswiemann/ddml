test_that("ddml_policy computes with a single model", {

  library(dplyr)

  load("tests/testthat/stuff.RData")

  nobs <- 500
  X <- matrix(rnorm(nobs * 2), nobs, 2)

  pscore <-

  D <- df$PRICEBIN_1

  choice <- 1 * (df %>% select(CHOICE) == 1)
  mc <- product_data$price_list$regular - product_data$margin
  y <- choice * (D - mc) * 100

  X <- as.matrix(df %>%
                  select(CHAIN_71,
                         PRICE_2, PRICE_3,
                         FEATURE_2, FEATURE_3, #FEATURE_1,
                         DISPLAY_2, DISPLAY_3, # DISPLAY_1,
                         INCOME, FAMSIZE, RETIRED, UNEMPLOYED, SINGLEMOM, AGE, HIGHSCHOOL, COLLEGE,
                         WHTCLR, CHILDREN, MARRIED, DOGS, CATS, RENTER, TVS,
                         LOYALTY_1, LOYALTY_2, LOYALTY_3,
                         FREQUENCY_1, FREQUENCY_2, FREQUENCY_3,
                         TOTAL_VARIETY, TOTAL_SPEND))

  policy = df$xgb_c

  policy_levels <- c(3, 4.79)

  cluster_variable = df$PANID


  learners = list(list(fun = ols))
  learners_DX = learners
  sample_folds = 3
  ensemble_type = "nnls"
  shortstack = FALSE
  cv_folds = 10
  custom_ensemble_weights = NULL
  custom_ensemble_weights_DX = custom_ensemble_weights
  subsamples_byD = NULL
  cv_subsamples_byD = NULL
  y_X_p_res_list = NULL
  p_X_res_list = NULL
  trim = 0.01
  silent = FALSE


  ddml_policy_fit <- ddml_policy(y, D, X,
                                 policy,
                                 policy_levels,
                                 learners = learners,
                                 cv_folds = 3,
                                 sample_folds = 10,
                                 cluster_variable = cluster_variable,
                                 shortstack = T,
                                 silent = T)



  policy_value <- mean(aggregate(ddml_policy_fit$psi_b, by = list(cluster_variable),
                            FUN = mean)[, 2])
  scores_i <- ddml_policy_fit$psi_b - policy_value

  scores <- aggregate(scores_i, by = list(cluster_variable),
            FUN = mean)[, 2]
  psi_a <- aggregate(ddml_policy_fit$psi_a, by = list(cluster_variable),
                     FUN = mean)[, 2]
  nobs <- length(unique(cluster_variable))
  sqrt(mean(scores^2) / nobs) / abs(mean(psi_a))


  ddml_policy_fit$policy_value
  summary(ddml_policy_fit)

  object <- ddml_policy_fit
  object$ensemble_type

  ddml_policy_fit$psi_b
  ddml_policy_fit$psi_a

  # Define arguments
  expect_warning({
    ddml_policy_fit <- ddml_policy(y, D, X,
                                   policy,
                                   policy_levels,
                             learners = learners,
                             cv_folds = 3,
                             sample_folds = 3,
                             silent = T)
  })
  # Check output with expectations
  expect_equal(length(ddml_ate_fit$ate), 1)
})#TEST_THAT
