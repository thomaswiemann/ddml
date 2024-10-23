test_that("ddml_policy computes with a single model", {
  # Simulate small dataset
  nobs <- 500
  X <- cbind(1, matrix(rnorm(nobs*39), nobs, 39))
  D_tld <-  X %*% runif(40) + rnorm(nobs)
  D <- 1 * (D_tld > mean(D_tld))
  y <- D + X %*% runif(40) + rnorm(nobs)
  policy <- sample(D, nobs)
  cluster_variable <- sample(1:100, nobs, replace = T)

  # Compute policy value
  learners = list(list(fun = ols))
  suppressWarnings({
    ddml_policy_fit <- ddml_policy(y, D, X,
                                   policy,
                                   policy_levels = c(0, 1),
                                   learners = learners,
                                   cv_folds = 3,
                                   sample_folds = 3,
                                   cluster_variable = cluster_variable,
                                   shortstack = T,
                                   silent = T)
  })#SUPPRESSWARNINGS

  # Compute policy value with re-scaled scores
  value_counts <- table(cluster_variable)
  count_vector <- value_counts[as.character(cluster_variable)]
  suppressWarnings({
    ddml_policy_fit_2 <- ddml_policy(y, D, X,
                                     policy,
                                     policy_levels = c(0, 1),
                                     learners = learners,
                                     cv_folds = 3,
                                     sample_folds = 3,
                                     cluster_variable = cluster_variable,
                                     shortstack = T,
                                     oos_pred = ddml_policy_fit$oos_pred,
                                     omega = 1/count_vector,
                                     silent = T)
  })#SUPPRESSWARNINGS

  # Compute rescaled score manually
  policy_value <- mean(aggregate(ddml_policy_fit$psi_b,
                                 by = list(cluster_variable),
                                 FUN = mean)[, 2])

  # Compute standard errors
  summary_res <- summary(ddml_policy_fit_2)

  # Compute standard errors manually for rescaled scores
  scores_i <- ddml_policy_fit$psi_b - policy_value
  scores <- aggregate(scores_i, by = list(cluster_variable),
                      FUN = mean)[, 2]
  psi_a <- aggregate(ddml_policy_fit$psi_a, by = list(cluster_variable),
                     FUN = mean)[, 2]
  nobs <- length(unique(cluster_variable))
  se <- sqrt(mean(scores^2) / nobs) / abs(mean(psi_a))

  # Compare point estimates and standard errors
  expect_equal(unname(round(ddml_policy_fit_2$policy_value, 2)),
               round(policy_value, 2))
  expect_equal(round(summary_res[2], 2), round(se, 2))
})#TEST_THAT
