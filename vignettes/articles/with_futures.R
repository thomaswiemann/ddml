# Example script demonstrating parallel cross-validation with base parallel

# Load necessary packages (assuming ddml is loaded/installed)
library(parallel) # For parallel execution functions used by crossval

library(devtools)
load_all()

# Ensure the crossval function and necessary learners (like ols) are available
# If running outside the ddml package context, you might need:
# source("R/crossval.R")
# source("R/learners.R") # Or wherever ols/base learners are defined

# --- Example Setup ---
set.seed(123) # for reproducibility

# Construct variables from the included Angrist & Evans (1998) data
y = AE98[, "worked"]
D = AE98[, "morekids"]
Z = AE98[, "samesex"]
X = AE98[, c("age","agefst","black","hisp","othrace","educ")]

# Define custom learner
library(gbm)

# Write gbm wrapper
mdl_gbm <- function(y, X, ...) {
  gbm_fit <- gbm::gbm.fit(x = X, y = y, ...) # fit gbm
  class(gbm_fit) <- c("mdl_gbm", class(gbm_fit)) # append class
  return(gbm_fit) # return fitted gbm object
}#MDL_GBM

# Write prediction method for gbm.wrapper
predict.mdl_gbm <- function(object, newdata, ...) {
  class(object) <- class(object)[-1] # remove mdl_gbm from the class list
  predict(object, newdata, type = "response", n.trees = object$n.trees)
}#MDL_GBM


# Define arguments for crossval
# Example learners: OLS on all vars, OLS on subset known to be relevant
# (In practice, subset selection would be part of the learning process)
learners <- list(
  list(fun = ols), # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 1,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 3,
                    verbose = FALSE)), # OLS using intercept + true non-zero vars, # OLS using all variables
  list(fun = mdl_gbm,
        args = list(distribution  = "bernoulli",
                    n.trees = 500,
                    interaction.depth = 5,
                    verbose = FALSE)) # OLS using intercept + true non-zero vars
)

# --- Sequential Execution ---
message("Running cross-validation sequentially...")
set.seed(456) # Set seed before sequential run
time_seq <- system.time({
  cv_res_seq <- crossval(y, X, Z = NULL,
                         learners = learners,
                         cv_folds = 2,
                         parallel = FALSE,
                         silent = FALSE, # Show progress
                         progress = "Sequential CV: ")
})
print(cv_res_seq$mspe)
print(time_seq)

# --- Parallel Execution using base `parallel` ---

message("\nRunning cross-validation in parallel using 2 cores...")
set.seed(456) # Reset seed to get same folds as sequential run
time_par <- system.time({
  cv_res_par <- crossval(y, X, Z = NULL,
                         learners = learners,
                         cv_folds = 2,
                         parallel = TRUE,
                         num.cores = 4,
                         parallel.export = c("mdl_gbm", "predict.mdl_gbm"),
                        #  parallel.packages = "gbm",
                         silent = FALSE # Progress reporting disabled
                         )
})
print(cv_res_par$mspe)
print(time_par)

# Verify results are numerically close (within tolerance)
print(all.equal(cv_res_seq$mspe, cv_res_par$mspe))
print(all.equal(cv_res_seq$oos_resid, cv_res_par$oos_resid))

# --- Clean up ---
# No explicit cleanup needed for stopCluster, handled by on.exit in crossval
message("\nParallel execution complete.")
