#' Estimator of the Mean Squared Prediction Error using Cross-Validation.
#'
#' @family utilities
#'
#' @description Estimator of the mean squared prediction error of
#'     different learners using cross-validation.
#'
#' @inheritParams ddml_plm
#' @param y The outcome variable.
#' @param X A (sparse) matrix of predictive variables.
#' @param Z Optional additional (sparse) matrix of predictive variables.
#' @param learners \code{learners} is a list of lists, each containing four
#'     named elements:
#'     \itemize{
#'         \item{\code{fun} The base learner function. The function must be
#'             such that it predicts a named input \code{y} using a named input
#'             \code{X}.}
#'         \item{\code{args} Optional arguments to be passed to \code{fun}.}
#'         \item{\code{assign_X} An optional vector of column indices
#'             corresponding to variables in \code{X} that are passed to
#'             the base learner.}
#'         \item{\code{assign_Z} An optional vector of column indices
#'             corresponding to variables in \code{Z} that are passed to
#'             the base learner.}
#'     }
#'     Omission of the \code{args} element results in default arguments being
#'     used in \code{fun}. Omission of \code{assign_X} (and/or \code{assign_Z})
#'     results in inclusion of all predictive variables in \code{X} (and/or
#'     \code{Z}).
#' @param cv_folds Number of folds used for cross-validation.
#' @param cv_subsamples List of vectors with sample indices for
#'     cross-validation.
#' @param parallel Logical indicating whether to use parallel processing.
#' @param num.cores Integer number of cores to use for parallel execution.
#'     Defaults to \code{NULL}, which uses \code{parallel::detectCores() - 1}.
#' @param parallel.export Optional character vector of names of objects
#'     (functions or data) in the calling environment to be explicitly
#'     exported to the parallel workers. This is necessary for globally defined
#'     helper functions or S3 methods used by custom learners.
#'     Example: `parallel.export = c("my_custom_learner", "predict.my_custom_learner")`.
#' @param parallel.packages Optional character vector of additional package names
#'     to be loaded on the parallel workers (using `library()`). Required if
#'     custom learners or their methods depend on functions from these packages.
#'     Example: `parallel.packages = "gbm"`.
#' @param progress String to print before learner and cv fold progress (only
#'     applies when `parallel = FALSE`).
#'
#' @section Parallel Processing using the `parallel` package:
#' When \code{parallel = TRUE}, this function uses the base R `parallel` package
#' via a PSOCK cluster to distribute the computation across all combinations of
#' learners and cross-validation folds. This requires explicit exporting of
#' necessary objects and loading of packages onto the worker nodes.
#' \itemize{
#'     \item **Number of Cores:** Controlled by `num.cores`. Defaults to one less
#'         than the detected number of cores.
#'     \item **Custom Objects/Functions:** User-defined learners, helper
#'         functions, or S3 methods needed by learners *must* be specified by name
#'         in the `parallel.export` argument. Ensure these objects exist in the
#'         environment where `crossval` is called (typically the global environment).
#'     \item **Packages:** Packages required by custom learners (e.g., `gbm` if
#'         using `gbm::gbm.fit`) *must* be specified in the `parallel.packages`
#'         argument. Base packages `stats` and `methods` are loaded automatically.
#' }
#' Reproducibility is handled using `parallel::clusterSetRNGStream`. The user
#' should set the global random seed *before* calling `crossval` for reproducible
#' results in both sequential and parallel modes. Output from parallel workers is
#' typically not displayed; progress reporting is disabled in parallel mode. If
#' parallel setup fails, execution falls back to sequential processing with a
#' warning.
#'
#' @return \code{crossval} returns a list containing the following components:
#'     \describe{
#'         \item{\code{mspe}}{A vector of MSPE estimates,
#'             each corresponding to a base learners (in chronological order).}
#'         \item{\code{oos_resid}}{A matrix of out-of-sample prediction errors,
#'             where rows correspond to observations (ordered by their original
#'             index) and columns correspond to base learners (in chronological
#'             order).}
#'         \item{\code{cv_subsamples}}{Pass-through of \code{cv_subsamples}.
#'             See above.}
#'     }
#' @export
#'
#' @examples
#' # Construct variables from the included Angrist & Evans (1998) data
#' y = AE98[, "worked"]
#' X = AE98[, c("morekids", "age","agefst","black","hisp","othrace","educ")]
#'
#' # Compare ols, lasso, and ridge using 4-fold cross-validation (sequentially)
#' set.seed(1) # Set seed for reproducibility
#' cv_res_seq <- crossval(y, X,
#'                        learners = list(list(fun = ols),
#'                                        list(fun = mdl_glmnet),
#'                                        list(fun = mdl_glmnet,
#'                                             args = list(alpha = 0))),
#'                        cv_folds = 4,
#'                        parallel = FALSE, # Explicitly sequential
#'                        silent = FALSE,
#'                        progress = "Sequential CV: ")
#' cv_res_seq$mspe
#'
#' # --- Example using parallel package ---
#'
#' # Define a custom learner (e.g., requires 'gbm' package and S3 method)
#' library(gbm)
#' mdl_gbm <- function(y, X, ...) {
#'   fit <- gbm::gbm.fit(x = X, y = y, ...)
#'   class(fit) <- c("mdl_gbm", class(fit)) # Assign custom class
#'   return(fit)
#' }
#' # Define predict method for the custom class
#' predict.mdl_gbm <- function(object, newdata, ...) {
#'   class(object) <- class(object)[-1L] # Remove custom class before calling base predict
#'   stats::predict(object, newdata = newdata, type = "response",
#'                  n.trees = object$n.trees, ...)
#' }
#'
#' learners_custom <- list(list(fun = ols), # ols from ddml package
#'                         list(fun = mdl_gbm, # Function defined globally
#'                              args = list(n.trees = 100, verbose = FALSE,
#'                                          distribution = "bernoulli")))
#'
#' # Run in parallel using 2 cores
#' # Export the globally defined functions, load the needed package 'gbm'
#' set.seed(1) # Reset seed for reproducibility with parallel RNG streams
#' cv_res_par <- crossval(y = y, X = X,
#'                        learners = learners_custom,
#'                        cv_folds = 4,
#'                        parallel = TRUE,
#'                        num.cores = 2,
#'                        parallel.export = c("mdl_gbm", "predict.mdl_gbm"),
#'                        parallel.packages = "gbm",
#'                        silent = TRUE) # Progress reporting disabled in parallel
#' print("Parallel MSPE:")
#' print(cv_res_par$mspe)
#'
#' # Compare with sequential (set seed again for identical fold generation)
#' set.seed(1)
#' cv_res_seq_custom <- crossval(y = y, X = X,
#'                        learners = learners_custom,
#'                        cv_folds = 4,
#'                        parallel = FALSE,
#'                        silent = TRUE)
#' print("Sequential MSPE:")
#' print(cv_res_seq_custom$mspe)
#' print(all.equal(cv_res_seq_custom$mspe, cv_res_par$mspe))
#' print(all.equal(cv_res_seq_custom$oos_resid, cv_res_par$oos_resid))
crossval <- function(y, X, Z = NULL,
                     learners,
                     cv_folds = 5,
                     cv_subsamples = NULL,
                     parallel = FALSE,
                     num.cores = NULL,
                     parallel.export = NULL,
                     parallel.packages = NULL,
                     silent = FALSE,
                     progress = NULL) {
  # Data parameters
  # Ensure y is a vector
  if (!is.vector(y) && !is.factor(y)) y <- as.vector(y)
  nobs <- length(y)
  nlearners <- length(learners)

  # Create cv sample fold tuple
  if (is.null(cv_subsamples)) {
    # Note: generate_subsamples uses the global RNG state
    # Ensure seed is set *before* calling crossval for reproducibility
    cv_subsamples <- generate_subsamples(nobs, cv_folds)
  }#IF
  cv_folds <- length(cv_subsamples)
  # Calculate total number of observations used across all folds
  # Handles potential variation if cv_subsamples provided manually
  nobs_per_fold <- lengths(cv_subsamples)
  nobs_used <- sum(nobs_per_fold)

  # Define the function to be applied over iterations
  # Adjusted to handle mapping from linear index to learner/fold
  cv_fun <- function(iter_idx) {
    # Select model and cv-fold for this job
    j <- ceiling(iter_idx / cv_folds) # jth model
    i <- (iter_idx - 1) %% cv_folds + 1 # ith CV fold (1-based)
    fold_x <- cv_subsamples[[i]]

    # Print progress (only if sequential and not silent)
    # This provides fine-grained progress within sequential chunks
    if (!parallel && !silent) {
      cat(paste0("\r", progress,
                 "learner ", j, "/", nlearners,
                 ", cv fold ", i, "/", cv_folds, "          "))
    }

    # Compute model for this fold
    crossval_compute(test_sample = fold_x,
                     learner = learners[[j]],
                     y = y, X = X, Z = Z)
  }

  # Setup parallel cluster (optional)
  cl <- NULL # Initialize cluster object
  if (parallel) {
    # Define essential variables to export from this function's environment
    vars_to_export <- c("y", "X", "learners", "cv_subsamples",
                        "crossval_compute", "cv_folds", "nlearners", "parallel")
    if (!is.null(Z)) vars_to_export <- c(vars_to_export, "Z")

    # Attempt to set up the cluster
    # User exports (`parallel.export`) are fetched from the parent frame
    cl <- setup_parallel_cluster(num.cores = num.cores,
                                 vars_to_export = vars_to_export,
                                 export_envir = environment(),
                                 user_export_vars = parallel.export,
                                 user_export_envir = parent.frame(), # Where crossval is called
                                 packages_to_load = parallel.packages)

    # If setup successful, register stopCluster on exit
    if (!is.null(cl)) {
      on.exit(parallel::stopCluster(cl), add = TRUE)
      if (!silent) message("Parallel backend set up using ", length(cl), " cores.")
    } else {
      # Fallback to sequential if setup failed (warning given by setup_parallel_cluster)
      parallel <- FALSE # Ensure sequential execution path is taken
    }
  }# IF parallel

  # --- Execute Calculation ---
  total_tasks <- cv_folds * nlearners
  results_list <- vector("list", total_tasks) # Pre-allocate list

  # Define the task indices (all tasks)
  task_indices <- 1:total_tasks

  if (!is.null(cl)) { # Parallel execution
    if (!silent) message("Starting parallel cross-validation...")
    par_results <- tryCatch({
      parallel::parSapply(cl, task_indices, cv_fun, simplify = FALSE)
    }, error = function(e) {
      # Ensure cluster is stopped immediately on error
      parallel::stopCluster(cl)
      on.exit(NULL) # Remove the original on.exit handler
      stop(sprintf("Error during parallel execution: %s", e$message), call. = FALSE)
    })

    # parSapply with simplify = FALSE returns a list
    results_list <- par_results
    if (!silent) message("Parallel cross-validation finished.")

  } else { # Sequential execution
    if (!silent && total_tasks > 0) cat(progress) # Initial progress message
    # Use lapply for sequential execution to always get a list back
    # Error handling within lapply loop is less clean, rely on cv_fun/crossval_compute
    seq_results <- tryCatch({
        lapply(task_indices, cv_fun)
    }, error = function(e) {
        # Add context to the error if it bubbles up from cv_fun/crossval_compute
        stop("Error during sequential execution: ", e$message, call. = FALSE)
    })
    results_list <- seq_results
    if (!silent && total_tasks > 0) cat("\n") # Final newline after sequential progress
  }

  # --- Compile Results --- (results_list contains residuals for each task)
  # Check if all elements are suitable before proceeding
  if (length(results_list) != total_tasks) {
    stop(paste("Internal error: Unexpected number of results returned.",
               "Expected:", total_tasks, "Got:", length(results_list)))
  }
  # Further checks can be added here (e.g., are all elements numeric vectors?)

  # Determine expected length per fold/learner based on cv_subsamples
  # Map linear index back to fold index i to get expected length
  expected_lengths <- nobs_per_fold[((task_indices - 1) %% cv_folds) + 1]

  # Verify lengths of returned residual vectors
  actual_lengths <- lengths(results_list)
  if (!all(actual_lengths == expected_lengths)) {
      mismatched_idx <- which(actual_lengths != expected_lengths)[1]
      stop(paste("Internal error: Mismatch in residual vector lengths.",
                 "Task", mismatched_idx, "Expected:", expected_lengths[mismatched_idx],
                 "Got:", actual_lengths[mismatched_idx]))
  }

  # Concatenate all residuals (already ordered by task: L1F1, L1F2.. L2F1, L2F2..)
  oos_resid_flat <- unlist(results_list, use.names = FALSE)

  # Reshape into matrix: rows = observations within folds, cols = learners
  # The order needs careful handling. We have blocks of folds per learner.
  # We need to interleave observations correctly to get obs x learner matrix.
  oos_resid_mat <- matrix(oos_resid_flat, nrow = nobs_used, ncol = nlearners)

  # The matrix currently has rows ordered like:
  # [obs_fold1, obs_fold2, ..., obs_foldK] for Learner 1
  # [obs_fold1, obs_fold2, ..., obs_foldK] for Learner 2
  # ...
  # We need to reorder rows to match the original observation indices.
  # 'original_indices_in_folds' lists the original obs index for each position
  # across the concatenated folds [fold1, fold2, ..., foldK].
  original_indices_in_folds <- unlist(cv_subsamples, use.names = FALSE)
  # 'row_order' gives the position in the concatenated list for each sorted original index.
  # Applying this order sorts the rows of oos_resid_mat according to the original obs index.
  row_order <- order(original_indices_in_folds)

  # Apply this order to the rows of the matrix
  oos_resid_ordered <- oos_resid_mat[row_order, , drop = FALSE]

  # Compute MSPE by learner (columns)
  mspe <- colMeans(oos_resid_ordered^2)

  # Organize and return output
  output <- list(mspe = mspe,
                 oos_resid = oos_resid_ordered,
                 cv_subsamples = cv_subsamples) # Return the folds used
  return(output)
}#CROSSVAL

# Complementary functions ======================================================
crossval_compute <- function(test_sample, learner,
                             y, X, Z = NULL) {

  # Check whether X, Z assignment has been specified. If not, include all.
  if (is.null(learner$assign_X)) learner$assign_X <- seq_len(ncol(X))
  if (!is.null(Z) && is.null(learner$assign_Z)) learner$assign_Z <- seq_len(ncol(Z))

  # Extract model arguments
  mdl_fun <- list(what = learner$fun, args = learner$args)
  assign_X <- learner$assign_X
  assign_Z <- learner$assign_Z # Will be NULL if Z is NULL or assign_Z not set

  # Combine X and Z for training and prediction data subsets
  train_idx <- setdiff(seq_along(y), test_sample)
  X_train <- X[train_idx, assign_X, drop = FALSE]
  X_test  <- X[test_sample, assign_X, drop = FALSE]
  Z_train <- if (!is.null(assign_Z)) Z[train_idx, assign_Z, drop = FALSE] else NULL
  Z_test  <- if (!is.null(assign_Z)) Z[test_sample, assign_Z, drop = FALSE] else NULL

  # Combine X and Z matrices for learner input
  combined_X_train <- cbind(X_train, Z_train)
  combined_X_test  <- cbind(X_test, Z_test)

  # Compute model for this fold
  mdl_fun$args$y <- y[train_idx]
  mdl_fun$args$X <- combined_X_train
  mdl_fit <- tryCatch({
    do.call(mdl_fun$what, mdl_fun$args)
  }, error = function(e) {
    # Append context to the error message before stopping
    stop("Error during base learner execution: ", e$message, call. = FALSE)
  })

  # Compute out of sample residuals
  # Ensure predict method is found (via exported objects/loaded packages)
  oos_fitted <- tryCatch({
    # Use stats::predict as a robust default, allowing S3 dispatch
    stats::predict(mdl_fit, combined_X_test)
  }, error = function(e) {
      # Append context to the error message
      stop("Error during predict method execution: ", e$message, call. = FALSE)
  })

  # Ensure predictions are conformable (vector or single-column matrix)
  if (is.matrix(oos_fitted) && ncol(oos_fitted) > 1) {
    warning("Predict method returned multiple columns; using only the first.")
    oos_fitted <- oos_fitted[, 1, drop = TRUE]
  }
  # Coerce to vector for consistent subtraction
  oos_fitted_vec <- as.vector(oos_fitted)

  # Check dimensions
  if (length(oos_fitted_vec) != length(test_sample)) {
      stop(paste("Prediction output has incorrect length. Expected:",
                 length(test_sample), "Got:", length(oos_fitted_vec)))
  }

  oos_resid <- y[test_sample] - oos_fitted_vec

  # Return residuals for this fold (as a vector)
  return(oos_resid)
}#CROSSVAL_COMPUTE

# --- Helper function for setting up parallel execution ---
#' @noRd
setup_parallel_cluster <- function(num.cores,
                                   vars_to_export, export_envir,
                                   user_export_vars, user_export_envir,
                                   packages_to_load) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' needed for parallel execution. Falling back to sequential.", call. = FALSE)
    return(NULL)
  }

  if (is.null(num.cores)) {
    ncores <- parallel::detectCores()
    if (is.na(ncores) || ncores <= 1) {
        ncores <- 1 # Fallback if detection fails or only 1 core
    } else {
        ncores <- ncores - 1L # Default N-1
    }
  } else {
    ncores <- max(1L, as.integer(num.cores))
  }

  if (ncores <= 1) {
      # Don't bother setting up cluster for 1 core, use sequential path
      # message("Only 1 core available or requested. Running sequentially.") # Optional message
      return(NULL)
  }

  # Create PSOCK cluster
  cl <- tryCatch({
    parallel::makeCluster(ncores, type = "PSOCK")
  }, error = function(e) {
    warning("Failed to create PSOCK cluster: ", e$message,
            "\nFalling back to sequential execution.", call. = FALSE)
    NULL
  })

  if (is.null(cl)) return(NULL)

  # Set RNG stream for reproducibility (uses global seed set *before* crossval call)
  # Seeding the stream ensures independent + reproducible random numbers on workers
  tryCatch({
    parallel::clusterSetRNGStream(cl, NULL) # NULL uses current global .Random.seed
  }, error = function(e) {
    warning("Failed to set RNG stream on cluster: ", e$message, call. = FALSE)
    parallel::stopCluster(cl)
    return(NULL)
  })

  # Export essential internal variables from the crossval environment
  tryCatch({
    parallel::clusterExport(cl, varlist = vars_to_export, envir = export_envir)
  }, error = function(e) {
     warning("Failed to export internal variables to cluster: ", e$message, call. = FALSE)
     parallel::stopCluster(cl)
     return(NULL)
  })

  # Export user-specified objects from the calling environment (parent.frame)
  if (!is.null(user_export_vars) && length(user_export_vars) > 0) {
    # Check if objects exist before attempting export
    missing_vars <- setdiff(user_export_vars, ls(envir = user_export_envir))
    if (length(missing_vars) > 0) {
        warning("Object(s) specified in parallel.export not found in calling environment: ",
                paste(missing_vars, collapse=", "), "\nFalling back to sequential.", call. = FALSE)
        parallel::stopCluster(cl)
        return(NULL)
    }
    tryCatch({
      parallel::clusterExport(cl, varlist = user_export_vars, envir = user_export_envir)
    }, error = function(e) {
      warning("Failed to export objects specified in parallel.export: ", e$message,
              "\nFalling back to sequential.", call. = FALSE)
      parallel::stopCluster(cl)
      return(NULL)
    })
  }

  # Load necessary packages on workers
  # Ensure base packages are included, avoid duplicates
  pkgs_to_load_unique <- unique(c("stats", "methods", packages_to_load))
  if (length(pkgs_to_load_unique) > 0) {
    tryCatch({
      # Use clusterCall to execute library() on each worker
      parallel::clusterCall(cl, function(pkgs) {
          # Invisible suppresses the output of lapply
          invisible(lapply(pkgs, library, character.only = TRUE,
                           warn.conflicts = FALSE, quietly = TRUE))
      }, pkgs_to_load_unique) # Pass packages as argument
    }, error = function(e) {
      warning("Failed to load packages specified in parallel.packages on workers: ",
              e$message, "\nFalling back to sequential.", call. = FALSE)
      parallel::stopCluster(cl)
      return(NULL)
    })
  }

  return(cl) # Return the cluster object if setup was successful
}