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
#'             corresponding to variables in \code{Z} that are passed to the
#'             base learner.}
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
#'     Defaults to \code{NULL}, which uses \code{parallel::detectCores()}.
#' @param parallel.export Optional character vector of names of objects
#'     (functions or data) in the calling environment to be explicitly
#'     exported to the parallel workers. This is necessary for globally defined
#'     helper functions or S3 methods used by custom learners.
#'     Example: `parallel.export = c("my_custom_learner", "predict.my_custom_learner")`.
#' @param parallel.packages Optional character vector of additional package names
#'     to be loaded on the parallel workers (using `library()`). Required if
#'     custom learners or their methods depend on functions from these packages.
#'     Example: `parallel.packages = "gbm"`.
#' @param progress String to print before learner and cv fold progress.
#'
#' @section Parallel Processing using the `parallel` package:
#' When \code{parallel = TRUE}, this function uses the base R `parallel` package
#' via a PSOCK cluster. This requires explicit exporting of necessary objects and
#' loading of packages onto the worker nodes.
#' \itemize{
#'     \item **Custom Objects/Functions:** User-defined learners, helper
#'         functions, or S3 methods needed by learners must be specified by name
#'         in the `parallel.export` argument. Ensure these objects exist in the
#'         environment where `crossval` is called.
#'     \item **Packages:** Packages required by custom learners (e.g., `gbm` if
#'         using `gbm::gbm.fit`) must be specified in the `parallel.packages`
#'         argument. Base packages `stats` and `methods` are loaded automatically.
#' }
#' Reproducibility is handled using `parallel::clusterSetRNGStream`. Output from
#' parallel workers is typically not displayed. If parallel setup fails, execution
#' falls back to sequential processing with a warning.
#'
#' @return \code{crossval} returns a list containing the following components:
#'     \describe{
#'         \item{\code{mspe}}{A vector of MSPE estimates,
#'             each corresponding to a base learners (in chronological order).}
#'         \item{\code{oos_resid}}{A matrix of out-of-sample prediction errors,
#'             each column corresponding to a base learners (in chronological
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
#' cv_res_seq <- crossval(y, X,
#'                        learners = list(list(fun = ols),
#'                                        list(fun = mdl_glmnet),
#'                                        list(fun = mdl_glmnet,
#'                                             args = list(alpha = 0))),
#'                        cv_folds = 4,
#'                        silent = TRUE)
#' cv_res_seq$mspe
#'
#' # --- Example using parallel package ---
#'
#' # Define a custom learner (e.g., requires 'gbm' package and S3 method)
#' library(gbm)
#' mdl_gbm <- function(y, X, ...) {
#'   fit <- gbm::gbm.fit(x = X, y = y, ...)
#'   class(fit) <- c("mdl_gbm_wrapper_class", class(fit))
#'   return(fit)
#' }
#' predict.mdl_gbm <- function(object, newdata, ...) {
#'   class(object) <- class(object)[-1L]
#'   stats::predict(object, newdata = newdata, type = "response",
#'                  n.trees = object$n.trees, ...)
#' }
#'
#' learners_custom <- list(list(fun = ols), # ols from ddml package
#'                         list(fun = mdl_gbm, # Function defined globally
#'                              args = list(n.trees = 100, verbose = FALSE)))
#'
#' # Run in parallel using 2 cores
#' # Export the globally defined functions, load the needed package 'gbm'
#' set.seed(123) # Set seed for reproducibility with clusterSetRNGStream
#' cv_res_par <- crossval(y = y, X = X,
#'                        learners = learners_custom,
#'                        cv_folds = 4,
#'                        parallel = TRUE,
#'                        num.cores = 2,
#'                        parallel.export = c("mdl_gbm", "predict.mdl_gbm"),
#'                        parallel.packages = "gbm",
#'                        silent = TRUE) # Progress reporting disabled
#' print(cv_res_par$mspe)
#'
#' # Compare with sequential (set seed again for identical fold generation)
#' set.seed(123)
#' cv_res_seq <- crossval(y = y, X = X,
#'                        learners = learners_custom,
#'                        cv_folds = 4,
#'                        parallel = FALSE,
#'                        silent = TRUE)
#' print(cv_res_seq$mspe)
#' print(all.equal(cv_res_seq$mspe, cv_res_par$mspe))
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
    # Define essential variables to export
    vars_to_export <- c("y", "X", "learners", "cv_subsamples",
                        "crossval_compute", "cv_folds", "nlearners")
    if (!is.null(Z)) vars_to_export <- c(vars_to_export, "Z")

    # Attempt to set up the cluster
    cl <- setup_parallel_cluster(num.cores = num.cores,
                                 vars_to_export = vars_to_export,
                                 export_envir = environment(),
                                 user_export_vars = parallel.export,
                                 user_export_envir = parent.frame(),
                                 packages_to_load = parallel.packages)

    # If setup successful, register stopCluster on exit
    if (!is.null(cl)) {
      on.exit(parallel::stopCluster(cl), add = TRUE)
    } else {
      # Fallback to sequential if setup failed
      parallel <- FALSE
    }
  }# IF parallel

  # --- Execute Calculation (Chunked by Learner) ---
  total_tasks <- cv_folds * nlearners
  results_list <- vector("list", total_tasks)

  for (j in 1:nlearners) {
    # Define task indices for the current learner j
    task_indices <- ((j - 1) * cv_folds + 1):(j * cv_folds)

    if (!is.null(cl)) { # Parallel execution for this learner
      chunk_res <- tryCatch({
        parallel::parSapply(cl, task_indices, cv_fun)
      }, error = function(e) {
        parallel::stopCluster(cl)
        on.exit(NULL)
        stop(sprintf("Error during parallel execution for learner %d: %s", j, e$message), call. = FALSE)
      })

      # Store results, handling different return types from parSapply
      if (is.list(chunk_res)) {
        # Result is already a list (e.g., simplify=FALSE or inconsistent returns)
        results_list[task_indices] <- chunk_res
      } else if (is.matrix(chunk_res)) {
        # Result is a matrix (each column is a fold's residual vector)
        # Split columns into a list before assignment
        results_list[task_indices] <- lapply(seq_len(ncol(chunk_res)), function(i) chunk_res[, i])
      } else if (is.vector(chunk_res)) {
        # Result is a vector (unlikely if folds > 1 obs, but handle defensively)
        # This assumes the vector elements correspond to folds
        results_list[task_indices] <- as.list(chunk_res)
      } else {
        stop(sprintf("Unexpected result type from parallel worker for learner %d: %s", j, class(chunk_res)[1]))
      }

      # Print parallel progress update
      if (!silent) {
        cat(sprintf("\rProgress: Learner %d / %d completed.         ", j, nlearners))
      }

    } else { # Sequential execution for this learner
      # sapply usually returns matrix/vector, but can return list
      chunk_res <- sapply(task_indices, cv_fun)

      # Store results, handling different return types from sapply
      if (is.list(chunk_res)) {
        results_list[task_indices] <- chunk_res
      } else if (is.matrix(chunk_res)) {
        results_list[task_indices] <- lapply(seq_len(ncol(chunk_res)), function(i) chunk_res[, i])
      } else if (is.vector(chunk_res)) {
        results_list[task_indices] <- as.list(chunk_res)
      } else {
        stop(sprintf("Unexpected result type from sequential execution for learner %d: %s", j, class(chunk_res)[1]))
      }
      # Sequential progress is handled inside cv_fun
    }
  }# FOR loop over learners

  # Final newline after parallel progress reporting is done
  if (!is.null(cl) && !silent) cat("\n")

  # --- Compile Results --- (Results are now in results_list)
  # Ensure all elements are vectors/conformable before concatenating
  # (This should be handled by crossval_compute returning vectors)
  oos_resid_list <- results_list

  oos_resid <- do.call(c, oos_resid_list) # Concatenate all residuals
  expected_length <- tryCatch(nobs_used * nlearners, error = function(e) NA)
  if (is.na(expected_length) || !is.numeric(oos_resid) || length(oos_resid) != expected_length) {
    stop(paste("Internal error: Mismatch in expected number/structure of residuals.",
               "Expected length:", expected_length, "Got length:", length(oos_resid)))
  }
  oos_resid_mat <- matrix(oos_resid, nrow = nobs_used, ncol = nlearners)
  row_order <- order(unlist(cv_subsamples))
  oos_resid_ordered <- oos_resid_mat[row_order, , drop = FALSE]

  # Compute MSPE by learner
  mspe <- colMeans(oos_resid_ordered^2)

  # Organize and return output
  output <- list(mspe = mspe,
                 oos_resid = oos_resid_ordered,
                 cv_subsamples = cv_subsamples)
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
    stop("Error during base learner execution in fold: ", e$message, call. = FALSE)
  })

  # Compute out of sample residuals
  # Ensure predict method is found (via exported objects/loaded packages)
  oos_fitted <- tryCatch({
    stats::predict(mdl_fit, combined_X_test)
  }, error = function(e) {
      stop("Error during predict method execution in fold: ", e$message, call. = FALSE)
  })

  # Ensure predictions are conformable
  oos_fitted_mat <- tryCatch({
      methods::as(oos_fitted, "matrix")
  }, error = function(e){
      stop("Could not coerce predictions to matrix: ", e$message, call. = FALSE)
  })
  # Check dimensions
  if (nrow(oos_fitted_mat) != length(test_sample)) {
      stop("Prediction output has incorrect number of rows.")
  }

  oos_resid <- y[test_sample] - oos_fitted_mat

  # Return residuals for this fold
  return(oos_resid)
}#CROSSVAL_COMPUTE

# --- Helper function for setting up parallel execution ---
#' @noRd
setup_parallel_cluster <- function(num.cores,
                                   vars_to_export, export_envir,
                                   user_export_vars, user_export_envir,
                                   packages_to_load) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("Package 'parallel' needed for parallel execution. Falling back.", call. = FALSE)
    return(NULL)
  }

  if (is.null(num.cores)) {
    ncores <- parallel::detectCores()
    if (is.na(ncores)) ncores <- 1 # Fallback
    ncores <- max(1L, ncores - 1L) # Default N-1
  } else {
    ncores <- max(1L, as.integer(num.cores))
  }

  if (ncores <= 1) return(NULL) # Not worth setting up parallel for 1 core

  # Create PSOCK cluster
  cl <- tryCatch({
    parallel::makeCluster(ncores, type = "PSOCK")
  }, error = function(e) {
    warning("Failed to create PSOCK cluster: ", e$message,
            "\nFalling back to sequential execution.", call. = FALSE)
    NULL
  })

  if (is.null(cl)) return(NULL)

  # Set RNG stream for reproducibility (uses global seed)
  tryCatch({
    parallel::clusterSetRNGStream(cl, NULL)
  }, error = function(e) {
    warning("Failed to set RNG stream on cluster: ", e$message, call. = FALSE)
    parallel::stopCluster(cl)
    return(NULL)
  })

  # Export essential internal variables
  tryCatch({
    parallel::clusterExport(cl, varlist = vars_to_export, envir = export_envir)
  }, error = function(e) {
     warning("Failed to export internal variables to cluster: ", e$message, call. = FALSE)
     parallel::stopCluster(cl)
     return(NULL)
  })

  # Export user-specified objects
  if (!is.null(user_export_vars) && length(user_export_vars) > 0) {
    tryCatch({
      parallel::clusterExport(cl, varlist = user_export_vars, envir = user_export_envir)
    }, error = function(e) {
      warning("Failed to export objects in parallel.export: ", e$message,
              "\nEnsure they exist in the calling environment.", call. = FALSE)
      parallel::stopCluster(cl)
      return(NULL)
    })
  }

  # Load necessary packages on workers (more concise loading)
  pkgs <- unique(c("stats", "methods", packages_to_load)) # Base packages always needed
  if (length(pkgs) > 0) {
    tryCatch({
      parallel::clusterCall(cl, function(pkgs_to_load) {
          invisible(lapply(pkgs_to_load, library, character.only = TRUE))
      }, pkgs) # Pass pkgs as argument to the function run on nodes
    }, error = function(e) {
      warning("Failed to load packages in parallel.packages on workers: ",
              e$message, call. = FALSE)
      parallel::stopCluster(cl)
      return(NULL)
    })
  }

  return(cl)
}
