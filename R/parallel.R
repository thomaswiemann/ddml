# Internal helper for parallel cluster setup.
setup_parallel_cluster <- function(num_cores,
                                   parallel_export = NULL,
                                   parallel_packages = NULL) {
  cl <- parallel::makeCluster(num_cores, type = "PSOCK")
  parallel::clusterEvalQ(cl, library(ddml))
  if (!is.null(parallel_export)) {
    parallel::clusterExport(cl, varlist = parallel_export,
                            envir = globalenv())
  }#IF
  if (!is.null(parallel_packages)) {
    for (pkg in parallel_packages) {
      parallel::clusterCall(cl, library, pkg,
                            character.only = TRUE)
    }#FOR
  }#IF
  
  # Ensure reproducible parallel random number generation (L'Ecuyer-CMRG)
  parallel::clusterSetRNGStream(cl)
  
  cl
}#SETUP_PARALLEL_CLUSTER

# Internal helper to unpack the parallel list argument.
parse_parallel <- function(parallel) {
  if (is.null(parallel)) {
    return(list(num_cores = 1, export = NULL,
                packages = NULL))
  }#IF
  if (!is.list(parallel)) {
    stop("'parallel' must be a list or NULL.",
         call. = FALSE)
  }#IF
  list(
    num_cores = if (!is.null(parallel$cores)) {
      parallel$cores
    } else {
      1
    },
    export = parallel$export,
    packages = parallel$packages
  )
}#PARSE_PARALLEL

# Run a function over seq_len(njobs) with optional parallel cluster.
# Encapsulates parsing, cluster lifecycle, fallback, and pbapply.
with_parallel <- function(njobs, fun, parallel, silent) {
  p <- parse_parallel(parallel)
  cl <- NULL
  if (p$num_cores > 1) {
    cl <- tryCatch(
      setup_parallel_cluster(p$num_cores, p$export,
                             p$packages),
      error = function(e) {
        warning("Parallel setup failed: ",
                conditionMessage(e),
                ". Falling back to sequential.",
                call. = FALSE)
        NULL
      })
  }#IF
  on.exit({
    if (!is.null(cl)) parallel::stopCluster(cl)
  }, add = TRUE)
  if (silent) {
    op <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(op), add = TRUE)
  }#IF
  pbapply::pblapply(seq_len(njobs), fun, cl = cl)
}#WITH_PARALLEL

# Internal helper for silent-aware messages.
info_msg <- function(..., silent = FALSE) {
  if (!silent) message(...)
}#INFO_MSG
