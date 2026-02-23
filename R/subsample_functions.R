# Collection of subsampling functions

# Main entry point for constructing all sample-splitting indices.
get_sample_splits <- function(cluster_variable,
                              sample_folds = 10,
                              cv_folds = NULL,
                              D = NULL,
                              stratify = !is.null(D),
                              subsamples = NULL,
                              subsamples_byD = NULL,
                              cv_subsamples = NULL,
                              cv_subsamples_byD = NULL) {
  by_D <- !is.null(D)

  # Auto-merge subsamples_byD into subsamples when only byD is given
  if (is.null(subsamples) & !is.null(subsamples_byD)) {
    if (!by_D) stop("subsamples_byD requires D to be specified.")
    subsamples <- merge_subsamples_byD(subsamples_byD, D)
  }#IF

  # Get crossfitting indices
  cf_indx <- get_crossfit_indices(
    cluster_variable = cluster_variable,
    sample_folds = sample_folds,
    D = D, by_D = by_D,
    stratify = stratify,
    subsamples = subsamples,
    subsamples_byD = subsamples_byD)
  subsamples <- cf_indx$subsamples
  subsamples_byD <- cf_indx$subsamples_byD

  # Build CV indices if requested
  if (!is.null(cv_folds) && is.null(cv_subsamples)) {
    n_sf <- length(subsamples)
    cv_subsamples <- rep(list(NULL), n_sf)
    cv_byD_raw <- if (by_D) rep(list(NULL), n_sf)

    for (k in seq_len(n_sf)) {
      cl_k <- cluster_variable[-subsamples[[k]]]
      if (by_D) {
        D_k <- D[-subsamples[[k]]]
        cv_tmp <- get_crossfit_indices(
          cl_k, sample_folds = cv_folds,
          D = D_k, by_D = TRUE, stratify = FALSE)
        cv_subsamples[[k]] <- cv_tmp$subsamples
        cv_byD_raw[[k]] <- cv_tmp$subsamples_byD
      } else {
        cv_tmp <- get_crossfit_indices(
          cl_k, sample_folds = cv_folds)
        cv_subsamples[[k]] <- cv_tmp$subsamples
      }#IFELSE
    }#FOR

    if (by_D && is.null(cv_subsamples_byD)) {
      cv_subsamples_byD <- switch_list_levels(cv_byD_raw)
    }#IF
  }#IF

  # Build auxiliary indices when D is given
  aux_indx <- NULL
  if (by_D) {
    aux_indx <- get_auxiliary_indx(subsamples_byD, D)
  }#IF

  # Clean up: return NULL (not list-of-NULLs) when not applicable
  if (!by_D) {
    subsamples_byD <- NULL
    cv_subsamples_byD <- NULL
  }#IF

  # Return output as list
  list(subsamples = subsamples,
       subsamples_byD = subsamples_byD,
       cv_subsamples = cv_subsamples,
       cv_subsamples_byD = cv_subsamples_byD,
       aux_indx = aux_indx)
}#GET_SAMPLE_SPLITS

# Internal: merge by-D subsamples into full-sample subsamples
merge_subsamples_byD <- function(subsamples_byD, D) {
  D_levels <- sort(unique(D))
  nD_levels <- length(D_levels)
  nobs <- length(D)
  is_D <- lapply(seq_len(nD_levels),
                 function(d) which(D == D_levels[d]))
  sample_folds <- length(subsamples_byD[[1]])

  subsamples <- rep(list(NULL), sample_folds)
  for (k in seq_len(sample_folds)) {
    for (d in seq_len(nD_levels)) {
      subsamples[[k]] <- c(subsamples[[k]],
                           is_D[[d]][subsamples_byD[[d]][[k]]])
    }#FOR
    subsamples[[k]] <- sort(subsamples[[k]])
  }#FOR
  subsamples
}#MERGE_SUBSAMPLES_BYD

# Internal: dispatcher for crossfit index construction
get_crossfit_indices <- function(cluster_variable,
                                 sample_folds = 10,
                                 D = NULL,
                                 by_D = !is.null(D),
                                 stratify = by_D,
                                 subsamples = NULL,
                                 subsamples_byD = NULL) {

  # Check whether subsamples need to be constructed
  compute_cf <- FALSE
  if (is.null(subsamples) & is.null(subsamples_byD)) {
    compute_cf <- TRUE
  } else if (is.null(subsamples) & !is.null(subsamples_byD)) {
    stop(paste0("``subsamples`` must also be set when setting ",
                "``subsamples_byD``."))
  } else if (!is.null(subsamples) & is.null(subsamples_byD) & by_D) {
    stop(paste0("When ``by_D==TRUE``, ``subsamples_byD`` must also be set ",
                "when setting ``subsamples``."))
  }#IF

  # Compute crossfit indices
  if (compute_cf) {
    if (stratify) {
      if (!by_D)
        stop("Stratified sampling only works when ``by_D=TRUE``.")
      cl_folds <- get_cf_indices_stratified(
        cluster_variable = cluster_variable,
        sample_folds = sample_folds, D = D)
    } else {
      cl_folds <- get_cf_indices_simple(
        cluster_variable = cluster_variable,
        sample_folds = sample_folds,
        by_D = by_D, D = D)
    }#IFELSE
    subsamples <- cl_folds$subsamples
    subsamples_byD <- cl_folds$subsamples_byD
  }#IF

  # Return list of NULLs if by_D is FALSE
  if (!by_D) subsamples_byD <- rep(list(NULL), sample_folds)

  # Return output as list
  list(subsamples = subsamples, subsamples_byD = subsamples_byD)
}#GET_CROSSFIT_INDICES

# Checks size of crossfitting subsamples
check_subsamples <- function(subsamples, subsamples_byD, stratify,
                             D = NULL, cv = FALSE) {
  by_D <- !is.null(subsamples_byD)

  if (!cv) {
    type <- "crossfitting"
    fold_arg <- "``sample_folds``"
  } else {
    type <- "crossvalidation"
    fold_arg <- "``sample_folds`` and/or ``cv_folds``"
  }#IFELSE

  # Compute fold counts to check for balance
  sample_folds <- length(subsamples)
  fold_counts <- lengths(subsamples)
  names(fold_counts) <- paste("Fold", seq_len(sample_folds))
  fold_counts_byD <- NULL
  if (by_D) {
    D_levels <- sort(unique(D))
    fold_counts_byD <- do.call(rbind, lapply(subsamples_byD, lengths))
    colnames(fold_counts_byD) <- paste("Fold", seq_len(sample_folds))
    rownames(fold_counts_byD) <- paste0("D=", as.character(D_levels))
  }#IF

  # Throw a warning if the smallest training set uses less than 100 obs
  training_counts <- sum(fold_counts) - fold_counts
  throw_warning <- FALSE
  if (min(training_counts) < 100) throw_warning <- TRUE
  training_counts_byD <- NULL
  if (by_D) {
    training_counts_byD <- sweep(fold_counts_byD, 1,
                                 rowSums(fold_counts_byD),
                                 FUN = function(x, y) y - x)
    if (min(training_counts_byD) < 100) throw_warning <- TRUE
  }#IF
  if (throw_warning & !by_D) {
    warning_text <- paste0(
      "One of the ", type, " subsamples only uses ",
      min(training_counts), " observations for training. Consider ",
      "increasing ", fold_arg, " if possible.")
  } else if (throw_warning & by_D & !stratify) {
    warning_text <- paste0(
      "One of the ", type, " subsamples only uses ",
      min(training_counts_byD), " observations for training. Consider ",
      "setting ``stratify=TRUE`` and/or increasing ", fold_arg,
      " if possible.")
  } else if (throw_warning & by_D & stratify) {
    warning_text <- paste0(
      "One of the ", type, " subsamples only uses ",
      min(training_counts_byD), " observations for training. Consider ",
      "increasing ", fold_arg, " if possible.")
  }#IFELSE
  if (throw_warning) warning(warning_text)

  # Return training counts (invisible)
  invisible(list(training_counts = training_counts,
                 training_counts_byD = training_counts_byD))
}#CHECK_SUBSAMPLES

# Stratified crossfit indices construction
get_cf_indices_stratified <- function(cluster_variable, sample_folds,
                                      D) {

  # Data parameters
  nobs <- length(cluster_variable)
  D_levels <- sort(unique(D))
  nD_levels <- length(D_levels)
  is_D <- rep(list(NULL), nD_levels)
  nobs_byD <- rep(0, nD_levels)
  nclusters_byD <- rep(0, nD_levels)
  for (d in seq_len(nD_levels)) {
    is_D[[d]] <- which(D == D_levels[d])
    nobs_byD[[d]] <- length(is_D[[d]])
    nclusters_byD[d] <- length(unique(cluster_variable[is_D[[d]]]))
  }#FOR

  # Error if number of clusters is smaller than the number of sample folds
  if (min(nclusters_byD) < sample_folds)
    stop(paste0("Number of clusters for at least one treatment level is ",
                "smaller than the number of sample folds."))

  # Check for clustering
  cluster <- (length(unique(cluster_variable)) != nobs)

  if (cluster) {

    clusters <- unique(cluster_variable)
    nclusters <- length(clusters)

    if (nclusters > 10000)
      warning(paste0(
        "Stratified subsample construction can take a long time",
        " when there are many clusters. Check whether",
        " stratification is necessary if you're short on time."))

    # Map cluster ids to indices in cluster_variable
    cl_indx_list <- split(seq_along(cluster_variable),
                          cluster_variable)
    cl_indx_list_byD <- rep(list(NULL), nD_levels)
    for (d in seq_len(nD_levels)) {
      cluster_variable_d <- cluster_variable[D == D_levels[d]]
      cl_indx_list_byD[[d]] <-
        split(seq_along(cluster_variable_d), cluster_variable_d)
    }#FOR

    # Calculate total D counts across clusters
    cluster_D_values <- lapply(clusters, function(cl) {
      unique(D[cl_indx_list[[as.character(cl)]]])
    })#LAPPLY
    names(cluster_D_values) <- clusters
    total_D_counts <- table(factor(unlist(cluster_D_values),
                                   levels = D_levels))

    # Calculate D counts per cluster (binary: is D level present?)
    cluster_D_counts <- lapply(clusters, function(cl) {
      D_values <- unique(D[cl_indx_list[[as.character(cl)]]])
      counts <- as.numeric(D_levels %in% D_values)
      names(counts) <- as.character(D_levels)
      counts
    })
    names(cluster_D_counts) <- clusters

    # Create a data frame with clusters and their D counts
    cluster_info <- data.frame(
      cluster = clusters,
      total_observations = lengths(cl_indx_list)[as.character(clusters)],
      stringsAsFactors = FALSE)
    cluster_counts_matrix <- do.call(rbind, cluster_D_counts)
    cluster_info <- cbind(cluster_info, cluster_counts_matrix)
    cluster_info <- cluster_info[
      order(-cluster_info$total_observations), ]

    # Assign clusters to folds to balance D counts
    fold_counts_byD <- matrix(0, nrow = nD_levels,
                              ncol = sample_folds)
    rownames(fold_counts_byD) <- as.character(D_levels)
    fold_assign <- integer(nclusters)
    names(fold_assign) <- as.character(clusters)
    D_folds_list <- list()
    for (D_value in D_levels) {
      D_folds_list[[as.character(D_value)]] <-
        vector("list", sample_folds)
    }#FOR
    for (i in seq_len(nrow(cluster_info))) {
      cl <- cluster_info$cluster[i]
      D_counts <- as.numeric(cluster_info[i, -(1:2)])
      names(D_counts) <- colnames(cluster_info)[-(1:2)]

      imbalance <- vapply(seq_len(sample_folds), function(k) {
        new_fold_counts_byD <- fold_counts_byD
        new_fold_counts_byD[, k] <-
          new_fold_counts_byD[, k] + D_counts
        range_per_D <- apply(new_fold_counts_byD, 1,
                             function(x) max(x) - min(x))
        sum(range_per_D)
      }, FUN.VALUE = numeric(1))#VAPPLY

      best_fold <- which.min(imbalance)
      fold_assign[as.character(cl)] <- best_fold

      fold_counts_byD[, best_fold] <-
        fold_counts_byD[, best_fold] + D_counts
    }#FOR

    # Create subsamples from fold assignments
    subsamples <- rep(list(NULL), sample_folds)
    subsamples_byD <- rep(list(NULL), nD_levels)
    for (k in seq_len(sample_folds)) {
      subsamples[[k]] <-
        unlist(unname(
          cl_indx_list[names(which(fold_assign == k))]))
      for (d in seq_len(nD_levels)) {
        subsamples_byD[[d]][[k]] <-
          unlist(unname(
            cl_indx_list_byD[[d]][names(which(fold_assign == k))]))
      }#FOR
    }#FOR

  } else {

    # Non-clustered: create random subsamples by treatment level
    is_D <- rep(list(NULL), nD_levels)
    nobs_byD <- rep(0, nD_levels)
    subsamples_byD <- rep(list(NULL), nD_levels)
    for (d in seq_len(nD_levels)) {
      is_D[[d]] <- which(D == D_levels[d])
      nobs_byD[[d]] <- length(is_D[[d]])
      subsamples_byD[[d]] <- generate_subsamples(nobs_byD[d],
                                                 sample_folds)
    }#FOR

    # Merge subsamples across treatment levels
    subsamples <- rep(list(NULL), sample_folds)
    for (k in seq_len(sample_folds)) {
      for (d in seq_len(nD_levels)) {
        subsamples[[k]] <- c(subsamples[[k]],
          seq_len(nobs)[is_D[[d]]][subsamples_byD[[d]][[k]]])
      }#FOR
      subsamples[[k]] <- sort(subsamples[[k]])
    }#FOR
  }#IFELSE

  list(subsamples = subsamples,
       subsamples_byD = subsamples_byD)
}#GET_CF_INDICES_STRATIFIED

# Simple (non-stratified) crossfit indices construction
get_cf_indices_simple <- function(cluster_variable, sample_folds,
                                  by_D = FALSE, D = NULL) {

  nobs <- length(cluster_variable)
  cluster <- (length(unique(cluster_variable)) != nobs)

  if (cluster) {
    tmp_cluster <- as.numeric(factor(cluster_variable))
    cluster_map <- split(seq_along(tmp_cluster), tmp_cluster)
    subsamples_temp <- generate_subsamples(
      length(unique(tmp_cluster)), sample_folds)
    subsamples <- lapply(subsamples_temp, function(x) {
      unname(unlist(cluster_map[x]))
    })#LAPPLY
  } else {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IFELSE

  # Create subsamples_byD (optional)
  subsamples_byD <- NULL
  if (by_D) {
    D_levels <- sort(unique(D))
    nD_levels <- length(D_levels)
    subsamples_byD <- rep(list(NULL), nD_levels)
    for (d in seq_len(nD_levels)) {
      tmp_indx <- rep(NA, nobs)
      is_Dd <- which(D == D_levels[d])
      tmp_indx[is_Dd] <- seq_along(is_Dd)
      subsamples_byD[[d]] <- rep(list(NULL), sample_folds)
      for (k in seq_len(sample_folds)) {
        tmp_indx_k <- tmp_indx[subsamples[[k]]]
        subsamples_byD[[d]][[k]] <- tmp_indx_k[!is.na(tmp_indx_k)]
      }#FOR
    }#FOR
  }#IF

  list(subsamples = subsamples,
       subsamples_byD = subsamples_byD)
}#GET_CF_INDICES_SIMPLE

# Simple function to generate subsamples.
generate_subsamples <- function(nobs, sample_folds) {
  sampleframe <- rep(seq_len(sample_folds),
                     ceiling(nobs / sample_folds))
  sample_groups <- sample(sampleframe, size = nobs, replace = FALSE)
  vapply(seq_len(sample_folds),
         function(x) list(which(sample_groups == x)),
         FUN.VALUE = list(1))
}#GENERATE_SUBSAMPLES

# Utility to swap nesting levels of a nested list
switch_list_levels <- function(lst) {
  if (all(vapply(lst, is.null, FUN.VALUE = logical(1)))) {
    return(lst)
  }#IF
  K <- length(lst[[1]])
  result <- vector("list", K)
  for (i in seq_len(K)) {
    result[[i]] <- lapply(lst, `[[`, i)
  }#FOR
  return(result)
}#SWITCH_LIST_LEVELS

# Function to create indices for auxiliary X
get_auxiliary_indx <- function(subsamples_byD, D) {
  nobs <- length(D)
  D_levels <- sort(unique(D))
  nD_levels <- length(D_levels)
  is_D <- rep(list(NULL), nD_levels)
  for (d in seq_len(nD_levels)) is_D[[d]] <- which(D == D_levels[d])
  sample_folds <- length(subsamples_byD[[1]])

  auxiliary_indx_list <- rep(list(NULL), nD_levels)
  for (d in seq_len(nD_levels)) {
    auxiliary_indx_list[[d]] <- rep(list(NULL), sample_folds)
    for (k in seq_len(sample_folds)) {
      for (h in setdiff(seq_len(nD_levels), d)) {
        auxiliary_indx_list[[d]][[k]] <-
          c(auxiliary_indx_list[[d]][[k]],
            seq_len(nobs)[is_D[[h]]][subsamples_byD[[h]][[k]]])
      }#FOR
    }#FOR
  }#FOR
  auxiliary_indx_list
}#GET_AUXILIARY_INDX

# Function to get X for corresponding auxiliary subsample
get_auxiliary_X <- function(auxiliary_indx_d, X) {
  sample_folds <- length(auxiliary_indx_d)
  auxiliary_X <- rep(list(NULL), sample_folds)
  for (k in seq_len(sample_folds)) {
    auxiliary_X[[k]] <- X[auxiliary_indx_d[[k]], , drop = FALSE]
  }#FOR
  auxiliary_X
}#GET_AUXILIARY_X
