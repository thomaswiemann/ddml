# Collection of subsampling functions

# Function to get all crossfitting, crossvalidation, and auxiliary indices
get_all_indx <- function(cluster_variable,
                         sample_folds = 10,
                         cv_folds = 10,
                         D = NULL,
                         by_D = !is.null(D),
                         stratify = by_D,
                         balance_on = "clusters", #observations
                         subsamples = NULL,
                         subsamples_byD = NULL,
                         cv_subsamples = NULL,
                         cv_subsamples_byD = NULL,
                         compute_cv_indices = FALSE,
                         compute_aux_X_indices = FALSE) {

  # Get crossfitting indices
  cf_indx <- get_crossfit_indices(cluster_variable = cluster_variable,
                                  sample_folds = sample_folds,
                                  D = D,
                                  by_D = by_D,
                                  stratify = stratify,
                                  balance_on = balance_on,
                                  subsamples = subsamples,
                                  subsamples_byD = subsamples_byD)
  subsamples <- cf_indx$subsamples
  subsamples_byD <- cf_indx$subsamples_byD

  # Get crossvalidation indices
  if (compute_cv_indices) {
    cv_indx <- get_crossval_indices(subsamples = cf_indx$subsamples,
                                    cluster_variable = cluster_variable,
                                    cv_folds = cv_folds,
                                    D = D,
                                    by_D = by_D,
                                    stratify = stratify,
                                    balance_on = balance_on,
                                    cv_subsamples = cv_subsamples,
                                    cv_subsamples_byD = cv_subsamples_byD)
    cv_subsamples <- cv_indx$cv_subsamples
    cv_subsamples_byD <- cv_indx$cv_subsamples_byD
  }#IF

  # Get auxiliary indices for X matrix
  aux_indx <- NULL
  if (compute_aux_X_indices) {
    aux_indx <- get_auxiliary_indx(subsamples_byD, D)
  }#IF

  # Return output as list
  list(subsamples = subsamples, subsamples_byD = subsamples_byD,
       cv_subsamples = cv_subsamples, cv_subsamples_byD = cv_subsamples_byD,
       aux_indx = aux_indx)
}#GET_ALL_INDX

# Function to get crossvalidation indices given crossfitting subsamples
get_crossval_indices <- function(subsamples,
                                 cluster_variable,
                                 cv_folds = 10,
                                 D = NULL,
                                 by_D = !is.null(D),
                                 stratify = by_D,
                                 balance_on = "clusters", #observations
                                 cv_subsamples = NULL,
                                 cv_subsamples_byD = NULL) {
  # Data parameters
  sample_folds <- length(subsamples)
  if (by_D) {
    D_levels <- sort(unique(D))
    nD_levels <- length(D_levels)
  }#IF

  # Check whether subsamples need to be initialized
  compute_cv <- FALSE
  if (is.null(cv_subsamples) & is.null(cv_subsamples_byD)) {
    compute_cv <- TRUE
    cv_subsamples <- cv_subsamples_byD <- rep(list(NULL), sample_folds)
  } else if (is.null(cv_subsamples) & !is.null(cv_subsamples_byD)) {
    stop(paste0("``cv_subsamples`` must also be set when setting ",
         "``cv_subsamples_byD``."))
  } else if (!is.null(cv_subsamples) & is.null(cv_subsamples_byD) & by_D) {
    stop(paste0("When ``by_D==TRUE``, ``cv_subsamples_byD`` must also be set ",
                "when setting ``cv_subsamples``."))
  }#IFELSE

  # Construct cv subsamples
  if (compute_cv) {
    for (k in 1:sample_folds) {
      # Get subsample observations used for training in sample fold k
      cluster_variable_k <- cluster_variable[-subsamples[[k]]]
      D_k <- D[-subsamples[[k]]]
      # Compute subsamples
      cv_tmp_k <- get_crossfit_indices(cluster_variable_k,
                                       sample_folds = cv_folds,
                                       D = D_k,
                                       by_D = by_D,
                                       stratify = stratify,
                                       balance_on = balance_on,
                                       subsamples = cv_subsamples[[k]],
                                       subsamples_byD = cv_subsamples_byD[[k]])
      # Populate cv objects
      cv_subsamples[[k]] <- cv_tmp_k$subsamples
      cv_subsamples_byD[[k]] <- cv_tmp_k$subsamples_byD
    }#FOR

    # Switch level of the cv_subsamples_byD
    cv_subsamples_byD <- switch_list_levels(cv_subsamples_byD)
  }#IF

  # Return output as list
  list(cv_subsamples = cv_subsamples, cv_subsamples_byD = cv_subsamples_byD)
}#GET_CROSSVAL_INDICES

# Function to get subsampling indices
get_crossfit_indices <- function(cluster_variable,
                                 sample_folds = 10,
                                 D = NULL,
                                 by_D = !is.null(D),
                                 stratify = by_D,
                                 balance_on = "clusters", #observations
                                 subsamples = NULL,
                                 subsamples_byD = NULL) {

  # Check whether subsamples need to be constructed
  compute_cf <- FALSE
  if (is.null(subsamples) & is.null(subsamples_byD)) {
    compute_cf <- TRUE
  } else if (is.null(subsamples) & !is.null(subsamples_byD)) {
    stop("``subsamples`` must also be set when setting ``subsamples_byD``.")
  } else if (!is.null(subsamples) & is.null(subsamples_byD) & by_D) {
    stop(paste0("When ``by_D==TRUE``, ``subsamples_byD`` must also be set when",
                " setting ``subsamples``."))
  }#IF

  # Compute crossfit indices
  if (compute_cf) {
    if (stratify) {
      if (!by_D)
        stop("Stratified sampling only works when ``by_D=TRUE``.")
      cl_folds <- get_cf_indices_stratified(cluster_variable = cluster_variable,
                                            sample_folds = sample_folds, D = D,
                                            balance_on = balance_on)
    } else {
      cl_folds <- get_cf_indices_simple(cluster_variable = cluster_variable,
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
  # Check whether subsamples_byD was constructed
  by_D <- !is.null(subsamples_byD)

  # Check which type of warning should be printed
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
  names(fold_counts) <- paste("Fold", 1:sample_folds)
  fold_counts_byD <- NULL
  if (by_D) {
    D_levels <- sort(unique(D))
    fold_counts_byD <- do.call(rbind, lapply(subsamples_byD, lengths))
    colnames(fold_counts_byD) <- paste("Fold", 1:sample_folds)
    rownames(fold_counts_byD) <- paste0("D=", as.character(D_levels))
  }#IF

  # Throw a warning if the smallest crossfit/val iteration uses less than 100
  training_counts <- sum(fold_counts) - fold_counts
  throw_warning <- FALSE
  if (min(training_counts) < 100) throw_warning <- TRUE
  training_counts_byD <- NULL
  if (by_D) {
    training_counts_byD <- sweep(fold_counts_byD, 1, rowSums(fold_counts_byD),
                                 FUN = function(x, y) y - x)
    if (min(training_counts_byD) < 100) throw_warning <- TRUE
  }#IF
  if (throw_warning & !by_D) {
    warning_text <- paste0("One of the ", type, " subsamples only uses ",
                           min(training_counts), " observations for ",
                           "training. Consider increasing ", fold_arg, "",
                           "if possible.")
  } else if (throw_warning & by_D & !stratify) {
    warning_text <- paste0("One of the ", type, " subsamples only uses ",
                           min(training_counts_byD), " observations for ",
                           "training. Consider setting ``stratify=TRUE`` ",
                           "and/or increasing ", fold_arg, " if possible.")

  } else if (throw_warning & by_D & stratify) {
    warning_text <- paste0("One of the ", type, " subsamples only uses ",
                           min(training_counts_byD), " observations for ",
                           "training. Consider increasing ", fold_arg,
                           " if possible.")
  }#IFELSE
  if (throw_warning) warning(warning_text)

  # Return training counts
  list(training_counts = training_counts,
       training_counts_byD = training_counts_byD)
}#CHECK_SUBSAMPLES

# Stratified crossfit indices construction
get_cf_indices_stratified <- function(cluster_variable, sample_folds, D,
                                      balance_on = "observations") {

  # Data parameters
  nobs <- length(cluster_variable)
  D_levels <- sort(unique(D))
  nD_levels <- length(D_levels)
  is_D <- rep(list(NULL), nD_levels)
  nobs_byD <- rep(0, nD_levels)
  nclusters_byD <- rep(0, nD_levels)
  for (d in 1:nD_levels) {
    is_D[[d]] <- which(D == D_levels[d])
    nobs_byD[[d]] <- length(is_D[[d]])
    nclusters_byD <- length(unique(cluster_variable[is_D[[d]]]))
  }#FOR

  # Error if number of clusters is smaller than the number of sample folds
  if (min(nclusters_byD) < sample_folds)
    stop(paste0("The number of clusters is smaller than the number of sample",
                " folds when D=",
                D_levels[which(min(nclusters_byD) < sample_folds)]))

  # Check how crossfit indices should be created
  cluster <- (length(unique(cluster_variable)) != length(cluster_variable))

  # Create crossfit indices
  if (cluster) {

    #' The below code assigns clusters to the sample folds.
    #'
    #' The procedure is aimed to improve balance (wrt to clusters or
    #'     observations per fold) in settings where D might not be constant
    #'     within a cluster. In this case, clusters contribute to multiple
    #'     by-treatment folds.
    #'
    #' To achieve better balance in the by-treatment subsamples, fold assignment
    #'     is done in iterative  fashion, by first assigning clusters with
    #'     higher treatment counts (e.g., more exposed treatment levels or more
    #'     observations). At every step, cluster imbalance is calculated, and
    #'     clusters are assigned to the cluster that increases imbalance the
    #'     least.
    #'
    #' Once clusters are assigned to sample folds, the subsamples and
    #'     subsamples_byD lists are created by mapping the assigned folds to
    #'     the sample indices in the data.

    # Get the list of unique clusters
    clusters <- unique(cluster_variable)
    nclusters <- length(clusters)

    # Throw a warning when there is a large number of clusters
    if (nclusters > 10000)
      warning(paste0("Stratified subsample construction can take a long time",
                     " when there are many clusters. Check whether",
                     " stratification is necessary if you're short on time."))

    # Map cluster ids to indices in cluster_variable
    cl_indx_list <- split(seq_along(cluster_variable),
                          cluster_variable)
    cl_indx_list_byD <- rep(list(NULL), nD_levels)
    for (d in 1:nD_levels) {
      cluster_variable_d <- cluster_variable[D==D_levels[d]]
      cl_indx_list_byD[[d]] <-
        split(seq_along(cluster_variable_d), cluster_variable_d)
    }#FOR

    # Calculate total counts across all observations or clusters
    if (balance_on == "observations") {
      total_D_counts <- table(D)
    } else if (balance_on == "clusters") {
      # For each cluster, record unique D values they have, then get counts
      cluster_D_values <- lapply(clusters, function(cl) {
        D_values <- unique(D[cl_indx_list[[as.character(cl)]]])
        return(D_values)
      })#LAPPLY
      names(cluster_D_values) <- clusters
      total_D_counts <- table(factor(unlist(cluster_D_values),
                                     levels = D_levels))
    }#IFELSE

    # Calculate D counts per cluster
    if (balance_on == "observations") {
      # Counts of observations per D value per cluster
      cluster_D_counts <- lapply(clusters, function(cl) {
        indices <- cl_indx_list[[as.character(cl)]]
        table(factor(D[indices], levels = D_levels))
      })
    } else if (balance_on == "clusters") {
      # Presence (1) or absence (0) of each D value per cluster
      cluster_D_counts <- lapply(clusters, function(cl) {
        D_values <- unique(D[cl_indx_list[[as.character(cl)]]])
        counts <- as.numeric(D_levels %in% D_values)
        names(counts) <- as.character(D_levels)
        return(counts)
      })
    }#IFELSE
    names(cluster_D_counts) <- clusters

    # Create a data frame with clusters and their D counts
    cluster_info <- data.frame(cluster = clusters,
                               total_observations = lengths(cl_indx_list),
                               stringsAsFactors = FALSE)
    cluster_counts_matrix <- do.call(rbind, cluster_D_counts)
    cluster_info <- cbind(cluster_info, cluster_counts_matrix)
    cluster_info <- cluster_info[order(-cluster_info$total_observations), ]

    # Assign clusters to folds to balance D counts
    fold_counts_byD <- matrix(0, nrow = nD_levels, ncol = sample_folds)
    rownames(fold_counts_byD) <- as.character(D_levels)
    fold_assign <- integer(nobs)
    names(fold_assign) <- as.character(clusters)
    overall_folds <- vector("list", sample_folds)
    D_folds_list <- list()
    for (D_value in D_levels) {
      D_folds_list[[as.character(D_value)]] <- vector("list", sample_folds)
    }#FOR
    for (i in 1:nrow(cluster_info)) {
      cl <- cluster_info$cluster[i]
      # Get the cluster's D counts
      D_counts <- as.numeric(cluster_info[i, -(1:2)])
      names(D_counts) <- colnames(cluster_info)[-(1:2)]

      # Calculate imbalance if the cluster is assigned to each fold
      imbalance <- sapply(1:sample_folds, function(k) {
        new_fold_counts_byD <- fold_counts_byD
        new_fold_counts_byD[, k] <- new_fold_counts_byD[, k] + D_counts
        range_per_D <- apply(new_fold_counts_byD, 1,
                             function(x) max(x) - min(x))
        imbalance_value <- sum(range_per_D)
        return(imbalance_value)
      })#SAPPLY

      # Assign cluster to the fold that results in the least imbalance
      best_fold <- which.min(imbalance)
      fold_assign[as.character(cl)] <- best_fold

      # Update the counts for the assigned fold
      fold_counts_byD[, best_fold] <- fold_counts_byD[, best_fold] + D_counts
    }#FOR

    # Create subsamples from fold assignments
    subsamples <- rep(list(NULL), sample_folds)
    subsamples_byD <- rep(list(NULL), nD_levels)
    for (k in 1:sample_folds) {
      subsamples[[k]] <-
        unlist(unname(cl_indx_list[names(which(fold_assign == k))]))
      for (d in 1:nD_levels) {
        subsamples_byD[[d]][[k]] <-
          unlist(unname(cl_indx_list_byD[[d]][names(which(fold_assign == k))]))
      }#FOR
    }#FOR

  } else {

    #' When observations are not clustered, stratified sampling is much easier.
    #'     The below code first creates random subsamples by treatment level.
    #'     Then, the by-treatment level subsamples are combined to create the
    #'     overall subsamples.
    #'
    #'     Note: This approach is not possible in the clustered setting because
    #'     clusters are of varying lengths (affecting balance) and might be
    #'     associated with multiple treatment levels.

    # First, create subsamples by treatment level
    is_D <- rep(list(NULL), nD_levels)
    nobs_byD <- rep(0, nD_levels)
    subsamples_byD <- rep(list(NULL), nD_levels)
    for (d in 1:nD_levels) {
      is_D[[d]] <- which(D == D_levels[d])
      nobs_byD[[d]] <- length(is_D[[d]])
      subsamples_byD[[d]] <- generate_subsamples(nobs_byD[d], sample_folds)
    }#FOR

    # Then, merge subsamples across treatment levels
    subsamples <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      # Sample folds
      for (d in 1:nD_levels) {
        subsamples[[k]] <- c(subsamples[[k]],
                             (1:nobs)[is_D[[d]]][subsamples_byD[[d]][[k]]])
      }#FOR
      subsamples[[k]] <- sort(subsamples[[k]])
    }#FOR
  }#IFELSE

  # Return subsamples and balancing stats
  list(subsamples = subsamples,
       subsamples_byD = subsamples_byD)
}#GET_CF_INDICES_STRATIFIED

# Simple (non-stratified) crossfit indices construction
get_cf_indices_simple <- function(cluster_variable, sample_folds,
                                  by_D = FALSE, D = NULL) {

  # Check how crossfit indices should be created
  nobs <- length(cluster_variable)
  cluster <- (length(unique(cluster_variable)) != nobs)

  if (cluster) {
    # Create temp cluster variable to efficiently map clusters to sample indices
    tmp_cluster <- as.numeric(factor(cluster_variable))
    cluster_map <- split(seq_along(tmp_cluster), tmp_cluster)
    subsamples_temp <- generate_subsamples(length(unique(tmp_cluster)),
                                           sample_folds)
    subsamples <- lapply(subsamples_temp, function (x) {
      unname(unlist(cluster_map[x]))
    })#LAPPLY
  } else {
    subsamples <- generate_subsamples(nobs, sample_folds)
  }#IFELSE

  # Create subsamples_byD (optional)
  subsamples_byD <- fold_counts_byD <- NULL
  if (by_D) {
    D_levels <- sort(unique(D))
    nD_levels <- length(D_levels)
    subsamples_byD <- rep(list(NULL), nD_levels)
    # Iterate through levels of D and sample_folds, populating subsamples_byD
    for (d in 1:nD_levels) {
      tmp_indx <- rep(NA, nobs) # Populate with NAs to filter D!=d samples
      is_Dd <- which(D==D_levels[d]) # Observations with D==d
      tmp_indx[is_Dd] <- seq_along(is_Dd) # Get sample indices for D==d sample
      subsamples_byD[[d]] <- rep(list(NULL), sample_folds)
      for (k in 1:sample_folds) {
        tmp_indx_k <- tmp_indx[subsamples[[k]]] # Select observations in k fold
        subsamples_byD[[d]][[k]] <- tmp_indx_k[!is.na(tmp_indx_k)]
      }#FOR
    }#FOR
  }#IF

  # Return subsamples and balancing stats
  list(subsamples = subsamples,
       subsamples_byD = subsamples_byD)
}#GET_CF_INDICES_SIMPLE

# Simple function to generate subsamples.
generate_subsamples <- function(nobs, sample_folds) {
  sampleframe <- rep(1:sample_folds, ceiling(nobs/sample_folds))
  sample_groups <- sample(sampleframe, size=nobs, replace=F)
  subsamples <- sapply(1:sample_folds,
                       function(x) {which(sample_groups == x)},
                       simplify = F)
  subsamples
}#GENERATE_SUBSAMPLES

# Simple function to accomodate cv_subsamples_byD construction
switch_list_levels <- function(lst) {
  # Check if the list contains only NULL values
  if (all(sapply(lst, is.null))) {
    return(lst)
  }#IF

  # Determine the number of second-level elements (K)
  K <- length(lst[[1]])

  # Initialize an empty list to store results
  result <- vector("list", K)

  # Loop over each second-level element and extract across all first-level elements
  for (i in seq_len(K)) {
    result[[i]] <- lapply(lst, `[[`, i)
  }

  return(result)
}#SWITCH_LIST_LEVELS

# Function to create indices for auxiliary X
get_auxiliary_indx <- function(subsamples_byD, D) {
  # Data parameters
  nobs <- length(D)
  D_levels <- sort(unique(D))
  nD_levels <- length(D_levels)
  is_D <- rep(list(NULL), nD_levels)
  for (d in 1:nD_levels) is_D[[d]] <- which(D == D_levels[d])
  sample_folds <- length(subsamples_byD[[1]])

  auxiliary_indx_list <- rep(list(NULL), nD_levels)
  for (d in 1:nD_levels) {
    auxiliary_indx_list[[d]] <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      for (h in setdiff((1:nD_levels), d)) {
        auxiliary_indx_list[[d]][[k]] <-
          c(auxiliary_indx_list[[d]][[k]],
            (1:nobs)[is_D[[h]]][subsamples_byD[[h]][[k]]])
      }#FOR
    }#FOR
  }#FOR

  # Return output
  auxiliary_indx_list
}#GET_AUXILIARY_INDX

# Function to get X for corresponding auxiliary subsample, used in ddml_ate, etc
get_auxiliary_X <- function(auxiliary_indx_d, X) {
  # Data parameters
  sample_folds <- length(auxiliary_indx_d)

  # Populate auxiliary X list
  auxiliary_X <- rep(list(NULL), sample_folds)
  for (k in 1:sample_folds) {
    auxiliary_X[[k]] <- X[auxiliary_indx_d[[k]], , drop = FALSE]
  }#FOR

  # return output
  auxiliary_X
}#GET_AUXILIARY_X
