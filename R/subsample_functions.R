# Collection of subsampling functions
get_crossfit_indices <- function(cluster_variable,
                                 sample_folds = 10, cv_folds = 10,
                                 D = NULL,
                                 subsamples = NULL,
                                 cv_subsamples_list = NULL,
                                 subsamples_byD = NULL,
                                 cv_subsamples_byD = NULL) {

  # Data parameters
  nobs <- length(cluster_variable)
  cluster <- !identical(cluster_variable,  1:nobs)

  # Check whether indices should be constructed by D
  by_D <- !is.null(D)

  if (by_D) {

    # Argument check
    if (!is.null(subsamples) & is.null(subsamples_byD))
      stop("Must also supply subsamples_byD if supplying subsamples.")
    if (!is.null(cv_subsamples_byD) & is.null(subsamples_byD))
      stop("Must also supply subsamples_byD if supplying cv_subsamples_byD.")
    if (!is.null(cv_subsamples_list) & is.null(cv_subsamples_byD))
      stop("Must also supply cv_subsamples_byD if supplying cv_subsamples_list")

    # Data parameters
    D_levels <- sort(unique(D))
    n_D_levels <- length(D_levels)
    is_D <- rep(list(NULL), n_D_levels)
    nobs_byD <- rep(0, n_D_levels)
    for (d in 1:n_D_levels) {
      is_D[[d]] <- which(D == D_levels[d])
      nobs_byD[[d]] <- length(is_D[[d]])
    }#FOR

    if (is.null(subsamples_byD)) {
      # Create sample fold tuple by treatment levels
      subsamples_byD <- rep(list(NULL), n_D_levels)
      for (d in 1:n_D_levels) {
        if (cluster) {
          # Create temp cluster variable to efficiently map clusters to sample indices
          tmp_cl <- get_temp_cluster(cluster_variable[D==D_levels[d]])
          subsamples_temp <- generate_subsamples(tmp_cl$n_cluster,
                                                 sample_folds)
          subsamples_byD[[d]] <- lapply(subsamples_temp, function (x) {
            unname(unlist(tmp_cl$cluster_map[x]))
          })#LAPPLY
        } else {
          subsamples_byD[[d]] <- generate_subsamples(nobs_byD[d], sample_folds)
        }#IFELSE
      }#FOR
    }#IF
    sample_folds <- length(subsamples_byD[[1]])




    # Merge subsamples across treatment levels
    subsamples <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      # Sample folds
      for (d in 1:n_D_levels) {
        subsamples[[k]] <- c(subsamples[[k]],
                             (1:nobs)[is_D[[d]]][subsamples_byD[[d]][[k]]])
      }#FOR
      subsamples[[k]] <- sort(subsamples[[k]])
    }#FOR

    # Warning: CV subsample creation currently ignores dependence!

    # Create CV subsamples by treatment level
    if (is.null(cv_subsamples_byD)) {
      cv_subsamples_byD <- rep(list(NULL), n_D_levels)
      for (d in 1:n_D_levels) {
        cv_subsamples_byD[[d]] <- rep(list(NULL), sample_folds)
        for (k in 1:sample_folds) {
          nobs_d_k <- nobs_byD[[d]] - length(subsamples_byD[[d]][[k]])
          cv_subsamples_byD[[d]][[k]] <-
            generate_subsamples(nobs_d_k, cv_folds)
        }# FOR
      }#FOR
    }#IF
    cv_folds <- length(cv_subsamples_byD[[1]][[1]])

    # Merge cv_subsamples across treatment levels
    cv_subsamples_list <- rep(list(NULL), sample_folds)
    for (k in 1:sample_folds) {
      # CV folds
      cv_subsamples_list[[k]] <- rep(list(NULL), cv_folds)
      for (d in 1:n_D_levels) {
        is_d_k <- which(D[-subsamples[[k]]] == D_levels[d])
        for (j in 1:cv_folds) {
          cv_subsamples_list[[k]][[j]] <-
            c(cv_subsamples_list[[k]][[j]],
              is_d_k[cv_subsamples_byD[[d]][[k]][[j]]])
        }#FOR
      }#FOR
      for (d in 1:n_D_levels) {
        cv_subsamples_list[[k]][[j]] <- sort(cv_subsamples_list[[k]][[j]])
      }#FOR
    }#FOR

  } else {

    # Argument check
    if (!is.null(cv_subsamples_list) & is.null(subsamples))
      stop("Must also supply subsamples. if supplying cv_subsamples_list")

    # Create sample fold tuple
    if (is.null(subsamples)) {
      if (cluster) {
        # Create temp cluster variable to efficiently map clusters to sample indices
        tmp_cl <- get_temp_cluster(cluster_variable)
        subsamples_temp <- generate_subsamples(tmp_cl$n_cluster,
                                               sample_folds)
        subsamples <- lapply(subsamples_temp, function (x) {
          unname(unlist(tmp_cl$cluster_map[x]))
        })#LAPPLY
      } else {
        subsamples <- generate_subsamples(nobs, sample_folds)
      }#IFELSE
    }#IF
    sample_folds <- length(subsamples)

    # Warning: CV subsample creation currently ignores dependence!

    # Create cv-subsamples tuple
    if (is.null(cv_subsamples_list)) {
      cv_subsamples_list <- rep(list(NULL), sample_folds)
      for (k in 1:sample_folds) {
        nobs_k <- nobs - length(subsamples[[k]])
        cv_subsamples_list[[k]] <- generate_subsamples(nobs_k, cv_folds)
      }# FOR
    }#IF
    cv_folds <- length(cv_subsamples_list[[1]])

    # Create NULL values for unconstructed objects
    cv_subsamples_byD <- subsamples_byD <- NULL
  }#IFELSE

  # Return output as list
  output <- list(subsamples = subsamples,
                 cv_subsamples_list = cv_subsamples_list,
                 subsamples_byD = subsamples_byD,
                 cv_subsamples_byD = cv_subsamples_byD,
                 sample_folds = sample_folds,
                 cv_folds = cv_folds)
}#GET_CROSSFIT_INDICES

# Function to create indices for auxiliary X
get_auxiliary_indx <- function(subsamples_byD, D) {
    # Data parameters
    nobs <- length(D)
    D_levels <- sort(unique(D))
    n_D_levels <- length(D_levels)
    is_D <- rep(list(NULL), n_D_levels)
    for (d in 1:n_D_levels) is_D[[d]] <- which(D == D_levels[d])
    sample_folds <- length(subsamples_byD[[1]])

    #auxiliary_X_list <- rep(list(NULL), n_D_levels)
    auxiliary_indx_list <- rep(list(NULL), n_D_levels)
    for (d in 1:n_D_levels) {
      auxiliary_indx_list[[d]] <- rep(list(NULL), sample_folds)
      for (k in 1:sample_folds) {
        for (h in setdiff((1:n_D_levels), d)) {
          auxiliary_indx_list[[d]][[k]] <-
            c(auxiliary_indx_list[[d]][[k]],
              (1:nobs)[is_D[[h]]][subsamples_byD[[h]][[k]]])
        }#FOR
      }#FOR
    }#FOR

  # Return output
  auxiliary_indx_list
}#GET_AUXILIARY_INDX

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

# Simple function to generate subsamples.
generate_subsamples <- function(nobs, sample_folds) {
  sampleframe <- rep(1:sample_folds, ceiling(nobs/sample_folds))
  sample_groups <- sample(sampleframe, size=nobs, replace=F)
  subsamples <- sapply(1:sample_folds,
                       function(x) {which(sample_groups == x)},
                       simplify = F)
  subsamples
}#GENERATE_SUBSAMPLES

# Simple function to create a temp cluster variable for more efficient mapping
get_temp_cluster <- function(cluster_variable) {
  tmp_cluster <- as.numeric(factor(cluster_variable))
  cluster_map <- split(seq_along(tmp_cluster), tmp_cluster)
  list(tmp_cluster = tmp_cluster, cluster_map = cluster_map,
       n_cluster = length(unique(tmp_cluster)))
}#GET_TEMP_CLUSTER
