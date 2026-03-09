#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a ddml Object
#'
#' Extracts coefficient estimates, standard errors, test
#' statistics, and p-values from a \code{ddml} estimator in a
#' format compatible with \pkg{modelsummary} and the
#' \pkg{broom} ecosystem.
#'
#' @param x A \code{ddml} object.
#' @param ensemble_idx Integer index of the ensemble type to
#'     report. Defaults to 1 (first ensemble type). Set to
#'     \code{NULL} to return results for all ensemble types.
#' @param conf.int Logical. Include confidence interval
#'     columns? Default \code{FALSE}.
#' @param conf.level Confidence level for intervals.
#'     Default 0.95.
#' @param type Character string specifying the
#'     variance-covariance estimator. One of \code{"HC1"}
#'     (default), \code{"HC0"}, or \code{"HC3"}.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, and \code{ensemble_type}. If
#'     \code{conf.int = TRUE}, also \code{conf.low} and
#'     \code{conf.high}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' tidy(plm_fit)
#' tidy(plm_fit, conf.int = TRUE)
#' }
#'
#' @export
#' @method tidy ddml
tidy.ddml <- function(x, ensemble_idx = 1, conf.int = FALSE,
                      conf.level = 0.95,
                      type = "HC1", ...) {
  type <- match.arg(type, c("HC1", "HC0", "HC3"))

  s <- summary(x, type = type)
  inf <- s$coefficients
  nensb <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(ensemble_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    if (any(ensemble_idx < 1) || any(ensemble_idx > nensb)) {
      stop(sprintf("ensemble_idx must be between 1 and %d", nensb))
    }#IF
    j_seq <- ensemble_idx
  }#IFELSE

  # Pre-allocate rows
  n_rows <- length(j_seq) * p
  term <- rep(dimnames(inf)[[1]], length(j_seq))
  ensemble_type <- rep(dimnames(inf)[[3]][j_seq], each = p)
  
  estimate <- numeric(n_rows)
  std.error <- numeric(n_rows)
  statistic <- numeric(n_rows)
  p.value <- numeric(n_rows)

  idx <- 1
  for (j in j_seq) {
    for (k in seq_len(p)) {
      estimate[idx] <- inf[k, 1, j]
      std.error[idx] <- inf[k, 2, j]
      statistic[idx] <- inf[k, 3, j]
      p.value[idx] <- inf[k, 4, j]
      idx <- idx + 1
    }#FOR
  }#FOR

  res <- data.frame(
    term = term,
    estimate = estimate,
    std.error = std.error,
    statistic = statistic,
    p.value = p.value,
    ensemble_type = ensemble_type,
    stringsAsFactors = FALSE
  )

  if (conf.int) {
    z <- stats::qnorm((1 + conf.level) / 2)
    res$conf.low <- res$estimate - z * res$std.error
    res$conf.high <- res$estimate + z * res$std.error
  }#IF

  res
}#TIDY.DDML

#' Glance at a ddml Object
#'
#' Returns a one-row summary of model-level statistics,
#' compatible with \pkg{modelsummary} and the \pkg{broom}
#' ecosystem.
#'
#' @param x A \code{ddml} object.
#' @param ... Currently unused.
#'
#' @return A one-row \code{data.frame} with columns
#'     \code{nobs}, \code{sample_folds}, \code{shortstack},
#'     \code{ensemble_type}, and \code{model_type}.
#'
#' @examples
#' \donttest{
#' y = AE98[, "worked"]
#' D = AE98[, "morekids"]
#' X = AE98[, c("age","agefst","black","hisp","othrace")]
#' plm_fit = ddml_plm(y, D, X,
#'                     learners = list(what = ols),
#'                     sample_folds = 2, silent = TRUE)
#' glance(plm_fit)
#' }
#'
#' @seealso \code{\link{tidy.ddml}},
#'     \code{\link{summary.ddml}}
#'
#' @export
#' @method glance ddml
glance.ddml <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    sample_folds = x$sample_folds,
    shortstack = if (is.null(x$shortstack)) FALSE else x$shortstack,
    ensemble_type = paste(x$ensemble_type, collapse = ", "),
    model_type = class(x)[1],
    estimator_name = if (is.null(x$estimator_name)) class(x)[1] else x$estimator_name,
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML