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
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{term},
#'     \code{estimate}, \code{std.error}, \code{statistic},
#'     \code{p.value}, and \code{ensemble_type}. If
#'     \code{conf.int = TRUE}, also \code{conf.low} and
#'     \code{conf.high}.
#'
#' @family ddml
#' @export
tidy.ddml <- function(x, ensemble_idx = 1, conf.int = FALSE,
                      conf.level = 0.95, ...) {
  s <- summary(x)
  inf <- s$inf_results
  nensb <- dim(inf)[3]
  p <- dim(inf)[1]

  if (is.null(ensemble_idx)) {
    j_seq <- seq_len(nensb)
  } else {
    j_seq <- ensemble_idx
  }#IFELSE

  rows <- list()
  for (j in j_seq) {
    for (k in seq_len(p)) {
      row <- data.frame(
        term = dimnames(inf)[[1]][k],
        estimate = inf[k, 1, j],
        std.error = inf[k, 2, j],
        statistic = inf[k, 3, j],
        p.value = inf[k, 4, j],
        ensemble_type = dimnames(inf)[[3]][j],
        stringsAsFactors = FALSE
      )
      if (conf.int) {
        z <- stats::qnorm((1 + conf.level) / 2)
        row$conf.low <- inf[k, 1, j] - z * inf[k, 2, j]
        row$conf.high <- inf[k, 1, j] + z * inf[k, 2, j]
      }#IF
      rows[[length(rows) + 1]] <- row
    }#FOR
  }#FOR
  do.call(rbind, rows)
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
#' @family ddml
#' @export
glance.ddml <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    sample_folds = x$sample_folds,
    shortstack = if (is.null(x$shortstack)) FALSE else x$shortstack,
    ensemble_type = paste(x$ensemble_type, collapse = ", "),
    model_type = class(x)[1],
    stringsAsFactors = FALSE
  )
}#GLANCE.DDML
