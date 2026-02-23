#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' Tidy a ddml object
#'
#' @param x A `ddml` object.
#' @param conf.int Logical indicating whether or not to include a confidence interval.
#' @param conf.level The confidence level to use.
#' @param ... Additional arguments passed to summary method.
#'
#' @return A data.frame containing tidy coefficients per ensemble type.
#' @export
tidy.ddml <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary(x)
  inf <- s$inf_results
  nensb <- dim(inf)[3]
  p <- dim(inf)[1]

  rows <- list()
  for (j in seq_len(nensb)) {
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

#' Glance at a ddml object
#'
#' @param x A `ddml` object.
#' @param ... Additional arguments.
#'
#' @return A one-row data.frame with model metadata.
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
