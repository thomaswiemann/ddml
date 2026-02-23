.onLoad <- function(libname, pkgname) {
  if (requireNamespace("texreg", quietly = TRUE)) {
    for (cls in c("ddml_plm", "ddml_pliv", "ddml_fpliv",
                  "ddml_ate", "ddml_att", "ddml_late")) {
      methods::setMethod("extract", signature = cls,
                         definition = extract_ddml,
                         where = asNamespace(pkgname))
    }#FOR
  }#IF
}#.ONLOAD

# Internal extract function (works for all ddml classes)
extract_ddml <- function(model, ...) {
  s <- summary(model)
  inf <- s$inf_results

  # Use first ensemble type for texreg
  coef_names <- dimnames(inf)[[1]]
  co <- inf[, 1, 1]
  se <- inf[, 2, 1]
  pval <- inf[, 4, 1]

  gof <- c(s$nobs, s$sample_folds)
  gof.names <- c("Num. obs.", "Crossfit folds")
  gof.decimal <- c(FALSE, FALSE)

  texreg::createTexreg(
    coef.names = coef_names,
    coef = co,
    se = se,
    pvalues = pval,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal
  )
}#EXTRACT_DDML
