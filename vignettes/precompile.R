#!/usr/bin/env Rscript
# Precompile vignettes that depend on external packages, large data,
# or long-running computations.
#
# Usage (from package root):
#   Rscript vignettes/precompile.R              # knit all vignettes
#   Rscript vignettes/precompile.R sparse did   # knit selected only

library(knitr)
set.seed(54321)

VIGNETTES <- list(
  list(name = "ddml",
       src  = "vignettes/ddml.Rmd.txt",
       out  = "vignettes/ddml.Rmd",
       deps = "ddml",
       note = "takes too long for CRAN"),
  list(name = "example_BLP95",
       src  = "vignettes/articles/example_BLP95.Rmd.txt",
       out  = "vignettes/articles/example_BLP95.Rmd",
       deps = c("ddml", "AER", "hdm")),
  list(name = "example_401k",
       src  = "vignettes/articles/example_401k.Rmd.txt",
       out  = "vignettes/articles/example_401k.Rmd",
       deps = "ddml"),
  list(name = "sparse",
       src  = "vignettes/articles/sparse.Rmd.txt",
       out  = "vignettes/articles/sparse.Rmd",
       deps = c("ddml", "Matrix")),
  list(name = "stacking",
       src  = "vignettes/articles/stacking.Rmd.txt",
       out  = "vignettes/articles/stacking.Rmd",
       deps = "ddml"),
  list(name = "new_ml_wrapper",
       src  = "vignettes/articles/new_ml_wrapper.Rmd.txt",
       out  = "vignettes/articles/new_ml_wrapper.Rmd",
       deps = c("ddml", "gbm", "keras")),
  list(name = "did",
       src  = "vignettes/articles/did.Rmd.txt",
       out  = "vignettes/articles/did.Rmd",
       deps = c("ddml", "did"),
       note = "needs setwd for fig.path")
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  VIGNETTES <- Filter(function(v) v$name %in% args, VIGNETTES)
  if (length(VIGNETTES) == 0L) {
    known <- vapply(VIGNETTES, `[[`, character(1L), "name")
    stop("No matching vignettes. Known names: ",
         paste(known, collapse = ", "))
  }
}

message("=== Precompiling ", length(VIGNETTES), " vignette(s) ===\n")

for (v in VIGNETTES) {
  missing <- vapply(v$deps, function(pkg) {
    !requireNamespace(pkg, quietly = TRUE)
  }, logical(1L))
  if (any(missing)) {
    message("[SKIP] ", v$name, " -- missing: ",
            paste(v$deps[missing], collapse = ", "))
    next
  }

  message("[KNIT] ", v$name, " (", v$src, " -> ", v$out, ")")
  t0 <- proc.time()

  if (identical(v$name, "did")) {
    old_wd <- setwd("vignettes/articles")
    on.exit(setwd(old_wd), add = TRUE)
    knit("did.Rmd.txt", "did.Rmd")
    setwd(old_wd)
    on.exit(NULL)
  } else {
    knit(v$src, v$out)
  }

  elapsed <- (proc.time() - t0)[["elapsed"]]
  message("       done in ", round(elapsed, 1L), "s\n")
}

message("=== Session info ===")
print(utils::sessionInfo())
