# Articles that depend on other packages are precompiled
library(knitr)

# depends on keras
knit("vignettes/articles/new_ml_wrapper.Rmd.txt",
     "vignettes/articles/new_ml_wrapper.Rmd")

# depends on AER, hdm
knit("vignettes/articles/example_BLP95.Rmd.txt",
     "vignettes/articles/example_BLP95.Rmd")

# depends on SIPP data
knit("vignettes/articles/example_401k.Rmd.txt",
     "vignettes/articles/example_401k.Rmd")

# depends on AK91 data
knit("vignettes/articles/sparse.Rmd.txt",
     "vignettes/articles/sparse.Rmd")

#
knit("vignettes/articles/stacking.Rmd.txt",
     "vignettes/articles/stacking.Rmd")

# takes too long to run on cran...
knit("vignettes/ddml.Rmd.txt",
     "vignettes/ddml.Rmd")

# depends on did, change wd for figure references...
setwd("vignettes/articles")
knit("did.Rmd.txt",
     "did.Rmd")
