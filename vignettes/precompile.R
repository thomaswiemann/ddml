# Articles that depend on other packages are precompiled
library(knitr)

# depends on keras
knit("vignettes/articles/new_ml_wrapper.Rmd.txt",
     "vignettes/articles/new_ml_wrapper.Rmd")

# depends on AER, hdm
knit("vignettes/articles/example_BLP95.Rmd.txt",
     "vignettes/articles/example_BLP95.Rmd")

# depends on Matrix
knit("vignettes/articles/example_401k.Rmd.txt",
     "vignettes/articles/example_401k.Rmd")

# depends on Matrix
knit("vignettes/articles/sparse.Rmd.txt",
     "vignettes/articles/sparse.Rmd")
