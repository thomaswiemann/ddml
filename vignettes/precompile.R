# Articles that depend on other packages are precompiled
library(knitr)

# depends on keras
knit("vignettes/articles/new_ml_wrapper.Rmd.lol",
     "vignettes/articles/new_ml_wrapper.Rmd")

# depends on AER, hdm
knit("vignettes/articles/example_BLP95.Rmd.lol",
     "vignettes/articles/example_BLP95.Rmd")
