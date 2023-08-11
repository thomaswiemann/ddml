#' 401k data from the Survey of Income and Program Participation (SIPP) 1991.
#'     Data taken from thereproduction exercise in Chernozhukov et al. (2018).
#'     [link](https://academic.oup.com/ectj/article/21/1/C1/5056401s%7D)

library(readstata13)
SIPP91 <- read.dta13("vignettes/articles/data/sipp1991.dta")

# Save
saveRDS(SIPP91, file = "vignettes/articles/data/SIPP91.rds")
