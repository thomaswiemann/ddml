#' 401k data from
#'
#' 401k data

library(readstata13)
SIPP_1991 = read.dta13("data-raw/sipp1991.dta")

# Use data in package
usethis::use_data(SIPP_1991, overwrite = TRUE)
