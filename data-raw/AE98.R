#' Random subsample from the data of Angrist & Evans (1991).
#'
#' Random subsample of 10,000 observations from the Angrist & Evans (1991) data.

# Set seed for subsample selection
set.seed(9265001)

# Load full data
library(readstata13)
AE98 <- read.dta13("data-raw/ae98_fulldata.dta")
AE98 <- AE98[AE98[, "msample"]==1,]

# Select subsample and variables
subsample <- sample(1:nrow(AE98), 5000)
AE98 <- AE98[subsample, ]
AE98 <- AE98[, c("workedm","weeksm1","hourswm", "morekids", "samesex",
                          "agem","agefstm","blackm",
                          "hispm","othracem","educm",
                          "boy1st", "boy2nd")]
colnames(AE98) <- c("worked","weeksw","hoursw", "morekids", "samesex",
                    "age","agefst","black",
                    "hisp","othrace","educ",
                    "boy1st", "boy2nd")
AE98 <- sapply(AE98, as.numeric)

# Use data in package
usethis::use_data(AE98, overwrite = TRUE)
