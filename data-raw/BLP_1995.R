#' Data from Berry, Levinsohn, Pakes (1995)
#'
#' Data from Berry, Levinsohn, Pakes (1995) as retrived from the reproduction
#' excercise in Chernozhukov, Hansen, and Spindler (2015)
#'
#'
#' @format A data frame with 2217 rows and 17 variables
#' \describe{
#'   \item{model.id}{ }
#'   \item{id}{ }
#'   ...
#' }
#' @source \url{/https://www.aeaweb.org/articles?id=10.1257/aer.p20151022/}

# Load the data from CHS2015 matlab files
library(R.matlab)
dat <- df <- readMat("data-raw/BLP_data.mat")

# Create a single dataframe containing all variables
df$cdindex <- df$own.dummies <- NULL
BLP_1995 <- data.frame(matrix(unlist(df), nrow=2217, byrow=F))
colnames(BLP_1995) <- names(df)
BLP_1995$model.name <- readMat("data-raw/BLP_data_str.mat")$model.name
BLP_1995$own.dummies <- dat$own.dummies

# Use data in package
usethis::use_data(BLP_1995, overwrite = TRUE)
