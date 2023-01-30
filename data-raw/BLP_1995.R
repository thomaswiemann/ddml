#' Automobile market data from Berry, Levinsohn, Pakes (1995)
#'
#' Automobile market data from Berry, Levinsohn, Pakes (1995) as retrived from
#' the reproduction excercise in Chernozhukov, Hansen, and Spindler (2015).

# Load the data from CHS2015 matlab files
library(R.matlab)
dat <- df <- readMat("data-raw/BLP_data.mat")

# Create a single dataframe containing all variables
df$cdindex <- df$own.dummies <- NULL
BLP_1995 <- data.frame(matrix(unlist(df), nrow=2217, byrow=F))
colnames(BLP_1995) <- names(df)
BLP_1995$model.name <- readMat("data-raw/BLP_data_str.mat")$model.name
BLP_1995 <- subset(BLP_1995, select = -c(const, model.id, product))

# Use data in package
usethis::use_data(BLP_1995, overwrite = TRUE)
