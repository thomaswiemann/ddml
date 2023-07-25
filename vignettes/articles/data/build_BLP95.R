# Automobile market data from Berry, Levinsohn, Pakes (1995) as retrived from
# the reproduction excercise in Chernozhukov, Hansen, and Spindler (2015).
# [link](https://www.aeaweb.org/articles?id=10.1257/aer.p20151022)

library(R.matlab)
df <- readMat("vignettes/articles/data/BLP_data.mat")

# Create a single dataframe containing all variables
df$cdindex <- df$own.dummies <- NULL
BLP95 <- data.frame(matrix(unlist(df), nrow=2217, byrow=F))
colnames(BLP95) <- names(df)
BLP95$model.name <- readMat("vignettes/articles/data/BLP_data_str.mat")$model.name
BLP95 <- subset(BLP95, select = -c(const, model.id, product))

# Save
saveRDS(BLP95, file = "vignettes/articles/data/BLP95.rds")
