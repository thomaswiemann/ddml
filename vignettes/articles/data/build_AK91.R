# Data on the Returns to Schooling from Angrist & Krueger (1991)

library(readstata13)

# Import data and rename variables
dat <- read.dta13("vignettes/articles/data/ak91.dta")
colnames(dat)[c(1,2,4,5,6,9,10,11,12,13,16,17,18,19,20,21,24,25,27)] <-
  c('AGE','AGEQ','EDUC','ENOCENT','ESOCENT','LWKLYWGE','MARRIED','MIDATL',
    'MT','NEWENG','CENSUS','POB','QOB','RACE','SMSA','SOATL','WNOCENT',
    'WSOCENT','YOB')

# Only keep male 1930-1939 sample
men_1930_1939 <- dat[,"YOB"] >= 30 & dat[,"YOB"] <= 39
dat <- dat[men_1930_1939,]

# Only keep log-weekly wage, yrs of education, and QOB, YOB, and POB
AK91 <- data.frame(sapply(dat[, c("QOB","YOB","POB")], as.factor))
AK91[c("LWKLYWGE","EDUC")] <- dat[, c("LWKLYWGE","EDUC")]

# Save
saveRDS(AK91, file = "vignettes/articles/data/AK91.rds")
