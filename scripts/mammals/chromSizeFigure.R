#### PACKAGES ####
#### LOAD DATA ####
dat <- read.csv("../../data/mammals/chromSizeSAF.csv")[,4:6]
colnames(dat) <- c("fused","total","ratio")

#### CALCULATE RATIO ####
dat$ratio <- dat$fused/dat$total

hist(dat$ratio)
