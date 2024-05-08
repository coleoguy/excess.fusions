library(viridis)
library(ggplot2)

# read in data
dat <- read.csv("Sex Chromosome Fusion Data.csv")[,1:5]

# calculate fusion ratios
dat$fusind <- dat$fused.index/dat$total.autosomes

# generate a null expectation. each species has a limited range 
# that can be as low as 1/n to 1.
null <- c()
for(i in 1:1000){
  curvals <- c()
  for(j in 1:nrow(dat)){
    curvals[j] <- sample(x=1:(dat$total.autosomes[j]+1), 1)
  }
  null[i] <- mean(curvals/(dat$total.autosomes+1))
}

ggplot(dat, aes(x=fusind, fill=Family)) + geom_histogram(bins=25) +
  theme_bw() + labs(x="Fusion Ratio") + 
  geom_vline(xintercept = mean(null)) +
  xlim(0,1.05)

small <- sum(dat$fusind > mean(null))
binom.test(39, 51)
