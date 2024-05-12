dat <- read.csv("sizevsnumber.csv")
library(ggplot2)
ggplot(dat, aes(y=bp.number, x=Gene.number)) + 
  geom_point(aes(colour=species), position="jitter", alpha=0.5, size=3) + 
  coord_flip() + 
  theme_bw() + 
  scale_size(range=c(1, 3)) + xlab("Gene Number") + ylab("Chromosome Size")
summary(lm(dat$Gene.number~dat$bp.number))
