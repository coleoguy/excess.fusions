dat <- read.csv("sizevsnumber.csv")
library(ggplot2)
taxa <- unique(dat$species)
rsq <- c()
for(i in 1:length(taxa)){
  curdat <- dat[dat$species==taxa[i],]
  foo <- summary(lm(curdat$Gene.number~curdat$bp.number))
  rsq[i] <- foo$adj.r.squared
}
names(rsq) <- taxa
ggplot(dat, aes(y=bp.number, x=Gene.number)) + 
  geom_point(aes(colour=species), position="jitter", alpha=0.5, size=3) + 
  coord_flip() + 
  theme_bw() + 
  scale_size(range=c(1, 3)) + xlab("Gene Number") + ylab("Chromosome Size") +
  annotate("text", y=590000000, x=3800, label= "R-Squared") +
  annotate("text", y=590000000, x=3500, label= "blue whale = 0.54") +
  annotate("text", y=590000000, x=3200, label= "camel = 0.49") +
  annotate("text", y=590000000, x=2900, label= "cow = 0.44") +
  annotate("text", y=590000000, x=2600, label= "domestic cat = 0.71") +
  annotate("text", y=590000000, x=2300, label= "horse = 0.80") +
  annotate("text", y=590000000, x=2000, label= "human = 0.68") +
  annotate("text", y=590000000, x=1700, label= "mouse = 0.51") +
  annotate("text", y=590000000, x=1400, label= "opossum = 0.97")
  

summary(lm(dat$Gene.number~dat$bp.number))




