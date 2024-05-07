#### LIBRARY ####
library(viridis)
library(ggplot2)

#### LOAD DATA ####
dat <- read.csv("Sex Chromosome Fusion Data.csv")
dat$Family[which(dat$Family == "Phy lostimidae")] <- "Phyllostimidae"
dat$fusind <- dat$fused.index/dat$total.autosomes

#### PLOT ####
fusedIndexHist <- ggplot(dat, aes(x=fusind, fill=Family)) + geom_histogram(bins=25) +
  theme_bw() + labs(x="Fusion Ratio") + 
  geom_vline(xintercept = .5) +
  xlim(0,1.05) + 
  theme(axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
        axis.text.y = element_text(size = 7.5),
        axis.title = element_text(face = "bold",
                                  size = 7.5),
        legend.text = element_text(size=7.5),
        legend.key.height= unit(0.075,"inch"),
        legend.key.width = unit(0.075,"inch"))

#hist(dat$fused.index/dat$total.autosomes)

#### SAVE ####
ggsave(fusedIndexHist,
       filename = paste0("./figS2.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")
ggsave(fusedIndexHist,
       filename = paste0("./figS2.jpg"),
       width = 3.5,
       height = 3.5,
       units = "in")

