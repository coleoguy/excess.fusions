#### PACKAGES ####
library(ggplot2)

#### SET UP CLADES AND THEME ####
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primates")

theme_violin_bar <- theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background= element_blank(),
                          axis.line = element_line(colour = "black"),
                          axis.text.x = element_text(size = 7,angle=45,vjust=1,hjust=1),
                          axis.text.y = element_text(size = 7.5),
                          axis.title = element_text(face = "bold",
                                                    size = 7.5),
                          axis.title.x = element_blank(),
                          legend.title = element_blank(),
                          legend.position = c(0.25,0.85),
                          legend.key = element_blank(),
                          legend.box = element_blank(),
                          legend.box.background = element_blank(),
                          legend.text = element_text(size=7.5),
                          legend.key.height= unit(0.15,"inch"),
                          legend.key.width = unit(0.15,"inch"))

#### SET UP DATA ####
clades_observed <- data.frame()
clades_null <- as.data.frame(matrix(ncol=2,nrow=5))

#loop through csv files for clades
for(i in 1:5){
  #get raw data
  raw.dat <- read.csv(paste0("../../outputs/SAF_proportions/subtrees/proportions_",clades[i],"_raw.csv"))[,-1]
  
  clades_observed <- rbind(clades_observed,
                           cbind(raw.dat[1:100,1],
                                 clades[i]))
  
  clades_null[i,] <- c(mean(raw.dat[101:200,1]),clades[i])
}

colnames(clades_null) <- c("mean","clade")
colnames(clades_observed) <- c("prop","clade")
clades_observed$prop <- as.numeric(clades_observed$prop)
clades_null$mean <- as.numeric(clades_null$mean)

clades_observed$clade <- factor(clades_observed$clade,
                                   levels=unique(clades_observed$clade))
clades_null$clade <- factor(clades_null$clade,
                               levels=unique(clades_null$clade))
levels(clades_observed$clade) <- levels(clades_null$clade) <- c("Carnivora",
                                                                "Artiodactlya",
                                                                "Yangochiroptera",
                                                                "Rodentia",
                                                                "Primates")

#### PLOT ####
subtree_scatter <- ggplot()+
  geom_jitter(data=clades_observed,
              mapping = aes(x=clades_observed$clade,
                  y=clades_observed$prop,
                  col=clades_observed$clade),
              alpha=0.5,size=0.75,position=position_jitterdodge())+
  geom_point(mapping=aes(x=factor(clades_null$clade),y=clades_null$mean,
                        color=clades_null$clade),
            show.legend = F,shape=95,size=7.5)+
  scale_y_continuous(name="Proportion SAF")+
  xlab("Clade")+
  scale_color_viridis_d()+
  theme_violin_bar

ggsave(subtree_scatter,
       filename = paste0("./fig5.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")
ggsave(subtree_scatter,
       filename = paste0("./fig5.jpg"),
       width = 3.5,
       height = 3.5,
       units = "in")














