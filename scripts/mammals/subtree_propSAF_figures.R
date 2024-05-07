#### PACKAGES ####
library(ggplot2)

#### SET UP CLADES AND THEME ####
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#basic density theme
theme_density <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background= element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text.x = element_text(size = 7.5),
                       axis.text.y = element_text(size = 7.5),
                       axis.title = element_text(face = "bold",
                                                 size = 7.5),
                       legend.title = element_blank(),
                       legend.position = c(0.21,0.8),
                       legend.key = element_blank(),
                       legend.box = element_blank(),
                       legend.box.background = element_blank(),
                       legend.text = element_text(size=7.5),
                       legend.key.height= unit(0.1,"inch"),
                       legend.key.width = unit(0.1,"inch"))

for(i in 1:5){
  
  #### LOAD DATA ####
  hpd.intervals <- read.csv(paste0("../../outputs/mammals/subtrees/HPD_",clades[i],"_intervals.csv"))[,-1]
  raw.dat <- read.csv(paste0("../../outputs/mammals/subtrees/proportions_",clades[i],"_raw.csv"))[,-1]
  print(clades[i])
  print(paste0("Observed mean: ", mean(raw.dat$proportion[which(raw.dat$category == "Observed")])))
  print(paste0("Observed hpd: ", hpd.intervals$x[which(hpd.intervals$category == "Observed")]))
  print(paste0("Null mean: ",mean(raw.dat$proportion[which(raw.dat$category == "Null")])))
  print(paste0("Null hpd: ", hpd.intervals$x[which(hpd.intervals$category == "Null")]))
  
}
  #### PLOT ####
  SAF.overlap <- ggplot()+
    geom_density(aes(x=raw.dat$proportion,fill=raw.dat$category),
                 alpha=0.5)+
    geom_line(mapping=aes(x=hpd.intervals$x,y=hpd.intervals$y,
                          color=hpd.intervals$category),
              show.legend = F)+
    scale_y_continuous("Density")+
    scale_fill_viridis_d()+
    scale_color_viridis_d()+
    xlab("Proportion SAF")+
    labs(subtitle = paste0(clades[i]))+
    theme_density
  
  plot(SAF.overlap)
  
  #### SAVE PLOT ####
  ggsave(SAF.overlap,
         filename = paste0("../../figures/mammals/subtrees/observed_null_overlap_",clades[i],".pdf"),
         width = 3.5,
         height = 3.5,
         units = "in")
  
}



