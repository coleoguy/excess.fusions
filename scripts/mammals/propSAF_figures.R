#### PACKAGES ####
library(ggplot2)

#### LOAD DATA ####
hpd.intervals <- read.csv("../../outputs/mammals/HPD_intervals.csv")[,-1]
raw.dat <- read.csv("../../outputs/mammals/proportions_raw.csv")[,-1]

#### OVERLAP PLOTS ####
theme_density <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background= element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text.x = element_text(size = 7.5),
                       axis.text.y = element_text(size = 7.5),
                       axis.title = element_text(face = "bold",
                                                 size = 7.5),
                       legend.title = element_blank(),
                       legend.position = c(0.21,0.9),
                       legend.key = element_blank(),
                       legend.box = element_blank(),
                       legend.box.background = element_blank(),
                       legend.text = element_text(size=7.5),
                       legend.key.height= unit(0.1,"inch"),
                       legend.key.width = unit(0.1,"inch"),
                       legend.spacing.y = unit(0.03,"inch"))

SAF.overlap <- ggplot()+
  geom_density(aes(x=raw.dat$proportion,fill=raw.dat$category),
               alpha=0.5,adjust=1.5)+
  geom_line(mapping=aes(x=hpd.intervals$x,y=hpd.intervals$y,
                        color=hpd.intervals$category),
            show.legend = F,
            size=1)+
  xlim(0.01,0.11)+
  scale_y_continuous("Density")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  xlab("Proportion SAF")+
  guides(fill = guide_legend(byrow = TRUE))+
  theme_density

plot(SAF.overlap)

#### SAVE PLOT ####
ggsave(SAF.overlap,
       filename = paste0("../../figures/mammals/observed_null_overlap.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")
ggsave(SAF.overlap,
       filename = paste0("../../figures/mammals/observed_null_overlap.jpg"),
       width = 3.5,
       height = 3.5,
       units = "in")

