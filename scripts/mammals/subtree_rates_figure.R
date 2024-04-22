#### LOAD PACKAGES ####
library(ggplot2)
library(phytools)
library(coda)

#### SET UP CLADES AND THEME ####
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

cladesForFigure <- c("Carnivora",
                     "Artiodactlya",
                     "Yangochiroptera",
                     "Rodentia",
                     "Primatomorpha")

theme_density <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background= element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text.x = element_text(size = 7.5),
                       axis.text.y = element_text(size = 7.5),
                       axis.title = element_text(face = "bold",
                                                 size = 7.5),
                       legend.title = element_blank(),
                       legend.position = c(0.81,0.9),
                       legend.key = element_blank(),
                       legend.box = element_blank(),
                       legend.box.background = element_blank(),
                       legend.text = element_text(size=7.5),
                       legend.key.height= unit(0.1,"inch"),
                       legend.key.width = unit(0.1,"inch"),
                       legend.spacing.y = unit(0.03,"inch"))

#### GET DEPTH OF EACH CLADE ####
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/tree.nex"),method="extend")

#Define edges of interest
edges <- c(181,448,796,1285,1718)

#Set up list
treeDepths <- list()

#Loop through edges of interest
for(i in 1:length(edges)){
  
  treeDepths[[i]] <- max(branching.times(extract.clade(tree,
                              node=tree$edge[edges[i],2]))) * 100
  
}

#### PREPARE DATA ####
#set up data frames
rawProps <- as.data.frame(matrix(NA,0,3))
hpdInts <- as.data.frame(matrix(NA,0,4))

#load in data for each clade
for( i in 1:5){
  load(paste0("../../outputs/mammals/mcmc/subtrees/hapauto/",clades[i],".RData"))
  
  #append to df
  rawProps <- rbind(rawProps,
                    cbind(model.mcmc.postburn$asc1/treeDepths[[i]],
                          model.mcmc.postburn$desc1/treeDepths[[i]],
                          cladesForFigure[i]))
  
  hpdInts <- rbind(hpdInts,
                   cbind(HPDinterval(as.mcmc(model.mcmc.postburn$asc1/treeDepths[[i]]))[1],
                         HPDinterval(as.mcmc(model.mcmc.postburn$desc1/treeDepths[[i]]))[1],
                         0-(0.75*i),
                         cladesForFigure[i]),
                   cbind(HPDinterval(as.mcmc(model.mcmc.postburn$asc1/treeDepths[[i]]))[2],
                         HPDinterval(as.mcmc(model.mcmc.postburn$desc1/treeDepths[[i]]))[2],
                         0-(0.75*i),
                         cladesForFigure[i]))

}

#name columns
colnames(rawProps) <- c("asc1","desc1","clade")
colnames(hpdInts) <- c("Xasc1","Xdesc1","Y","clade")
rawProps$clade <- factor(rawProps$clade,
                                levels=unique(rawProps$clade))
hpdInts$clade <- factor(hpdInts$clade,
                            levels=unique(hpdInts$clade))

#### FIGURE ####
fissionOverlap <- ggplot()+
  geom_density(aes(x=as.numeric(rawProps$asc1),fill=rawProps$clade),
               alpha=0.5,adjust=1.5)+
  geom_line(mapping=aes(x=as.numeric(hpdInts$Xasc1),y=as.numeric(hpdInts$Y),
                        color=hpdInts$clade),
            show.legend = F,
            size=1)+
  #xlim(0.01,0.11)+
  scale_y_continuous("Density")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  xlab("Fission rate (per MY)")+
  guides(fill = guide_legend(byrow = TRUE))+
  theme_density

fusionOverlap <- ggplot()+
  geom_density(aes(x=as.numeric(rawProps$desc1),fill=rawProps$clade),
               alpha=0.5,adjust=1.5)+
  geom_line(mapping=aes(x=as.numeric(hpdInts$Xdesc1),y=as.numeric(hpdInts$Y),
                        color=hpdInts$clade),
            show.legend = F,
            size=1)+
  #xlim(0.01,0.11)+
  scale_y_continuous("Density")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  xlab("Fusion rate (per MY)")+
  guides(fill = guide_legend(byrow = TRUE))+
  theme_density

#### SAVE FIGURES ####
ggsave(fissionOverlap,
       filename = paste0("../../figures/mammals/subtrees/subtree_fission_overlap.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")
ggsave(fusionOverlap,
       filename = paste0("../../figures/mammals/subtrees/subtree_fusion_overlap.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")