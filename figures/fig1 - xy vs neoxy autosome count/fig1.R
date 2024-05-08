# This script compares haploid autosome number of different mammlian families
# A plot is generated which visualizes these comparisons
# This corresponds with figure 1 in the manuscript *INSERT MANUSCRIPT TITLE *

#### PACKAGES ####
library(ggplot2)

#### LOAD DATA ####
dat <- read.csv("../../data/chromes/mammal_chroms_working.csv",
                as.is=T)[,c(2,3,4,5,12,13)]
dat <- dat[-which(dat$codedSCS == "unclear"),]

#### GATHER CHROMOSOMES ####
#factor columns
dat$Family <- factor(dat$Family)
dat$codedSCS <- as.factor(dat$codedSCS)

#determine subset of families to check
families <- names(which(summary(dat$Family) >= 10))
families <- families[-length(families)]

#create df
meanHapauto <- as.data.frame(matrix(NA,nrow = 2*length(families),ncol=3))
index <- 1

#loop through families
for(i in 1:length(families)){
  #subset dataset
  subDat <- subset(dat,Family==families[i])
  
  #skip if no fused scs
  if(length(which(subDat$codedSCS == "fused")) == 0){
    next()
  }
  
  #fill dataset
  meanHapauto[c(index,index+1),1] <- families[i]
  meanHapauto[c(index,index+1),2] <- c("XY","NeoXY")
  meanHapauto[c(index,index+1),3] <- c(mean(subDat$hapauto[which(subDat$codedSCS == "unfused")],
                                            na.rm=T),
                                       mean(subDat$hapauto[which(subDat$codedSCS == "fused")],
                                            na.rm=T))
  #increment index
  index <- index+2
}

#remove NA rows
meanHapauto <- meanHapauto[-which(is.na(meanHapauto[,1])),]

#name/factor columns
colnames(meanHapauto) <- c("family","scs","hapauto")
meanHapauto$family <- factor(meanHapauto$family)
meanHapauto$scs <- factor(meanHapauto$scs,levels=c("XY","NeoXY"))

#### DETERMINE LINETYPE ####

#determine if fused less than or greater than unfused
meanHapauto$line <- NA

for(i in levels(meanHapauto$family)){
  if(meanHapauto$hapauto[which(meanHapauto$family==i & meanHapauto$scs=="NeoXY")] <=
     meanHapauto$hapauto[which(meanHapauto$family==i & meanHapauto$scs=="XY")]){
    meanHapauto$line[which(meanHapauto$family==i)] <- "filled"
  } else {
    meanHapauto$line[which(meanHapauto$family==i)] <- "dashed"
  }
}

meanHapauto$line <- factor(meanHapauto$line,levels=c("filled","dashed"))

#### PLOT ####

meanHapPlot <- ggplot(data=meanHapauto,aes(x=scs,y=hapauto)) +
  geom_point(aes(group=family,color=family)) +
  geom_line(aes(group=family,color=family,linetype=line),show.legend = F) + 
  scale_y_continuous("Mean haploid autosome number")+
  scale_x_discrete("Sex chromosome system") +
  scale_color_viridis_d()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 7.5),
        axis.text.y = element_text(size = 7.5),
        axis.title = element_text(face = "bold",
                                  size = 7.5),
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size=7.5),
        legend.key.height= unit(0.15,"inch"),
        legend.key.width = unit(0.03,"inch"))

#### SAVE PLOT ####
ggsave(meanHapPlot,
       filename = paste0("./fig1.pdf"),
       width = 3.5,
       height = 3.5,
       units = "in")
ggsave(meanHapPlot,
       filename = paste0("./fig1.jpg"),
       width = 3.5,
       height = 3.5,
       units = "in")

         
         
         
         
         