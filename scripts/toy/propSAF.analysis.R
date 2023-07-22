# This script calculates the observed and null proportions of SA-fusions
# A plot is produced which compares the two distributions

#### PACKAGES ####

library(phytools)
library(evobiR)
library(coda)
library(ggplot2)
source("../functions.R")

#### LOAD IN DATA ####
#replace toy with your clade directory 
dat <- read.csv("../../data/toy/dat.csv",
                as.is=T)
load(file="../../outputs/toy/hists.summarized.RData") 
trans <- hists.summarized$count

#### EXTRACT COLUMNS ASSOCIATED WITH TRANSITIONS OF INTEREST ####

#Get column names
names <- colnames(trans[,2:ncol(trans)])

#Split by commas
split.names <- strsplit(names,",")

#Create vectors of rows
SAF.cols <- c()
AAF.cols <- c()

for(i in 1:length(split.names)){
  if(as.numeric(split.names[[i]][1]) == as.numeric(split.names[[i]][2])-(max(dat$hapauto) - min(dat$hapauto)) &&
     as.numeric(split.names[[i]][2]) >= (max(dat$hapauto) - min(dat$hapauto) + 2) &&
     as.numeric(split.names[[i]][1]) <= (max(dat$hapauto) - min(dat$hapauto) + 1)){
    SAF.cols <- c(SAF.cols,i+1)
  } else if(as.numeric(split.names[[i]][1]) == as.numeric(split.names[[i]][2])+1 &&
            as.numeric(split.names[[i]][1]) != (max(dat$hapauto) - min(dat$hapauto) + 2)){
    AAF.cols <- c(AAF.cols,i+1)
  }
}

#### SUM SAF AND AA COUNTS ####
#Get rowSums
SAF.counts <- rowSums(trans[,SAF.cols])
AAF.counts <- rowSums(trans[,AAF.cols])

#Divide
obspropSAF <- SAF.counts/(AAF.counts + SAF.counts)

#### NULL SAF ####

#null proportions vector
expSA <- c()

for(i in 1:100){
  
  #print iterations
  print(paste0("Tree ",i))
  
  #get current times
  times <- hists.summarized$times[i,1:(ncol(hists.summarized$times) - 1)]/
    hists.summarized$times[i,ncol(hists.summarized$times)]
  
  #vector for current map proportions
  expSA.cur <- c()
  
  #loop through states
  for(j in 1:length(times)){
    if(j > length(times)/2){
      #convert state to diploid autosome
      Da <- (j - max(dat$hapauto) + min(dat$hapauto) - 1) * 2
      
      #Calculate proportion SAF
      expSA.cur[j] <- Pfsa2(Da = Da,
                            scs = "XY", # change "XY" to your fused scs
                            mud = 0.5) * times[j]
    } else {
      #convert state to diploid autosome
      Da <- (j + min(dat$hapauto) - 1) * 2
      
      #Calculate proportion SAF
      expSA.cur[j] <- Pfsa2(Da = Da,
                            scs = "XO", # change "XO" to your unfused scs
                            mud = 0.5) * times[j]
      
    }
    names(expSA.cur)[j] <- Da
  }
  
  expSA[i] <- sum(expSA.cur)
}

#### SUMMARIZE AND PLOT ####
#Dataframe for raw proportions
raw.props <- as.data.frame(matrix(NA,
                                  nrow=200,
                                  ncol=2,
                                  byrow = F))

colnames(raw.props) <- c("proportion","category")

raw.props[1:100,] <- cbind(obspropSAF,
                           rep("Observed",100))

raw.props[101:200,] <- cbind(expSA,
                             rep("Null",100))

raw.props$proportion <- as.numeric(raw.props$proportion)

#Dataframe for results
hpd.intervals <- as.data.frame(matrix(NA,
                                      nrow=4,
                                      ncol=3,
                                      byrow = F))

colnames(hpd.intervals) <- c("x",
                             "y",
                             "category")
hpd.intervals$category <- c(rep("Observed",2),
                            rep("Null",2))
#adjust these y values as needed
hpd.intervals$y <- c(rep(-0.5,2),
                     rep(-0.75,2))

hpd.intervals$x <- c(HPDinterval(as.mcmc(obspropSAF)),
                     HPDinterval(as.mcmc(expSA)))

# plot
theme_density <- theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.background= element_blank(),
                       axis.line = element_line(colour = "black"),
                       axis.text.x = element_text(size = 15),
                       axis.text.y = element_text(size = 15),
                       axis.title = element_text(face = "bold",
                                                 size = 15),
                       legend.title = element_blank(),
                       plot.title = element_text(face = "bold",
                                                 size = 17,
                                                 hjust=0.5))

SAF.overlap <- ggplot()+
  geom_density(aes(x=raw.props$proportion,fill=raw.props$category),
               alpha=0.5)+
  geom_line(mapping=aes(x=hpd.intervals$x,y=hpd.intervals$y,
                        color=hpd.intervals$category),
            show.legend = F)+
  scale_y_continuous("Density")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  xlab("Proportion SAF")+
  theme_density

plot(SAF.overlap)

#### SAVE ####
write.csv(hpd.intervals,
          paste0("../../outputs/toy/HPD_intervals.csv"),
          quote=F,
          row.names=T)
write.csv(raw.props,
          paste0("../../outputs/toy/proportions_raw.csv"),
          quote=F,
          row.names=T)
ggsave(SAF.overlap,
       filename = paste0("../../figures/toy/observed_null_overlap.pdf"),
       width = 7,
       height = 7,
       units = "in")










