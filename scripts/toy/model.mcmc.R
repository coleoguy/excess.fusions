#### PACKAGES ####
library(phytools)
library(ape)
library(chromePlus)
library(diversitree)

#### BUILD TEST TREE AND TEST DATASET ####
# tree <- rcoal(20)
# dat <- as.data.frame(matrix(data=c(tree$tip.label,
#                                    sample(x=1:10,size=20,replace=T),
#                                    sample(c("fused","unfused"),size=20,replace=T)),
#                             ncol=3,
#                             byrow = F)) 
# colnames(dat) <- c("tree.name","hapauto","scs")
# write.nexus(tree,file="../../data/toy/tree.nex")
# write.csv(dat,file="../../data/toy/dat.csv",row.names = F)

########## START HERE #############
dat <- read.csv("../../data/toy/dat.csv",
                as.is=T)
tree <- force.ultrametric(read.nexus("../../data/toy/tree.nex"))

#make sure your tip states and tree are pruned to one another here, and that your
#tip states are in the same order that your tree tip labels are in

#### MODEL ####

#Convert codedSCS column to right format
for(i in 1:nrow(dat)){
  if(dat$scs[i] ==  "fused"){
    dat$scs[i] <- 0
  } else {
    dat$scs[i] <- 1
  }
}

dat$scs <- as.numeric(dat$scs)

#Convert to data matrix
data.matrix <- datatoMatrix(dat,
                            c(min(dat$hapauto),max(dat$hapauto)),
                            hyper = T)

#Make mkn model
model <- make.mkn(tree,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#Constrain model
model.con <- constrainMkn(data.matrix,
                          model,
                          hyper = T,
                          polyploidy = F,
                          verbose = T,
                          constrain = list(saf.model=T,
                                           sym.hyperstates=T))

parMat <- maat <- model.con$`parameter matrix`



