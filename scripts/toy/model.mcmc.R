# This script builds our model for karyotypic evolution and estimates
# rates for the parameters using an MCMC

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
# tree$edge.length <- tree$edge.length/max(branching.times(tree))
# write.nexus(tree,file="../../data/toy/tree.nex")
# write.csv(dat,file="../../data/toy/dat.csv",row.names = F)

########## START HERE #############
#change to your tree and csv files, change toy directory to your clade directory
dat <- read.csv("../../data/toy/dat.csv",
                as.is=T)
#make sure to check your tree is utlrametric
tree <- force.ultrametric(read.nexus("../../data/toy/tree.nex"))

#the follwing lines convert the edge lengths of the tree to unit length and 
#save the unit length tree
tree$edge.length <- tree$edge.length/max(branching.times(tree))
write.nexus(tree,file="../../data/toy/tree.nex") #change the directory accordingly

#make sure your tip states and tree are pruned to one another here, and that your
#tip states are in the same order that your tree tip labels are in

#### BUILD MODEL ####

#Convert codedSCS column to right format
# replace "fused" with whatever sex chromosome system you think is the fused state
# e.g, XY for XO/XY and NeoXY for XY/NeoXY
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

#extract parameter matrix
parMat <- mat <- model.con$`parameter matrix`

#### ESTIMATE RATES ####
#test run MCMC
temp.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                               x.init=runif(4,0,1),
                               prior=make.prior.exponential(r=0.5),
                               #upper=c(100,100,100,100),
                               nsteps = 100,
                               w=1)

#extract tuning params
temp.mcmc <- temp.mcmc[-c(1:50), ]
w <- diff(sapply(temp.mcmc[2:5],
                 quantile, c(.05, .95)))

#run MCMC
model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                x.init=runif(4,0,1),
                                prior=make.prior.exponential(r=0.5),
                                #upper=c(100,100,100,100),
                                nsteps = 500,
                                w=w)

#### BUILD QMATRIX ####
model.mcmc.postburn <- model.mcmc[450:500,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1),
            mean(model.mcmc.postburn$tranSAF),
            mean(model.mcmc.postburn$tranRo))

names(params) <- colnames(model.mcmc[,2:5])

#substitute rates into matrix
parMat[parMat == "asc1"] <- params[1]
parMat[parMat == "desc1"] <- params[2]
parMat[parMat == "tranSAF"] <- params[3]
parMat[parMat == "tranRo"] <- params[4]
mat[mat == "asc1"] <- 1
mat[mat == "desc1"] <- 2
mat[mat == "tranSAF"] <- 3
mat[mat == "tranRo"] <- 4

#convert to numeric
parMat <- as.data.frame(parMat)
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)
mat <- as.data.frame(mat)
mat <- sapply(mat[,1:ncol(mat)],as.numeric)

#set diagonals
diag(parMat) <- diag(mat) <- 0
diag(parMat) <- -rowSums(parMat)

#### SAVE MATRICS ####
#change toy to your clade
write.csv(parMat,
          paste0("../../data/toy/Q_matrix.csv"),
          row.names=F,quote=F)
write.csv(mat,
          paste0("../../data/toy/transition_matrix.csv"),
          row.names=F,quote=F)

