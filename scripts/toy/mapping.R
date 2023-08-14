# This script performs stochastic mapping using the karyotypic tip states,
# phylogeny, and Qmatrix

#### PACKAGES ####
library(phytools)
source("../functions.R")

#### LOAD DATA ####
#replace "toy" here with your clade directory
dat <- read.csv("../../data/toy/dat.csv",
                as.is=T)
tree <- read.tree("../../data/toy/tree.nex")
mat <- as.matrix(read.csv("../../data/toy/transition_matrix.csv",
                          as.is=T,header = T))
Qmat <- as.matrix(read.csv("../../data/toy//Q_matrix.csv",
                           as.is=T,header=T))

#### CONVERT HAPAUTO + SCS INTO SIMULAION STATE ####

#make sure your tip states and tree are pruned to one another here, and that your
#tip states are in the same order that your tree tip labels are in
hapauto.min <- min(dat$hapauto)
hapauto.max <- max(dat$hapauto)

#transform hapauto and scs
# replace "fused" with whatever sex chromosome system you think is the fused state
# e.g, XY for XO/XY and NeoXY for XY/NeoXY
for(i in 1:nrow(dat)){
  if(dat[i,3] == "fused"){
    dat$sim.state[i] <- dat$hapauto[i] - (hapauto.min -1) + hapauto.max
  } else {
    dat$sim.state[i] <- dat$hapauto[i] - (hapauto.min -1)
  }
}

#### BUILD DATA MATRIX ####
#data.matrix
data.matrix <- matrix(0,nrow=nrow(dat),
                      ncol=ncol(mat))

colnames(data.matrix) <- 1:ncol(mat)
rownames(data.matrix) <- dat$tree.name

#fill data matrix
for(i in 1:nrow(data.matrix)){
  data.matrix[i,dat$sim.state[i]] <- 1
}

#column and rownames for mat/Qmat
rownames(mat) <-colnames(mat) <- rownames(Qmat) <- colnames(Qmat) <- 1:ncol(Qmat)

#### STOCHASTIC MAPPING ####
hists <- make.simmap2(tree = tree,
                      x = data.matrix,
                      model = mat,
                      nsim = 100,
                      Q = Qmat,
                      rejmax = 10000000,
                      rejint = 1000000,
                      pi="fitzjohn",
                      monitor=T)

#### FIX SIMMAP ####
dat.for.fixing <- dat[,-c(2,3)]
hists <- fix.simmap(hists,dat.for.fixing,mat)

#### SUMARIZE AND SAVE OUPUTS ####
hists.summarized <- describe.simmap2(hists)
save(hists, file="../../outputs/toy/hists.RData")
save(hists.summarized, file = "../../outputs/toy/hists.summarized.RData")




