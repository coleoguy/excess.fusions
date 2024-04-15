#### LOAD PACKAGES ####
library(ape)
library(phytools)
library(viridis) 
source("../functions.R")

#### LOAD IN DATA ####
dat <- read.csv("../../data/mammals/chromes/dat.csv",
                as.is=T)[,c(1,3)]
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/tree.nex"),method = "extend")
mat <- as.matrix(read.csv("../../data/mammals/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header = T))
# For analysis of adjusted prior rates
# Qmat <- as.matrix(read.csv("../../data/mammals/transition_matrix/Q_matrix_hapauto_adjusted_prior.csv",
#                            as.is=T,header=T)) 

# For analysis of final rates
# Qmat <- as.matrix(read.csv("../../data/mammals/transition_matrix/Q_matrix_hapauto_final.csv",
#                            as.is=T,header=T)) 

#scale tree to unit rate
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#### BUILD DATA MATRIX ####

#data.matrix
data.matrix <- matrix(0,nrow=nrow(dat),
                      ncol=max(dat$hapauto) - 1)

colnames(data.matrix) <- 1:49
rownames(data.matrix) <- dat$tree.name

#fill data matrix
for(i in 1:nrow(data.matrix)){
  data.matrix[i,as.numeric(dat$hapauto[i]) - 1] <- 1
}

#column and rownames for mat/Qmat
rownames(mat) <- 1:49
colnames(mat) <- 1:49
rownames(Qmat) <- 1:49
colnames(Qmat) <- 1:49

####PERFORM TEST STOCHASTIC MAP ####
testHist <- make.simmap2(tree = tree,
                      x = data.matrix,
                      model = mat,
                      nsim = 1,
                      Q = Qmat,
                      rejmax = 1000000,
                      rejint = 100000,
                      pi="fitzjohn",
                      monitor=T)

#fix stochastic maps
dat$sim.state <- dat$hapauto - 1
datForFixing <- dat[,-2]
testHistFixed <- fix.simmap(testHist,datForFixing,mat)


#### SUBSET TREE TO JUST ARTIODACTYLS ####

#subset map
hist.subset <- extract.clade.simmap(tree = testHistFixed,
                                    node=testHistFixed$edge[796,2])

dat <- dat[which(dat$tree.name %in% hist.subset$tip.label),]

#plot
cols <- c(viridis(49))
names(cols) <- c(1:49)
plotSimmap(hist.subset,col=cols)

#### GET NODE STATES ####

#vector to store
node.states <- c()

for(i in 117:231){
  
  #get edge which descends from node
  desc.edge <- which(hist.subset$edge[,1] == i)[1]
  
  #get initial state of descendent edge
  node.state <- names(hist.subset$maps[[desc.edge]][1])
  
  #check that it matches with last state of leading edge
  if(i != 117){
    
    lead.edge <- which(hist.subset$edge[,2] == i)
    
    lead.state <- names(hist.subset$maps[[lead.edge]][length(hist.subset$maps[[lead.edge]])])
    
    if(lead.state != node.state){
      print("warning: node state does not match last state of leading edge")
    }
  }
  node.states[i-116] <- node.state
}

#### PLOT ####
plotSimmap(hist.subset,col=cols,fsize=0.1,lwd=0.1,pts=T)
plotSimmap(hist.subset,col=cols,fsize=0.0000001,lwd=0.1)
nodelabels(node.states,frame="none",cex=0.6)
tiplabels(dat$hapauto,frame="none",cex=0.6)
save(testHist, file = "../../outputs/mammals/testHist.RData")
save(testHistFixed, file = "../../outputs/mammals/testHistFixed.RData")




