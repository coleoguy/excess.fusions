#### PACKAGES ####
library(phytools)
library(doSNOW)
library(viridis)
library(doSNOW)
source("functions.R")

#### LOAD DATA ####
dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,3)]
tree <- force.ultrametric(read.nexus("../data/trees/tree.nex"),method="extend")
mat <- as.matrix(read.csv("../data/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header = T))
Qmat <- as.matrix(read.csv("../data/transition_matrix/Q_matrix_hapauto_final.csv",
                           as.is=T,header=T))

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

#### SET UP PARALLELIZATION ####

# Define number of clusters
nClust <- 100

# Set up clusters, print intermediates to 
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

#### STOCHASTIC MAPPING ####
hists <- foreach(i=1:100,
                .verbose = T,
                .packages = c("phytools","maps","ape")) %dopar% {
                  make.simmap2(tree = tree,
                     x = data.matrix,
                     model = mat,
                     nsim = 1,
                     Q = Qmat,
                     rejmax = 1000000,
                     rejint = 100000,
                     pi="fitzjohn",
                     monitor=T,
                     parallel=c(i,100))
                }

#Close cluster connection
stopCluster(cl)

#Set class of sim.out
class(hists) <- c("multiSimmap","multiPhylo")

#### FIX STOCHASTIC MAPS ####
dat$sim.state <- dat$hapauto - 1
dat.for.fixing <- dat[,-2]
hists.fixed <- fix.simmap(hists,dat.for.fixing,mat)

#### SUMMARIZE STOCHASTIC MAPS ####
cols <- c(viridis(49))
names(cols) <- c(1:49)
plotSimmap(hists.fixed[[1]],col=cols,fsize = 0.05,lwd=1)
hists.summarized <- describe.simmap2(hists.fixed)

#### SAVE OUTPUTS ####
save(hists.fixed, file = "../outputs/hapauto_maps/hists.fixed.RData")
save(hists, file="../outputs/hapauto_maps/hists.RData")
save(hists.summarized, file = "../outputs/hapauto_maps/hists.summarized.RData")

