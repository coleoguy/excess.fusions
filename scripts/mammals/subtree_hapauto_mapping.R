#### PACKAGES ####

library(phytools)
library(doSNOW)
library(viridis)
source("../functions.R")

#### LOAD DATA ####
dat <- read.csv("../../data/mammals/chromes/dat.csv",
                as.is=T)[,c(1,3)]
mat <- as.matrix(read.csv("../../data/mammals/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header = T))
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/tree.nex"),method="extend")
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primatomorpha")

#Define edges of interest
edges <- c(181,448,796,1285,1718)

#Set up list
trees <- list()

#Loop through edges of interest
for(i in 1:length(edges)){
  
  trees[[i]] <- extract.clade(tree,
                              node=tree$edge[edges[i],2])
  
}

for(i in 1:5){
  # load in subtree and Qmatrix
  split.tree <- force.ultrametric(read.nexus(paste0("../../data/mammals/trees/subtrees/tree_",clades[i],".nex")),method="extend")
  Qmat <- as.matrix(read.csv(paste0("../../data/mammals/transition_matrix/subtree_matrices/hapauto/matrix_",clades[i],".csv"),
                            as.is=T,header = T))
  #get tree depth
  treeDepth <- max(branching.times(trees[[i]]))*100
  #calculate and print rates
  autoFus <- Qmat[2,1]/treeDepth
  autoFis <- Qmat[1,2]/treeDepth
  print(clades[i])
  print(paste0("Fusion: ", autoFus, "     Fission: ", autoFis))
  #Subset data
  split.data <- subset(dat, tree.name %in% split.tree$tip.label)
  
  #### BUILD DATA MATRIX ####
  #data.matrix
  data.matrix <- matrix(0,nrow=nrow(split.data),
                        ncol=ncol(mat))
  
  colnames(data.matrix) <- 1:49
  rownames(data.matrix) <- split.data$tree.name
  
  #fill data matrix
  for(j in 1:nrow(data.matrix)){
    data.matrix[j,as.numeric(split.data$hapauto[j]) - 1] <- 1
  }
  
  #change matri names
  rownames(mat) <- colnames(mat) <- rownames(Qmat) <- colnames(Qmat) <-1:49
  
  #### SET UP PARALLELIZATION ####
  nClust <- 100

  # Set up clusters, print intermediates to
  cl <- makeCluster(nClust, outfile = "")
  registerDoSNOW(cl)

  #### PERFORM STOCHASTIC MAPPING ####
  hists <- foreach(j=1:100,
          .verbose = T,
          .packages = c("phytools","maps","ape")) %dopar% {
            make.simmap2(tree = split.tree,
                         x = data.matrix,
                         model = mat,
                         nsim = 1,
                         Q = Qmat,
                         rejmax = 1000000,
                         rejint = 100000,
                         pi="fitzjohn",
                         monitor=T,
                         parallel=c(j,100))
          }

  stopCluster(cl)
  rm(cl)
  
  save(hists, file=paste0("../../outputs/mammals/hapauto_maps/subtrees/hists.",clades[i],".RData"))
  
  #Set class of sim.out
  class(hists) <- c("multiSimmap","multiPhylo")
  
  #### FIX STOCHASTIC MAPS ####
  split.data$sim.state <- split.data$hapauto - 1
  dat.for.fixing <- split.data[,-2]
  hists.fixed <- fix.simmap(hists,dat.for.fixing,mat)
  
  #### SUMMARIZE STOCHASTIC MAPS ####
  hists.summarized <- describe.simmap2(hists.fixed)
  
  #### SAVE OUTPUTS ####
  save(hists.summarized, file = paste0("../../outputs/mammals/hapauto_maps/subtrees/hists.",clades[i],".summarized.RData"))
  save(hists, file=paste0("../../outputs/mammals/hapauto_maps/subtrees/hists.",clades[i],".RData"))
  save(hists.fixed, file=paste0("../../outputs/mammals/hapauto_maps/subtrees/hists.",clades[i],".fixed.RData"))
}






