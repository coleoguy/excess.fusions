#### PACKAGES ####
library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../../data/mammals/chromes/dat.csv",
                 as.is=T)[,c(1,3,4)]
tree <- force.ultrametric(tree = read.tree("../../data/mammals/Alvarez-Carretero_etal_SI/timetrees/01_step2/4705sp_mean.nwk"),
                          method = "extend")

#### PRUNE TREE ####

#Trim phylogeny and chromosome data
dat <- dat[dat$tree.name %in% tree$tip.label, ]
rownames(dat) <- 1:nrow(dat)
tree <- keep.tip(tree, dat$tree.name)

#### SAVE TREE ####
write.nexus(tree,file="../../data/mammals/trees/tree.nex")
