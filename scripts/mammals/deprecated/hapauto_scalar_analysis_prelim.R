#### PACKAGES ####
library(phytools)
library(evobiR)

#### LOAD DATA ####
dat <- read.csv("../../data/mammals/chromes/dat.csv",
                as.is=T)
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/tree.nex"),method="extend")
mat <- as.matrix(read.csv("../../data/mammals/transition_matrix/transition_matrix_hapauto.csv",
                          as.is=T,header=F))[-1,]
Qmat <- as.matrix(read.csv("../../data/mammals/transition_matrix/Q_matrix_hapauto_prelim.csv",
                           as.is=T,header=T))

#scale tree to unit rate
tree$edge.length <- tree$edge.length/max(branching.times(tree))

#### PREPARE DATA ####
colnames(mat) <- 1:49
colnames(Qmat) <- 1:49
tip.states <- dat$hapauto - 1
names(tip.states) <- dat$tree.name

#### RUN SCALING ANALYSIS ####
scaled.tree <- scaleTreeRates(tree = tree,
                              tip.states = tip.states,
                              max.ratio = 2,
                              nbins=10,
                              model = mat,
                              fixedQ = Qmat,
                              pi="fitzjohn")

#### PLOT SCALED TREE AND DETERMINE CUTOFF ####
plot.phyloscaled(scaled.tree,cex=0.05,edge.width = 0.1) #SAVE OUTPUT HERE

plot(tree,cex=0.05,edge.width = 0.1)
edgelabels(cex=0.1,frame="none")

#cut.tree <- extract.clade(tree,node=tree$edge[1120,2])
#dropClade1 <- extract.clade(tree,node=tree$edge[890,2])$tip.label
dropClade1 <- extract.clade(tree,node=tree$edge[1285,2])$tip.label
#dropClade2 <- extract.clade(tree,node=tree$edge[1027,2])$tip.label
cut.tree <- drop.tip(tree, c(dropClade1))

plot(cut.tree,cex=0.05,edge.width = 0.1)

#### SAVE SCALED TREE AND CUT TREE ####
write.nexus(cut.tree,file="../../data/mammals/trees/cut.tree.nex")
save(scaled.tree,file="../../outputs/mammals/scaled.tree.prelim.RData")


