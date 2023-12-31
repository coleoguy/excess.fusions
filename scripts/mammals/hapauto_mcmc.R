#### LOAD PACKAGES ####
library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../../data/mammals/chromes/dat.csv",
                as.is=T)[,c(1,3)]
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/cut.tree.nex"))

#subset dat to only include tips in cut tree
dat <- dat[which(dat$tree.name %in% tree$tip.label),]

#### MODEL ####

#Convert to data matrix
data.matrix <- datatoMatrix(dat,
                            c(2,50),
                            hyper = F)

#Make mkn model
model <- make.mkn(tree,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#constrain model
model.con <- constrainMkn(data.matrix,
                          model,
                          hyper = F,
                          polyploidy = F,
                          verbose = T,
                          constrain = list(drop.poly=T,
                                           drop.demi=T))

#test run MCMC
temp.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                x.init=runif(2,0,1),
                                prior=make.prior.exponential(r=10),
                                #upper=c(100,100,100,100),
                                nsteps = 100,
                                w=1)

#extract tuning params
temp.mcmc <- temp.mcmc[-c(1:50), ]
w <- diff(sapply(temp.mcmc[2:3],
                 quantile, c(.05, .95)))

#run MCMC
model.mcmc <- diversitree::mcmc(lik=model.con$`likelihood function`,
                                x.init=runif(2,0,1),
                                prior=make.prior.exponential(r=10),
                                #upper=c(100,100,100,100),
                                nsteps = 500,
                                w=w)

#### BUILD QMATRIX ####
model.mcmc.postburn <- model.mcmc[450:500,]

#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1))

names(params) <- colnames(model.mcmc[,2:3])

#Sub into matrix
parMat <- model.con$`parameter matrix`
parMat[parMat == "asc1"] <- params[1]
parMat[parMat == "desc1"] <- params[2]

#convert to numeric
parMat <- as.data.frame(parMat)
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

#calculate diagonal
diag(parMat) <- -rowSums(parMat)

#### SAVE QMATRIX ####
write.csv(parMat,
          paste0("../../data/mammals/transition_matrix/Q_matrix_hapauto.csv"),
          row.names=F,quote=F)

