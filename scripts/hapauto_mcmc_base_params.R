#### LOAD PACKAGES ####
library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../data/chromes/dat.csv",
                as.is=T)[,c(1,3)]
tree <- force.ultrametric(read.nexus("../data/trees/tree.nex"),method = "extend")

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
                               prior=make.prior.exponential(r=0.5),
                               #upper=c(100,100,100,100),
                               nsteps = 100,
                               w=1)

#extract tuning params
temp.mcmc <- temp.mcmc[-c(1:50), ]
w <- diff(sapply(temp.mcmc[2:3],
                 quantile, c(.05, .95)))

#run 4 replicate MCMC
model.mcmc <- list()
model.mcmc.postburn <- as.data.frame(matrix(NA,nrow=0,ncol=4))
for(i in 1:4){
  model.mcmc[[i]]<- diversitree::mcmc(lik=model.con$`likelihood function`,
                                      x.init=runif(2,0,1),
                                      prior=make.prior.exponential(r=0.5),
                                      #upper=c(100,100,100,100),
                                      nsteps = 500,
                                      w=w)
  model.mcmc.postburn <- rbind(model.mcmc.postburn,model.mcmc[[i]][451:500,])
}

#### BUILD QMATRIX ####
#Get mean params
params <- c(mean(model.mcmc.postburn$asc1),
            mean(model.mcmc.postburn$desc1))

names(params) <- colnames(model.mcmc.postburn[,2:3])

#Sub into matrix
parMat <- model.con$`parameter matrix`
mat <- model.con$`parameter matrix`
parMat[parMat == "asc1"] <- params[1]
parMat[parMat == "desc1"] <- params[2]
mat[mat == "asc1"] <- 1
mat[mat == "desc1"] <- 2

#convert to numeric
parMat <- as.data.frame(parMat)
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

#calculate diagonal
diag(parMat) <- -rowSums(parMat)

#### SAVE QMATRIX AND MCMC####
write.csv(parMat,
          paste0("../data/transition_matrix/Q_matrix_hapauto_base.csv"),
          row.names=F,quote=F)
write.csv(mat,
          paste0("../data/transition_matrix/transition_matrix_hapauto.csv"),
          row.names=F,quote=F)
save(model.mcmc,model.mcmc.postburn,file="../outputs/mcmc/hapauto_base.RData")

