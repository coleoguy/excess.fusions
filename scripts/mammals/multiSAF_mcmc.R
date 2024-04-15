#### LOAD PACKAGES ####
library(phytools)
library(chromePlus)
library(diversitree)

#### LOAD DATA ####

dat <- read.csv("../../data/mammals/chromes/dat.csv",
                as.is=T)[,c(1,4)]
tree <- force.ultrametric(read.nexus("../../data/mammals/trees/tree.nex"),method="extend")

#### BUILD + CONSTRAIN MODEL ####

#build matrix for multiSAF model
parMat <- matrix("0",
                 nrow=max(dat$codedSCS) + 1,
                 ncol=max(dat$codedSCS) + 1)
colnames(parMat) <- 0:7
rownames(parMat) <- 0:7
mat <- parMat

for(i in 1:nrow(parMat)){
  if(i != nrow(parMat)){
    parMat[i,i+1] <- "SAF"
    mat[i,i+1] <- 1
  }
  if(i != 1){
    parMat[i,1] <- "Ro"
    mat[i,1] <- 2
  }
}

#build rate table
rate.table <- as.data.frame(matrix(, nrow(parMat) * ncol(parMat), 
                                   3))
rate.table[, 1] <- rep(as.character(1:nrow(parMat)), 
                       each = ncol(parMat))
rate.table[, 2] <- rep(as.character(1:ncol(parMat)), nrow(parMat))
rate.table[, 3] <- as.character(c(t(parMat)))
rate.table <- rate.table[rate.table[, 1] != rate.table[, 
                                                       2], ]
#make data matrix
data.matrix <- matrix(0,nrow=nrow(dat),
                        ncol=ncol(parMat))

colnames(data.matrix) <- colnames(parMat)
rownames(data.matrix) <- dat$tree.name

#fill data matrix
for(i in 1:nrow(data.matrix)){
  data.matrix[i,as.numeric(dat$codedSCS[i]) + 1] <- 1
}

#Make mkn model
model <- make.mkn(tree,
                  data.matrix,
                  ncol(data.matrix),
                  strict=F,
                  control=list(method="ode"))

#get formula for constraining
formulae <- c()
for (i in 1:nrow(rate.table)) {
  formulae[i] <- paste("q", rate.table[i, 1], rate.table[i, 
                                                         2], " ~ ", rate.table[i, 3], collapse = "", sep = "")
}

extras <- c("SAF", "Ro")

#constrain model
model.con <- constrain(model, formulae = formulae, extra = extras)

#### ESTIMATE RATES ####

#test run MCMC
temp.mcmc <- diversitree::mcmc(lik=model.con,
                                x.init=runif(2,0,1),
                                prior=make.prior.exponential(r=0.5),
                                #upper=c(100,100,100,100),
                                nsteps = 100,
                                w=1)
#extract tuning params
temp.mcmc <- temp.mcmc[-c(1:50), ]
w <- diff(sapply(temp.mcmc[2:3],
                 quantile, c(.05, .95)))


#run MCMC
model.mcmc <- diversitree::mcmc(lik=model.con,
                                x.init=runif(2,0,1),
                                prior=make.prior.exponential(r=0.5),
                                #upper=c(100,100,100,100),
                                nsteps = 500,
                                w=w)

#check for convergence
plot(x=model.mcmc$i,y=model.mcmc$p,type="l")

#### BUILD QMATRIX ####

#Extract post burn portion
model.mcmc.postburn <- model.mcmc[450:500,]

#Get mean params
params <- c(mean(model.mcmc.postburn$SAF),
            mean(model.mcmc.postburn$Ro))

names(params) <- colnames(model.mcmc[,2:3])

#Sub into matrix
parMat[parMat == "SAF"] <- params[1]
parMat[parMat == "Ro"] <- params[2]

#convert to numeric
parMat <- as.data.frame(parMat)
parMat <- sapply(parMat[,1:ncol(parMat)],as.numeric)

#calculate diagonal
diag(parMat) <- -rowSums(parMat)

#### SAVE QMATRIX AND MCMC####
save(model.mcmc,file="../../outputs/mammals/mcmc/multiSAF.RData")
write.csv(parMat,
          paste0("../../data/mammals/transition_matrix/Q_matrix_SAF.csv"),
          row.names=F,quote=F)
write.csv(mat,
          paste0("../../data/mammals/transition_matrix/transition_matrix_SAF.csv"),
          row.names=F,quote=F)






