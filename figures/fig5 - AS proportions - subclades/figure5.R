library(beeswarm)
library(coda)
library(viridis)
clades <- c("carnivora",
            "artiodactyla",
            "yangochiroptera",
            "anomaluromorpha_castorimorpha_myomorpha",
            "primates")
dat <- as.data.frame(matrix(NA,100, 5))
null <-  as.data.frame(matrix(NA,100, 5))
colnames(dat) <- colnames(null)<- clades
for(i in 1:5){
  #get raw data
  raw.dat <- read.csv(paste0("../../outputs/SAF_proportions/subtrees/proportions_",clades[i],"_raw.csv"))[,-1]
  dat[,i] <- raw.dat[1:100, 1]
  null[,i] <- raw.dat[101:200, 1]
}
n.hpds <- as.data.frame(matrix(NA, 5, 2))
colnames(n.hpds) <- c("low","high")
clades <- c("Carnivora",
            "Artiodactyla",
            "Yangochiroptera",
            "Rodentia",
            "Pimates")
row.names(n.hpds) <- clades
o.hpds <- n.hpds
for(i in 1:5){
  n.hpds[i, 1:2] <- HPDinterval(as.mcmc(null[,i]))[1:2]
  o.hpds[i, 1:2] <- HPDinterval(as.mcmc(dat[,i]))[1:2]
}

beeswarm(dat, method="center", corral = "wrap", col=rgb(.1,.1,.1,.2),
         pch=16, cex=.5, xaxt="n", yaxt="n", ylim=c(-0.01,0.23),ylab="Proportion SAF")
axis(side=2,at=c(0,0.05,0.10,0.15,0.20),labels=c("0.00","0.05","0.10","0.15","0.20"),cex.axis=.8,las=2)
text(x=1:5, y=rep(-0.005, 5), clades, cex=.7)
for(i in 1:5){
  if(i <5){
    lines(x=c(i,i), y=n.hpds[i,1:2],lwd=6, col=viridis(10)[7])
    lines(x=c(i,i), y=o.hpds[i,1:2],lwd=6, col=viridis(10)[2])
  }else{
    lines(x=c(i,i)+.04, y=n.hpds[i,1:2],lwd=6, col=viridis(10)[7])
    lines(x=c(i,i)-.04, y=o.hpds[i,1:2],lwd=6, col=viridis(10)[2])
  }
  
}
points(x=.5,y=.22,pch=16,col=viridis(10)[2])
points(x=.5,y=.21,pch=16,col=viridis(10)[7])
text(x=.35,y=.23,pos=4,"Highest Posterior Density",cex=.7)
text(x=.5,y=.22,pos=4,"observed",cex=.7)
text(x=.5,y=.21,pos=4,"expected",cex=.7)
# export PDF 
