library(tidyverse)

#Read in bronchial ITLN1 network eigengenes (brown) and tracheal ITLN1 network eigengenes (purple) 
#calculated in each of the tracheal and bronchial datasets
eig<-read.table("Data/Tracheal_Bronchial_ITLN1_network_eigengenes.txt",header=T,row.names=1,strings=F)

#Make box plots
pdf("FigSX.pdf",width=8,height=5)
par(bty="l",mar=c(8,4,2,2),mfrow=c(1,2))
boxplot(eigs$eig[grep("brown",eigs$comparison)]~droplevels(eigs$comparison[grep("brown",eigs$comparison)]),las=1,ylab="ITLN1 network eigengene",names=rep("",4),
	xlab="",col=c("white","grey"),outline=F,ylim=c(min(eigs$eig),0.3))
beeswarm(eigs$eig[grep("brown",eigs$comparison)]~droplevels(eigs$comparison[grep("brown",eigs$comparison)]),method="swarm", corral="random",add=TRUE,cex=0.8)
text(x=1:4,y=par()$usr[3]-0.06*(par()$usr[4]-par()$usr[3]),labels=c("Bronchial CTRL","Bronchial IL-13",
	"Tracheal CTRL", "Tracheal IL-13"),srt=45,adj=1,xpd=T)

boxplot(eigs$eig[grep("purple",eigs$comparison)]~droplevels(eigs$comparison[grep("purple",eigs$comparison)]),las=1,ylab="ITLN1 network eigengene",names=rep("",4),
	xlab="",col=c("white","grey"),outline=F,ylim=c(min(eigs$eig),0.3))
beeswarm(eigs$eig[grep("purple",eigs$comparison)]~droplevels(eigs$comparison[grep("purple",eigs$comparison)]),method="swarm", corral="random",add=TRUE,cex=0.8)
text(x=1:4,y=par()$usr[3]-0.06*(par()$usr[4]-par()$usr[3]),labels=c("Bronchial CTRL","Bronchial IL-13",
	"Tracheal CTRL", "Tracheal IL-13"),srt=45,adj=1,xpd=T)
dev.off()















