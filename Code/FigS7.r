library(tidyverse)
library(openxlsx)
library(beeswarm)



#============================================================================#
#   Read in T2 network eigengenes for bronchial and tracheal datasets        #
#============================================================================#

#Bring in brown T2 network eigengenes (bronch) for the bronch dataset
eig_brown_bronch<-read.table("Data/eig_brown_bronch.txt",header=T,row.names=1)

#Bring in brown T2 network eigengenes (bronch) for the trach dataset
eig_brown_trach<-read.table("Data/eig_brown_trach.txt",header=T,row.names=1)

#Bring in purple T2 network eigengenes (trach) for the bronch dataset
eig_purple_bronch<-read.table("Data/eig_purple_bronch.txt",header=T,row.names=1)

#Bring in purple T2 network eigengenes (trach) for the trach dataset
eig_purple_trach<-read.table("Data/eig_purple_trach.txt",header=T,row.names=1)

#Combine into one table
eigs<-rbind(data.frame("eig"=eig_brown_bronch[,1],"comparison"=paste("Bronch_brown",sapply(strsplit(rownames(eig_brown_bronch),"_"),
	function(x)x[2]),sep="_"),row.names=rownames(eig_brown_bronch)),
	data.frame("eig"=eig_brown_trach[,1],"comparison"=paste("Trach_brown",gsub("CTL","BSA",sapply(strsplit(rownames(eig_brown_trach),"\\."),
	function(x)x[2])),sep="_"),row.names=rownames(eig_brown_trach)),
	data.frame("eig"=eig_purple_bronch[,1],"comparison"=paste("Bronch_purple",sapply(strsplit(rownames(eig_purple_bronch),"_"),
	function(x)x[2]),sep="_"),row.names=rownames(eig_purple_bronch)),
	data.frame("eig"=eig_purple_trach[,1],"comparison"=paste("Trach_purple",gsub("CTL","BSA",sapply(strsplit(rownames(eig_purple_trach),"\\."),
	function(x)x[2])),sep="_"),row.names=rownames(eig_purple_trach)))
eigs$comparison<-factor(eigs$comparison,levels=c("Bronch_brown_BSA","Bronch_brown_IL13","Trach_brown_BSA","Trach_brown_IL13",
	"Bronch_purple_BSA","Bronch_purple_IL13","Trach_purple_BSA","Trach_purple_IL13"))





#======================#
#   Make box plots     #
#======================#

#Make box plots
pdf("FigS7ab.pdf",width=8,height=5)
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










