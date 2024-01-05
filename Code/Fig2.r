library(tidyverse)
library(plotrix)
library(beeswarm)
library(openxlsx)
library(lmerTest)


#====================================#
#              Figure 2c             #
#====================================#

#Read in the raw MMC speeds and displacements
distTab<-readRDS("Data/distTab.rds")

#Get geometric mean for each video
treatVec<-paste(distTab$wash,distTab$ko_status,distTab$treatment,distTab$donor,distTab$insert,distTab$rep,sep="_")
distTabMean<-data.frame()
for(i in 1:length(unique(treatVec))){
	currData<-distTab[which(treatVec == unique(treatVec)[i]),]
	meanData<-cbind("total_disp"=round(exp(mean(log(currData$total_disp + 1))) - 1,4),
		"speed_totalavg"=round(exp(mean(log(currData$speed_totalavg + 1))) - 1,4),
		"wash"=as.character(currData$wash[1]),"ko_status"=as.character(currData$ko_status[1]),
		"treatment"=as.character(currData$treatment[1]),
		"donor"=as.character(currData$donor[1]),"insert"=as.character(currData$insert[1]),
		"rep"=as.character(currData$rep[1]))
	distTabMean<-rbind(distTabMean,meanData)
}
distTabMean$total_disp<-as.numeric(distTabMean$total_disp)
distTabMean$speed_totalavg<-as.numeric(distTabMean$speed_totalavg)
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
distTabMean$treatment<-factor(distTabMean$treatment,levels=c("BSA","IL13"))
distTabMean$wash<-factor(distTabMean$wash,levels=c("NoWash","PBS","DTT","ATP-DTT"))
distTabMean$donor<-factor(distTabMean$donor,levels=c("T72","T73","T74"))
#Make insert separate for each donor/stimulation/insert number
distTabMean$insert<-paste(distTabMean$wash,distTabMean$ko_status,distTabMean$treatment,distTabMean$donor,distTabMean$insert,sep="_")

#Make box plots
donorColVec<-color.scale(as.numeric(distTabMean$donor),extremes=c("dodgerblue2","forestgreen","saddlebrown"),color.spec="rgb")
pdf("Fig2c.pdf",width=9,height=5)
par(bty="l",mar=c(8,5,2,3))
factorVec<-factor(paste(distTabMean$wash,distTabMean$ko_status,distTabMean$treatment,sep="_"),
	levels=c("NoWash_scrb_BSA","NoWash_scrb_IL13","NoWash_KO_BSA","NoWash_KO_IL13",
	"PBS_scrb_BSA","PBS_scrb_IL13","PBS_KO_BSA","PBS_KO_IL13",
	"DTT_scrb_BSA","DTT_scrb_IL13","DTT_KO_BSA","DTT_KO_IL13",
	"ATP-DTT_scrb_BSA","ATP-DTT_scrb_IL13","ATP-DTT_KO_BSA","ATP-DTT_KO_IL13"))
boxplot(log(distTabMean$speed_totalavg)~factorVec,las=1,col="white",outline=F,ylab="Log average speed",
	names=rep("",length(levels(factorVec))),at=c(1:4,6:9,11:14,16:19),xlab="")
beeswarm(log(distTabMean$speed_totalavg)~factorVec,pch=16,cex=0.7,
	pwcol=donorColVec,method="swarm", corral="random",add=TRUE,at=c(1:4,6:9,11:14,16:19))
text(x=c(1:4,6:9,11:14,16:19),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=rep(c("scrb, control","scrb, IL-13","ITLN1 KO, control","ITLN1 KO, IL-13"),4),
	srt=45,adj=1,xpd=T)
dev.off()


#Model speed
#Prepare table to store everything
lmerRes<-data.frame(matrix(nrow=8,ncol=8))
colnames(lmerRes)<-c("Wash","Comparison","Estimate","St. error","df","t value","p-value","interaction p-value")
lmerRes$Wash<-c("NoWash","NoWash","PBS","PBS","DTT","DTT","ATP+DTT","ATP+DTT")
lmerRes$Comparison<-c(rep(c("IL13, scrb","IL13, KO"),4))

########## None
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "NoWash"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[1,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[1,8]<-summary(currlmer)$coefficients[4,5]
#Interaction p-value = 0.477

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "NoWash"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status +	(1 | currData$donor) + (1| currData$insert))
lmerRes[2,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[2,8]<-summary(currlmer)$coefficients[4,5]


########## PBS
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "PBS"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[3,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[3,8]<-summary(currlmer)$coefficients[4,5]
#Interaction p-value = 0.0357

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "PBS"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[4,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[4,8]<-summary(currlmer)$coefficients[4,5]


########## DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "DTT"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[5,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[5,8]<-summary(currlmer)$coefficients[4,5]
#Interaction p-value = 0.109

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "DTT"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[6,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[6,8]<-summary(currlmer)$coefficients[4,5]


########## ATP + DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "ATP-DTT"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[7,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[7,8]<-summary(currlmer)$coefficients[4,5]
#Interaction p-value = 0.0375

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "ATP-DTT"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + 
	(1 | currData$donor) + (1| currData$insert))
lmerRes[8,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[8,8]<-summary(currlmer)$coefficients[4,5]












#====================================#
#              Figure 2d             #
#====================================#

#Read in the raw CBF values
CBF<-readRDS("Data/CBF.rds")

#Plot
pdf("Fig2d.pdf",height=3,width=3)
par(bty="l",mar=c(5,4,2,2))
boxplot(CBF$value~CBF$batch_treatment,las=1,at=c(1,2,4,5),col=c("white","grey"),medcol="black",names=c(rep("",4)),
	ylab="Ciliary beat frequency (Hz)",medlwd=2,cex.lab=0.9,cex.axis=0.9,outline=F,outcex=0.2,xlab="")
#beeswarm(CBF$value~CBF$batch_treatment,pwcol=colVec,method="swarm", corral="random",cex=0.04,add=TRUE,at=c(1,2,4,5))
text(x=c(1,2,4,5),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=levels(CBF$batch_treatment),srt=45,adj=1,xpd=T,cex=0.8)	
dev.off()

#Model
lmer_CBF<-lmer(CBF$value~CBF$batch * CBF$treatment + (1|CBF$donor) + (1|CBF$insert_donor_batch_treatment))

#Test for IL-13 effect in KO
IL13_effect_ko<-contest1D(lmer_CBF, rbind(c(0,0,1,1)), joint=T) #0.0232

#Test for IL-13 effect in scramble
IL13_effect_scramble<-contest1D(lmer_CBF, rbind(c(0,0,1,0)), joint=T) #1.62e-08









