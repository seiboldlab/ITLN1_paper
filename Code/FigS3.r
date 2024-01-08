library(plotrix)
library(beeswarm)
library(openxlsx)
library(lmerTest)




#====================================#
#              Figure S3a            #
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




##################### BOXPLOTS

####Displacement
donorColVec<-color.scale(as.numeric(distTabMean$donor),extremes=c("dodgerblue2","forestgreen","saddlebrown"),color.spec="rgb")
pdf("FigS2a.pdf",width=7,height=5)
par(bty="l",mar=c(8,5,2,3))
factorVec<-factor(paste(distTabMean$wash,distTabMean$ko_status,distTabMean$treatment,sep="_"),
	levels=c("NoWash_scrb_BSA","NoWash_scrb_IL13","NoWash_KO_BSA","NoWash_KO_IL13",
	"PBS_scrb_BSA","PBS_scrb_IL13","PBS_KO_BSA","PBS_KO_IL13",
	"DTT_scrb_BSA","DTT_scrb_IL13","DTT_KO_BSA","DTT_KO_IL13",
	"ATP-DTT_scrb_BSA","ATP-DTT_scrb_IL13","ATP-DTT_KO_BSA","ATP-DTT_KO_IL13"))
boxplot(log(distTabMean$total_disp)~factorVec,las=1,outline=F,ylab="Log total movement",
	names=rep("",length(levels(factorVec))),at=c(1:4,6:9,11:14,16:19),xlab="",col=c("white","grey"))
beeswarm(log(distTabMean$total_disp)~factorVec,pch=16,cex=0.7,
	pwcol=donorColVec,method="swarm", corral="random",add=TRUE,at=c(1:4,6:9,11:14,16:19))
text(x=c(1:4,6:9,11:14,16:19),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=rep(c("scrb, control","scrb, IL-13","ITLN1 KO, control","ITLN1 KO, IL-13"),5),
	srt=45,adj=1,xpd=T)
dev.off()





#################### MODELS

###Displacement

#Prepare table to store everything
lmerRes<-data.frame(matrix(nrow=8,ncol=8))
colnames(lmerRes)<-c("Wash","Comparison","Estimate","St. error","df","t value","p-value","interaction p-value")
lmerRes$Wash<-c("NoWash","NoWash","PBS","PBS","DTT","DTT","ATP+DTT","ATP+DTT")
lmerRes$Comparison<-c(rep(c("IL13, scrb","IL13, KO"),4))


########## None
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "NoWash"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[1,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[1,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -2.057, standerror = 0.2015
#Interaction p-value = 0.426

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "NoWash"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[2,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[2,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.803, standerror = 0.222
#Interaction p-value = 0.426


########## PBS
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "PBS"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[3,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[3,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -2.239, standerror = 0.231
#Interaction p-value = 0.0374

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "PBS"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[4,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[4,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.530, standerror = 0.230
#Interaction p-value = 0.0374


########## DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "DTT"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[5,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[5,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.617, standerror = 0.163
#Interaction p-value = 0.0768

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "DTT"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[6,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[6,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.191, standerror = 0.163
#Interaction p-value = 0.0768


########## ATP + DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","KO"))
currData<-distTabMean[which(distTabMean$wash == "ATP-DTT"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[7,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[7,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.622, standerror = 0.156
#Interaction p-value = 0.057

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("KO","scrb"))
currData<-distTabMean[which(distTabMean$wash == "ATP-DTT"),]
currlmer<-lmer(log(currData$total_disp)~currData$treatment * currData$ko_status + (1 | currData$donor) + (1| currData$insert))
lmerRes[8,3:7]<-summary(currlmer)$coefficients[2,1:5]
lmerRes[8,8]<-summary(currlmer)$coefficients[4,5]
# estimate = -1.187, standerror = 0.155
#Interaction p-value = 0.057












#====================================#
#              Figure S3b-c          #
#====================================#

#Read in the raw MMC speeds and displacements
distTab<-readRDS("Data/distTab_OLD.rds")

#Get geometric mean for each video
treatVec<-paste(distTab$wash,distTab$ko_status,distTab$treatment,distTab$donor,distTab$rep,sep="_")
distTabMean<-data.frame()
for(i in 1:length(unique(treatVec))){
	currData<-distTab[which(treatVec == unique(treatVec)[i]),]
	meanData<-cbind("total_disp"=round(exp(mean(log(currData$total_disp + 0.001))) - 0.001,4),
		"speed_totalavg"=round(exp(mean(log(currData$speed_totalavg + 0.00001))) - 0.00001,4),
		"pointwise_disp"=round(exp(mean(log(currData$pointwise_disp))),4),
		"speed_pointwiseavg"=round(exp(mean(log(currData$speed_pointwiseavg))),4),
		"wash"=as.character(currData$wash[1]),"ko_status"=as.character(currData$ko_status[1]),
		"treatment"=as.character(currData$treatment[1]),
		"donor"=as.character(currData$donor[1]),"rep"=as.character(currData$rep[1]))
	distTabMean<-rbind(distTabMean,meanData)
}
distTabMean$total_disp<-as.numeric(distTabMean$total_disp)
distTabMean$speed_totalavg<-as.numeric(distTabMean$speed_totalavg)
distTabMean$pointwise_disp<-as.numeric(distTabMean$pointwise_disp)
distTabMean$speed_pointwiseavg<-as.numeric(distTabMean$speed_pointwiseavg)
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
distTabMean$treatment<-factor(distTabMean$treatment,levels=c("bsa","il13"))
distTabMean$wash<-factor(distTabMean$wash,levels=c("none","pbs","dtt"))
distTabMean$donor<-factor(distTabMean$donor,levels=c("d72","d73","d74"))




#Make box plots
#Speed
donorColVec<-color.scale(as.numeric(distTabMean$donor),extremes=c("dodgerblue2","forestgreen","saddlebrown"),color.spec="rgb")
pdf("FigS2b.pdf",width=7,height=5)
par(bty="l",mar=c(8,5,2,3))
factorVec<-factor(paste(distTabMean$wash,distTabMean$ko_status,distTabMean$treatment,sep="_"),
	levels=c("none_scrb_bsa","none_scrb_il13","none_ko_bsa","none_ko_il13",
	"pbs_scrb_bsa","pbs_scrb_il13","pbs_ko_bsa","pbs_ko_il13",
	"dtt_scrb_bsa","dtt_scrb_il13","dtt_ko_bsa","dtt_ko_il13"))
boxplot(log(distTabMean$speed_totalavg)~factorVec,las=1,col=c("white","grey"),outline=F,ylab="Log total speed",
	names=rep("",length(levels(factorVec))),at=c(1:4,6:9,11:14),xlab="")
beeswarm(log(distTabMean$speed_totalavg)~factorVec,pch=16,cex=0.8,
	pwcol=donorColVec,method="swarm", corral="random",add=TRUE,at=c(1:4,6:9,11:14))
text(x=c(1:4,6:9,11:14),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=rep(c("scrb, control","scrb, IL-13","ITLN1 KO, control","ITLN1 KO, IL-13"),3),
	srt=45,adj=1,xpd=T)
dev.off()


#Test
########## None
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "none"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -2.219, standerror = 0.192
#Interaction p-value = 0.10

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "none"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.774, standerror = 0.187


########## PBS
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "pbs"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -2.498, standerror = 0.237
#Interaction p-value = 0.021

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "pbs"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.709, standerror = 0.238

########## DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "dtt"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.552, standerror = 0.147
#Interaction p-value = 0.008

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "dtt"),]
currlmer<-lmer(log(currData$speed_totalavg)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -0.985, standerror = 0.147



#Displacement
donorColVec<-color.scale(as.numeric(distTabMean$donor),extremes=c("dodgerblue2","forestgreen","saddlebrown"),color.spec="rgb")
pdf("FigS2c.pdf",width=7,height=5)
par(bty="l",mar=c(8,5,2,3))
factorVec<-factor(paste(distTabMean$wash,distTabMean$ko_status,distTabMean$treatment,sep="_"),
	levels=c("none_scrb_bsa","none_scrb_il13","none_ko_bsa","none_ko_il13",
	"pbs_scrb_bsa","pbs_scrb_il13","pbs_ko_bsa","pbs_ko_il13",
	"dtt_scrb_bsa","dtt_scrb_il13","dtt_ko_bsa","dtt_ko_il13"))
boxplot(log(distTabMean$total_disp)~factorVec,las=1,col=c("white","grey"),outline=F,ylab="Log total average movement",
	names=rep("",length(levels(factorVec))),at=c(1:4,6:9,11:14),xlab="")
beeswarm(log(distTabMean$total_disp)~factorVec,pch=16,cex=0.8,
	pwcol=donorColVec,method="swarm", corral="random",add=TRUE,at=c(1:4,6:9,11:14))
text(x=c(1:4,6:9,11:14),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=rep(c("scrb, control","scrb, IL-13","ITLN1 KO, control","ITLN1 KO, IL-13"),3),
	srt=45,adj=1,xpd=T)
dev.off()



#Test
########## None
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "none"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -2.219, standerror = 0.192
#Interaction p-value = 0.10

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "none"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.774, standerror = 0.187


########## PBS
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "pbs"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -2.498, standerror = 0.237
#Interaction p-value = 0.021

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "pbs"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.709, standerror = 0.238


########## DTT
distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("scrb","ko"))
currData<-distTabMean[which(distTabMean$wash == "dtt"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -1.552, standerror = 0.147
#Interaction p-value = 0.008

distTabMean$ko_status<-factor(distTabMean$ko_status,levels=c("ko","scrb"))
currData<-distTabMean[which(distTabMean$wash == "dtt"),]
currlmer<-lmer(log(currData$total_dsip)~currData$treatment * currData$ko_status + (1 | currData$donor))
summary(currlmer)$coefficients # estimate = -0.985, standerror = 0.147













