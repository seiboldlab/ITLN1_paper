library(plotrix)
library(beeswarm)
library(openxlsx)
library(lmerTest)



#====================================#
#              Figure S4a            #
#====================================#

#Read in the raw CBF values
CBF<-readRDS("Data/CBF.rds")

#Create factor for stratification
CBF$batch_treatment_donor_insert<-paste(CBF$batch,CBF$treatment,CBF$donor,CBF$insert,sep="_")
CBF$batch_treatment_donor_insert<-factor(CBF$batch_treatment_donor_insert,
	levels=c("scramble_BSA_T72_I1","scramble_BSA_T72_I2","scramble_BSA_T72_I3","scramble_IL13_T72_I1","scramble_IL13_T72_I2","scramble_IL13_T72_I3",
	"scramble_BSA_T73_I1","scramble_BSA_T73_I2","scramble_BSA_T73_I3","scramble_IL13_T73_I1","scramble_IL13_T73_I2","scramble_IL13_T73_I3",
	"scramble_BSA_T74_I1","scramble_BSA_T74_I2","scramble_BSA_T74_I3","scramble_IL13_T74_I1","scramble_IL13_T74_I2","scramble_IL13_T74_I3",	
	"KO_BSA_T72_I2","KO_BSA_T72_I3","KO_IL13_T72_I2","KO_IL13_T72_I3",
	"KO_BSA_T73_I1","KO_BSA_T73_I2","KO_BSA_T73_I3","KO_IL13_T73_I1","KO_IL13_T73_I2","KO_IL13_T73_I3",
	"KO_BSA_T74_I1","KO_BSA_T74_I2","KO_BSA_T74_I3","KO_IL13_T74_I1","KO_IL13_T74_I2","KO_IL13_T74_I3"))

#Make box plot
pdf("FigS4a.pdf",height=4.5,width=9)
par(bty="l")
boxplot(CBF$value~CBF$batch_treatment_donor_insert,las=1,
	at=c(1,2,3,4,5,6,9,10,11,12,13,14,17,18,19,20,21,22,25,26,27,28,31,32,33,34,35,36,39,40,41,42,43,44),
	col=c(rep(c(rep("white",3),rep("grey",3)),3),rep("white",2),rep("grey",2),rep(c(rep("white",3),rep("grey",3)),2)),
	medcol="black",names=c(rep("",34)),ylab="Ciliary beat frequency (Hz)",medlwd=2,cex.lab=0.9,
	cex.axis=0.9,outline=F,outcex=0.3,xlab="",ylim=c(0,30))
dev.off()


#Run interaction tests for each donor separately
## IL-13 effects in scramble
currData<-CBF[which(CBF$donor == "T72"),]
lmer_CBF<-lmer(currData$value~currData$batch * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=0.000101

currData<-CBF[which(CBF$donor == "T73"),]
lmer_CBF<-lmer(currData$value~currData$batch * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=0.0978

currData<-CBF[which(CBF$donor == "T74"),]
lmer_CBF<-lmer(currData$value~currData$batch * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=2.81e-06


#IL-13 effects in KO
currData<-CBF[which(CBF$donor == "T72"),]
lmer_CBF<-lmer(currData$value~factor(currData$batch,levels=c("KO","scramble")) * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=0.0299)

currData<-CBF[which(CBF$donor == "T73"),]
lmer_CBF<-lmer(currData$value~factor(currData$batch,levels=c("KO","scramble")) * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=0.121)

currData<-CBF[which(CBF$donor == "T74"),]
lmer_CBF<-lmer(currData$value~factor(currData$batch,levels=c("KO","scramble")) * currData$treatment + (1|currData$insert_donor_batch_treatment))
summary(lmer_CBF) # p=0.0956)











#====================================#
#              Figure S4b            #
#====================================#

#Bring in CBF data
CBF<-readRDS("Data/CBF_OLD.rds")

#Do box plots not broken out
#Add column combining treatment and batch
CBF$batch_treatment<-paste(CBF$batch,CBF$treatment,sep="_")
CBF$batch_treatment<-factor(CBF$batch_treatment,levels=c("scramble_BSA","scramble_IL13","KO_BSA","KO_IL13"))

#Plot
pdf("FigS4b.pdf",height=3.5,width=2.6)
par(bty="l")
boxplot(CBF$value~CBF$batch_treatment,las=1,at=c(1,2,3,4),col=c("white","grey"),medcol="black",names=c(rep("",4)),
	ylab="Ciliary beat frequency (Hz)",medlwd=2,cex.lab=0.9,cex.axis=0.9,outline=T,outcex=0.4,
	ylim=range(CBF$value))
dev.off()

#Test
CBF$batch<-factor(CBF$batch,levels=c("scramble","KO"))
lmer_CBF<-lmer(CBF$value~CBF$batch * CBF$treatment + (1 | CBF$donor)) 

#Test for IL-13 effect in KO
IL13_effect_ko<-contest1D(lmer_CBF, rbind(c(0,0,1,1)), joint=T) #0.9327

#Test for IL-13 effect in scramble
IL13_effect_scramble<-contest1D(lmer_CBF, rbind(c(0,0,1,0)), joint=T) #9.40e-26

#Test for interaction 
interaction_test<-contest1D(lmer_CBF, rbind(c(0,0,0,1)), joint=T) #2.60e-18











#====================================#
#              Figure S4c            #
#====================================#

#Bring in CBF data
CBF<-readRDS("Data/CBF_OLD.rds")

#Do box plots
#Add column combining treatment and batch
CBF$batch_treatment_donor_insert<-paste(CBF$batch,CBF$treatment,CBF$donor,CBF$insert,sep="_")
CBF$batch_treatment_donor_insert<-factor(CBF$batch_treatment_donor_insert,
	levels=c("scramble_BSA_T72_i2a","scramble_IL13_T72_i1a","scramble_BSA_T73_i2a","scramble_IL13_T73_i1a",
	"KO_BSA_T72_i2b","KO_BSA_T72_i4","KO_IL13_T72_i1b","KO_IL13_T72_i3",
	"KO_BSA_T73_i2b","KO_BSA_T73_i4","KO_IL13_T73_i1b","KO_IL13_T73_i3"))

pdf("FigS4c.pdf",height=3.5,width=6)
par(bty="l",mar=c(7,4,2,2))
boxplot(CBF$value~CBF$batch_treatment_donor_insert,las=1,at=c(1,2,4,5,7,8,9,10,12,13,14,15),
	col=c("white","grey","white","grey","white","white","grey","grey","white","white","grey","grey"),
	medcol="black",names=c(rep("",12)),ylab="Ciliary beat frequency (Hz)",medlwd=2,cex.lab=0.9,
	cex.axis=0.9,outline=T,outcex=0.3,xlab="",ylim=c(min(CBF$value),max(CBF$value)))
text(x=c(1,2,4,5,7,8,9,10,12,13,14,15),y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
	labels=c("insert 1","insert 1","insert 1",
	"insert 1","insert 1","insert 2","insert 1",
	"insert 2", "insert 1","insert 2","insert 1",
	"insert 2"),srt=45,adj=1,xpd=T,cex=0.8)	
dev.off()


#test for each donor and ko/scramble factor separately
currData<-CBF[which(CBF$donor == "T72" & CBF$batch == "scramble"),]
lmer_CBF<-lm(currData$value~currData$treatment)
summary(lmer_CBF) #0.062

currData<-CBF[which(CBF$donor == "T73" & CBF$batch == "scramble"),]
lmer_CBF<-lm(currData$value~currData$treatment)
summary(lmer_CBF) #6.29e-41

currData<-CBF[which(CBF$donor == "T72" & CBF$batch == "KO"),]
lmer_CBF<-lmer(currData$value~currData$treatment + (1 | currData$insertDonor))
summary(lmer_CBF) #0.309

currData<-CBF[which(CBF$donor == "T73" & CBF$batch == "KO"),]
lmer_CBF<-lmer(currData$value~currData$treatment + (1 | currData$insertDonor))
summary(lmer_CBF) #0.279













