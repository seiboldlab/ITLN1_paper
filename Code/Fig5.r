library(plotrix)
library(beeswarm)
library(openxlsx)
library(lmerTest)
library(ggbeeswarm)
library(cowplot)
library(DESeq2)
library(edgeR)
library(tidyverse)



#====================================#
#              Figure 5a             #
#====================================#

#Read in raw basolateral RNA-seq count data
secRna_353<-readRDS("Data/secRna_353.rda")

#Filter lowly expressed genes
good_genes <- rowSums(secRna_353 >= 10) >= floor(ncol(secRna_353) * .1)
secRna_353_culled<-secRna_353[good_genes,]

#Specify design
subject <- sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[1])
treatment <- sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[2]) #BSA versus IL13
status <- sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[3])  #Asthmatic versus healthy
design <- data.frame(row.names = colnames(secRna_353_culled), subject = as.factor(subject), treatment = treatment,
	status = status)

#Normalize data
dds_Rna <- DESeqDataSetFromMatrix(
  countData = secRna_353_culled,
  colData = design,
  design = ~subject + treatment)
dds_Rna$treatment <- relevel(dds_Rna$treatment,"BSA")
dds_Rna_sizeFactors<-estimateSizeFactors(dds_Rna)
secRna_353_norm<-counts(dds_Rna_sizeFactors, normalized=T)
secRna_353_norm_log<-log2(secRna_353_norm + 1)

#Bring in the genotypes
itln_genotypes<-read.table("Data/secretome_metadata.txt",sep="\t",header=T,stringsAsFactors=F)

#Cull the genotype list to match the donors we have
itln_genotypes<-itln_genotypes[which(itln_genotypes$Sample %in% sapply(strsplit(colnames(secRna_353_norm_log),"_"),function(x)x[1])),]
itln_genotypes<-itln_genotypes[order(itln_genotypes$Sample),]
itln_genotypes<-itln_genotypes[rep(1:nrow(itln_genotypes),each=2),] 
itln_genotypes$rs4656959<-factor(itln_genotypes$rs4656959,levels=c("AA","AG","GG"))

#Re-set design matrix
subject <- sapply(strsplit(colnames(secRna_353_norm_log), split ="_"),function(x)x[1])
treatment <- sapply(strsplit(colnames(secRna_353_norm_log), split ="_"),function(x)x[2]) #BSA versus IL13
status <- itln_genotypes$rs4656959  #AA vs GG
design <- data.frame(row.names = colnames(secRna_353_norm_log), subject = as.factor(subject),treatment=treatment, status = status)

#Make box plots for all the data
pdf("Fig5a.pdf",height=3,width=3.5)
par(bty="l")
boxplot(secRna_353_norm_log["ITLN1",]~factor(paste(design$treatment,design$status,sep="_"),
	levels=c("BSA_AA","BSA_AG","BSA_GG","IL13_AA","IL13_AG","IL13_GG")),las=1,
	ylab="Log normalizedexpression",cex.main=1,medlwd=2,cex.lab=0.8,cex.axis=0.8,cex=0.5,
	outline=F,col=c(rep("white",3),rep("grey",3)),medcol="black",at=c(1:3,5:7),
	xlab="rs4656959 genotype",names=rep(c(""),6),ylim=c(min(secRna_353_norm_log["ITLN1",]),14))
beeswarm(secRna_353_norm_log["ITLN1",]~factor(paste(design$treatment,design$status,sep="_"),
	levels=c("BSA_AA","BSA_AG","BSA_GG","IL13_AA","IL13_AG","IL13_GG")),
	method="swarm", corral="random",add=TRUE,cex=0.5,corralWidth=3,pch=16,at=c(1:3,5:7))
text(x=c(1.4:3.4,5.4:7.4),y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=rep(c("A/A","A/G","G/G"),2),srt=0,adj=1,xpd=T,cex=0.7)
dev.off()

# Differential expression analysis
#IGG and AA, BSA
design_temp<-design[which(design$treatment == "BSA" & itln_genotypes$rs4656959 != "AG"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=secRna_353_culled[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(secRna_353_culled[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#logFC    logCPM        F      PValue       FDR
#ITLN1 -2.79462 0.7018145 11.47249 0.003749978 0.9986488

#GG and AG, BSA 
design_temp<-design[which(design$treatment == "BSA" & itln_genotypes$rs4656959 != "AA"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=secRna_353_culled[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(secRna_353_culled[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#   logFC    logCPM        F     PValue      FDR
#2.780744 0.5737707 7.262262 0.01856075 0.997537

#GG and AA, IL13
design_temp<-design[which(design$treatment == "IL13" & itln_genotypes$rs4656959 != "AG"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=secRna_353_culled[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(secRna_353_culled[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#    logFC   logCPM        F       PValue       FDR
#-7.034811 8.185806 30.69541 5.916891e-05 0.5110375

#GG and AG, IL13
design_temp<-design[which(design$treatment == "IL13" & itln_genotypes$rs4656959 != "AA"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=secRna_353_culled[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(secRna_353_culled[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#    logFC   logCPM        F       PValue        FDR
#-6.344925 7.420357 65.76217 1.540819e-06 0.02334649








#====================================#
#              Figure 5b             #
#====================================#

######### APICAL

#Bring in raw apical secretome counts
apical<-read.table("Data/secretome_count_matrix.txt")

#Next, bring in the genotypes
itln_genotypes<-read.table("Data/secretome_metadata.txt",sep="\t",header=T,stringsAsFactors=F)

#Cull the genotypelist to match the donors we have
itln_genotypes<-itln_genotypes[which(gsub("HBEC","X",itln_genotypes$Sample) %in% 
	sapply(strsplit(colnames(apical),"_"),function(x)x[1])),]
itln_genotypes<-itln_genotypes[order(itln_genotypes$Sample),c(1,6)]
itln_genotypes<-itln_genotypes[rep(seq_len(nrow(itln_genotypes)), each=2),]

#Specify design
subject <- sapply(strsplit(colnames(apical), split ="_"),function(x)x[1])
treatment <- factor(sapply(strsplit(colnames(apical), split ="_"),function(x)x[2]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(itln_genotypes$rs4656959,levels=c("AA","AG","GG"))
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_AA","IL13_AA","BSA_AG","IL13_AG","BSA_GG","IL13_GG"))
design <- data.frame(row.names = colnames(apical), subject = subject, treatment = treatment,
	status = status, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = apical,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
apical_norm<-counts(dds_sizeFactors, normalized=T)

#Box plots
pdf("Fig5b1.pdf",width=2.5,height=3.1)
par(bty="l")
boxplot(apical_norm["ITLN1_ITLN1",which(treatment=="IL13")]~status[which(treatment=="IL13")],
	las=1,ylab="Normalized abundance",cex.main=1,medlwd=2,cex.lab=0.8,cex.axis=0.8,
	outline=F,col="grey",xlab="rs4656959 genotype",
	ylim=c(min(apical_norm["ITLN1_ITLN1",which(treatment=="IL13")]),13))
beeswarm(apical_norm["ITLN1_ITLN1",which(treatment=="IL13")]~status[which(treatment=="IL13")],
	col="black",method="swarm", corral="random",ylab="",xlab="",pch=16,cex=0.7,add=T)
dev.off()


#Differential expression analysis
design_temp<-design[which(design$treatment == "IL13"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=apical[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(apical[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#    logFC   logCPM        F       PValue          FDR
#-6.138874 12.72661 27.45055 1.654912e-07 0.0001143544





######### MUCUS

#Bring in the apical mucus secretome counts
mucus<-read.table("Data/secretome_count_matrix_newMucus.txt")
colnames(mucus)<-sapply(strsplit(colnames(mucus), split ="_"),function(x)paste(x[2],x[3],x[4],x[5],sep="_"))
mucus<-mucus[,order(colnames(mucus))]

#Next, bring in the genotypes
itln_genotypes<-read.table("Data/secretome_metadata.txt",sep="\t",header=T,stringsAsFactors=F)

#Cull the genotype list to match the donors we have
itln_genotypes<-itln_genotypes[which(gsub("HBEC","",itln_genotypes$Sample) %in% 
	sapply(strsplit(colnames(mucus),"_"),function(x)x[1])),]
itln_genotypes<-itln_genotypes[order(itln_genotypes$Sample),c(1,6)]
itln_genotypes<-itln_genotypes[rep(seq_len(nrow(itln_genotypes)), each=2),]

#Specify design
subject <- sapply(strsplit(colnames(mucus), split ="_"),function(x)x[1])
treatment <- factor(sapply(strsplit(colnames(mucus), split ="_"),function(x)x[4]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(itln_genotypes$rs4656959,levels=c("AA","AG","GG"))
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_AA","IL13_AA","BSA_AG","IL13_AG","BSA_GG","IL13_GG"))
design <- data.frame(row.names = colnames(mucus), subject = subject, treatment = treatment,
	status = status, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = mucus,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
mucus_norm<-counts(dds_sizeFactors, normalized=T)

#Box plots
pdf("Fig5b1.pdf",width=2.5,height=3.1)
par(bty="l")
boxplot(mucus_norm["ITLN1",which(treatment=="IL13")]~status[which(treatment=="IL13")],
	las=1,ylab="Normalized abundance",cex.main=1,medlwd=2,cex.lab=0.8,cex.axis=0.8,
	outline=F,col="grey",xlab="rs4656959 genotype",
	ylim=c(min(mucus_norm["ITLN1",which(treatment=="IL13")]),9))
beeswarm(mucus_norm["ITLN1",which(treatment=="IL13")]~status[which(treatment=="IL13")],
	col="black",method="swarm", corral="random",ylab="",xlab="",pch=16,cex=0.7,add=T)
dev.off()

#Differential expression analysis
##GG vs AA
design_temp<-design[which(design$treatment == "IL13"),]
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=mucus[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(mucus[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#    logFC   logCPM        F       PValue        FDR
#-5.833219 9.117397 127.1691 1.409585e-05 0.03748086

##GG vs AG+AA
design_temp<-design[which(design$treatment == "IL13"),]
design_temp$status<-plyr::mapvalues(design_temp$status,from="AG",to="AA")
design_temp$status<-droplevels(design_temp$status)
y <- DGEList(counts=mucus[,rownames(design_temp)], group = design_temp$status)
y <- calcNormFactors(y)
design.mat<-model.matrix(~status, design_temp)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
##GG vs AA
qlf <- glmQLFTest(fit, coef='statusGG')
res<-topTags(qlf, n=nrow(mucus[,rownames(design_temp)]))
res<-res[[1]]
res["ITLN1",]
#    logFC   logCPM        F       PValue         FDR
#-5.752943 9.117397 145.9287 2.543639e-06 0.006763535











#====================================#
#              Figure 5c             #
#====================================#

#Bring in tracheal expression of ITLN1
exptab<-read.table("Data/ITLN1_exp.txt",sep="\t",header=T,stringsAsFactors=F)[,1:3]

# AA
#Test for difference 
exptab_AA<-exptab[grep("AA",exptab$comparison),]
currlmer<-lmer(log10(exptab_AA$exp)~exptab_AA$comparison + (1 | exptab_AA$Donor))
summary(currlmer)$coefficients #2.06085e-05
# GG
#Test for difference 
exptab_GG<-exptab[grep("GG",exptab$comparison),]
currlmer<-lmer(log10(exptab_GG$exp)~exptab_GG$comparison + (1 | exptab_GG$Donor))
summary(currlmer)$coefficients #4.115948e-02




#Bring in tracheal expression of MUC5AC
exptab<-read.table("Data/MUC5AC_exp.txt",sep="\t",header=T,stringsAsFactors=F)[,1:3]

# AA
#Test for difference 
exptab_AA<-exptab[grep("AA",exptab$comparison),]
currlmer<-lmer(log10(exptab_AA$exp)~exptab_AA$comparison + (1 | exptab_AA$Donor))
summary(currlmer)$coefficients #0.0002799360
# GG
#Test for difference 
exptab_GG<-exptab[grep("GG",exptab$comparison),]
currlmer<-lmer(log10(exptab_GG$exp)~exptab_GG$comparison + (1 | exptab_GG$Donor))
summary(currlmer)$coefficients #0.001123401





#Box plots
exptab_I<-read.table("Data/ITLN1_exp.txt",sep="\t",header=T,stringsAsFactors=F)[,1:3]
exptab_M<-read.table("Data/MUC5AC_exp.txt",sep="\t",header=T,stringsAsFactors=F)[1:20,1:3]
exptab<-rbind(data.frame(exptab_M,"gene"="MUC5AC"),data.frame(exptab_I,"gene"="ITLN1"))
treat<-factor(paste(exptab$comparison,exptab$gene,sep="_"),
	levels=c(unique(paste(exptab$comparison,exptab$gene,sep="_"))[c(1:2,5:6,3:4,7:8)]))
pdf("Fig5c.pdf",width=5,height=3.5)
par(bty="l",mar=c(6,5,2,2))
boxplot(log10(exptab$exp)~treat,las=1,col=c("white","gray"),cex=0.7,ylim=c(-1.7,1.5),outline=F,medlwd=2,
	yaxt="n",names=rep("",8),xlab="",at=c(1,2,4,5,7,8,10,11),ylab="Log normalized expression")
axis(2,at=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),las=1)
beeswarm(log10(exptab$exp)~treat,method="swarm", corral="random",add=TRUE,cex=0.7,
	corralWidth=3,pch=16,at=c(1,2,4,5,7,8,10,11))
text(x=c(1.2:2.4,4.2,5.4,7.2,8.4,10.2,11.4),y=par()$usr[3]-0.12*(par()$usr[4]-par()$usr[3]),
	labels=levels(treat),srt=45,adj=1,xpd=T,cex=0.7)
dev.off()

















