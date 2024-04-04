library(openxlsx)
library(tidyverse)
library(lmerTest)
library(DESeq2)
library(beeswarm)
library(edgeR)





#====================================#
#              Figure S1a            #
#====================================#

#Read in raw basolateral RNA-seq count data
secRna_353<-readRDS("Data/secRna_353.rda")

#Filter lowly expressed genes
good_genes <- rowSums(secRna_353 >= 10) >= floor(ncol(secRna_353) * .1)
secRna_353_culled<-secRna_353[good_genes,]

#Specify design
subject <- sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[1])
treatment <- sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[2]) #BSA versus IL13
status <- factor(sapply(strsplit(colnames(secRna_353_culled), split ="_"),function(x)x[3]),
	levels=c("H","A"))  #Asthmatic versus healthy
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_H","IL13_H","BSA_A","IL13_A"))
design <- data.frame(row.names = colnames(secRna_353_culled), subject = as.factor(subject), treatment = treatment,
	status = status, treatment_status = treatment_status)

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
#Add to design
design$genotype<-itln_genotypes$rs4656959

#Also make version of the dataset with GG individuals removed
design_noGG<-design %>% filter(genotype != "GG")
secRna_353_noGG<-secRna_353_culled[,rownames(design_noGG)]
secRna_353_norm_log_noGG<-secRna_353_norm_log[,rownames(design_noGG)]

#Make plots
pdf("FigS1a.pdf",height=3.6,width=5.7)
#1
par(bty="n",mfrow=c(1,2))
boxplot(simplify2array(secRna_353_norm_log["ITLN1",])~design$treatment_status,
	ylim = c(0,14),names=c("","","",""),xlab="",
	las=1,ylab="Log normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(secRna_353_norm_log["ITLN1",]~design$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
#2
boxplot(simplify2array(secRna_353_norm_log_noGG["ITLN1",])~design_noGG$treatment_status, 
	ylim = c(0,14),names=c("","","",""),xlab="",
	las=1,ylab="Log normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(secRna_353_norm_log_noGG["ITLN1",]~design_noGG$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
dev.off()

#Test IL13 effect
#IL13 effects stratified by asthma status, removing GG samples
design_noGG_temp<-design_noGG[order(design_noGG$status),]
subject_nested<-factor(c(rep(1:3,each=2),rep(1:13,each=2)))
design_noGG_temp<-cbind(design_noGG_temp,subject_nested)
design_noGG_temp<-design_noGG_temp[order(design_noGG_temp$subject),]

#Define model matrix
design_noGG_temp.mat<-model.matrix(~status + status:subject_nested + status:treatment, design_noGG_temp)

##Note that this nesting creates all possible nested variables
#Remove variables that are all zero
y <- DGEList(counts=secRna_353_noGG, group = design_noGG_temp$treatment)
y <- calcNormFactors(y)
design_noGG_temp.mat<-model.matrix(~status + status:subject_nested + status:treatment, design_noGG_temp)
idx <- which(colSums(design_noGG_temp.mat == 0) == nrow(design_noGG_temp.mat))
design_noGG_temp.mat <- design_noGG_temp.mat[,-idx]

#Perform Differential Expression
y <- DGEList(counts=secRna_353_noGG, group = design_noGG_temp$treatment)
y <- calcNormFactors(y)
y <- estimateDisp(y, design_noGG_temp.mat, robust=TRUE)
fit <- glmQLFit(y, design_noGG_temp.mat, robust=TRUE)

#Contrast for IL-13 effect in healthy
qlf <- glmQLFTest(fit, coef='statusH:treatmentIL13')
res_mucus<-topTags(qlf, n=nrow(secRna_353_noGG))
res_mucus<-res_mucus[[1]]
res_mucus["ITLN1",] #p-value = 2.03e-11, FDR = 6.14e-08, LFC = 7.46

#Contrast for IL-13 effect in asthmatic
qlf <- glmQLFTest(fit, coef='statusA:treatmentIL13')
res_mucus<-topTags(qlf, n=nrow(secRna_353_noGG))
res_mucus<-res_mucus[[1]]
res_mucus["ITLN1",] #p-value = 1.83e-17, FDR = 6.67e-14, LFC = 7.56	









#====================================#
#              Figure S1b            #
#====================================#

######### APICAL

#Bring in raw apical secretome counts
apical<-read.table("Data/secretome_count_matrix.txt")

#Bring in the genotypes
itln_genotypes<-read.table("Data/secretome_metadata.txt",sep="\t",header=T,stringsAsFactors=F)

#Cull the genotype list to match the donors we have
itln_genotypes<-itln_genotypes[which(gsub("HBEC","X",itln_genotypes$Sample) %in% 
	sapply(strsplit(colnames(apical),"_"),function(x)x[1])),]
itln_genotypes<-itln_genotypes[order(itln_genotypes$Sample),c(1,6)]
itln_genotypes<-itln_genotypes[rep(seq_len(nrow(itln_genotypes)), each=2),]

#Specify design
subject <- sapply(strsplit(colnames(apical), split ="_"),function(x)x[1])
treatment <- factor(sapply(strsplit(colnames(apical), split ="_"),function(x)x[2]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(sapply(strsplit(colnames(apical), split ="_"),function(x)x[4]),levels=c("H","A")) #BSA versus IL13
genotype <- factor(itln_genotypes$rs4656959,levels=c("AA","AG","GG"))
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_H","IL13_H","BSA_A","IL13_A"))
design <- data.frame(row.names = colnames(apical), subject = subject, treatment = treatment,
	status = status, genotype = genotype, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = apical,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
apical_norm<-counts(dds_sizeFactors, normalized=T)

#Also make versions with GG individuals removed
design_noGG<-design %>% filter(genotype != "GG")
apical_noGG<-apical[,rownames(design_noGG)]
apical_norm_noGG<-apical_norm[,rownames(design_noGG)]

#Make plots
pdf("FigS1b.pdf",height=3.6,width=5.7)
#1
par(bty="n",mfrow=c(1,2))
boxplot(simplify2array(apical_norm["ITLN1_ITLN1",])~design$treatment_status,
	ylim = c(0,14),names=c("","","",""),xlab="",
	las=1,ylab="Normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(apical_norm["ITLN1_ITLN1",]~design$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
#2
boxplot(simplify2array(apical_norm_noGG["ITLN1_ITLN1",])~design_noGG$treatment_status, 
	ylim = c(0,14),names=c("","","",""),xlab="",
	las=1,ylab="Normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(apical_norm_noGG["ITLN1_ITLN1",]~design_noGG$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
dev.off()

#Test
design_justGG<-design_noGG[which(design_noGG$status == "A"),]
d1_justGG<-apical_noGG[,rownames(design_justGG)]
y <- DGEList(counts=d1_justGG, group = design_justGG$treatment)
y <- calcNormFactors(y)
design_noGG_temp.mat<-model.matrix(~subject + treatment, design_justGG)
y <- estimateDisp(y, design_noGG_temp.mat, robust=TRUE)
fit <- glmQLFit(y, design_noGG_temp.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='treatmentIL13')
res_mucus_IL13<-topTags(qlf, n=nrow(d1_justGG))
res_mucus_IL13<-res_mucus_IL13[[1]]
res_mucus_IL13["ITLN1_ITLN1",] #p-value = 4.45E-22, FDR = 7.68E-20, LFC = 6.15












#====================================#
#              Figure S1c            #
#====================================#

######### MUCUS

#Bring in the apical mucus secretome counts
mucus<-read.table("Data/secretome_count_matrix_newMucus.txt")
colnames(mucus)<-sapply(strsplit(colnames(mucus), split ="_"),function(x)paste(x[2],x[3],x[4],x[5],sep="_"))
mucus<-mucus[,order(colnames(mucus))]

#Bring in the genotypes
itln_genotypes<-read.table("Data/secretome_metadata.txt",sep="\t",header=T,stringsAsFactors=F)

#Cull the genotype list to match the donors we have
itln_genotypes<-itln_genotypes[which(gsub("HBEC","",itln_genotypes$Sample) %in% 
	sapply(strsplit(colnames(mucus),"_"),function(x)x[1])),]
itln_genotypes<-itln_genotypes[order(itln_genotypes$Sample),c(1,6)]
itln_genotypes<-itln_genotypes[rep(seq_len(nrow(itln_genotypes)), each=2),]

#Specify design
subject <- sapply(strsplit(colnames(mucus), split ="_"),function(x)x[1])
treatment <- factor(sapply(strsplit(colnames(mucus), split ="_"),function(x)x[4]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(sapply(strsplit(colnames(mucus), split ="_"),function(x)x[3]),levels=c("H","A")) #BSA versus IL13
genotype <- factor(itln_genotypes$rs4656959,levels=c("AA","AG","GG"))
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_H","IL13_H","BSA_A","IL13_A"))
design <- data.frame(row.names = colnames(mucus), subject = subject, treatment = treatment,
	status = status, genotype = genotype, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = mucus,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
mucus_norm<-counts(dds_sizeFactors, normalized=T)

#Also make versions with GG individuals removed
design_noGG<-design %>% filter(genotype != "GG")
mucus_noGG<-mucus[,rownames(design_noGG)]
mucus_norm_noGG<-mucus_norm[,rownames(design_noGG)]

#Make plots
pdf("FigS1c.pdf",height=3.6,width=5.7)
#1
par(bty="n",mfrow=c(1,2))
boxplot(simplify2array(mucus_norm["ITLN1",])~design$treatment_status,
	ylim = c(0,8),names=c("","","",""),xlab="",
	las=1,ylab="Normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(mucus_norm["ITLN1",]~design$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
#2
boxplot(simplify2array(mucus_norm_noGG["ITLN1",])~design_noGG$treatment_status, 
	ylim = c(0,8),names=c("","","",""),xlab="",
	las=1,ylab="Normalized expression",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.9,outline=F,cex.out=0.5)
beeswarm(mucus_norm_noGG["ITLN1",]~design_noGG$treatment_status, method="swarm", corral="random",ylab="",xlab="",
	pch=16,cex=0.9,add=T)
text(x=1.3:4.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=levels(design$treatment_status),srt=45,adj=1,xpd=T,cex=0.8)
dev.off()


#Test
#IL13 effects stratified by asthma status, removing GG samples
design_noGG_temp<-design_noGG[order(design_noGG$status),]
subject_nested<-factor(c(rep(1:2,each=2),rep(1:4,each=2)))
design_noGG_temp<-cbind(design_noGG_temp,subject_nested)
design_noGG_temp$treatment <- factor(design_noGG_temp$treatment,levels=c("BSA","IL13"))
design_noGG_temp<-design_noGG_temp[order(design_noGG_temp$subject),]

#Define model matrix
design_noGG_temp.mat<-model.matrix(~status + status:subject_nested + status:treatment, design_noGG_temp)

##Note that this nesting creates all possible nested variables
#Remove variables that are all zero
y <- DGEList(counts=mucus_noGG, group = design_noGG_temp$treatment)
y <- calcNormFactors(y)
design_noGG_temp.mat<-model.matrix(~status + status:subject_nested + status:treatment, design_noGG_temp)
idx <- which(colSums(design_noGG_temp.mat == 0) == nrow(design_noGG_temp.mat))
design_noGG_temp.mat <- design_noGG_temp.mat[,-idx]

#Perform Differential Expression
y <- DGEList(counts=mucus_noGG, group = design_noGG_temp$treatment)
y <- calcNormFactors(y)
y <- estimateDisp(y, design_noGG_temp.mat, robust=TRUE)
fit <- glmQLFit(y, design_noGG_temp.mat, robust=TRUE)

#Contrast for IL-13 effect in healthy
qlf <- glmQLFTest(fit, coef='statusH:treatmentIL13')
res_mucus<-topTags(qlf, n=nrow(mucus_noGG))
res_mucus<-res_mucus[[1]]
res_mucus["ITLN1",]
#   logFC   logCPM        F      PValue       FDR
#5.619198 8.878016 161.9776 0.005856378 0.8651172

#Contrast for IL-13 effect in asthmatic
qlf <- glmQLFTest(fit, coef='statusA:treatmentIL13')
res_mucus<-topTags(qlf, n=nrow(mucus_noGG))
res_mucus<-res_mucus[[1]]
res_mucus["ITLN1",]
#   logFC   logCPM        F      PValue       FDR
#5.789402 8.878016 391.3461 0.002412796 0.2491384	









#====================================#
#              Figure S1d            #
#====================================#

#Bring in cell fractions from CIBERSORTx
fracts<-read.table("Data/CIBERSORTx_results_SmodeBatchCorrection_0.5threshold.txt",sep="\t",header=T,strings=F)

#Bring in brown module eigengenes
eig<-read.table("Data/WGCNA.eigengenes.txt",header=T,strings=F)

#Correlate ITLN1 network expression with fractions
correlations<-sapply(fracts[,c(2:7)],function(x)cor(x,eig$MEbrown)) %>% sort(decreasing=T)
ci_lower<-sapply(fracts[,c(2:7)],function(x)cor.test(x,eig$MEbrown)$conf.int[1]) %>% sort(decreasing=T)
ci_upper<-sapply(fracts[,c(2:7)],function(x)cor.test(x,eig$MEbrown)$conf.int[2]) %>% sort(decreasing=T)

#Now plot
pdf("FigS1d.pdf",width=5,height=5.5)
cellTypeLabelsForPlot<-names(correlations)
par(mar=c(12,7,2,2),bty="l")
plot(x=1:length(correlations),y=correlations,ylim=c(-1,1),pch=16,las=1,xlab="",xaxt='n',ylab="Pearson correlation")
for(i in 1:length(correlations)){
	lines(x=c(i,i),y=c(ci_lower[i],ci_upper[i]))
}
abline(h=0,lty=2)
axis(1,at=1:length(correlations),label=rep("",length(correlations)))
text(x=1:length(correlations),y=par()$usr[3]-0.07*(par()$usr[4]-par()$usr[3]),labels=cellTypeLabelsForPlot,
		srt=45,adj=1,xpd=T)
dev.off()







