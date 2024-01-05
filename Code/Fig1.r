library(tidyverse)
library(WGCNA)
library(DESeq2)
library(edgeR)
library(beeswarm)
library(openxlsx)

source("Helper_functions.r")


#====================================#
#              Figure 1b             #
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

# Make box plot of ITLN1,stratified by IL-13 vs BSA
pdf("Fig1b.pdf",height=3.5,width=2.5)
boxplot(secRna_353_norm_log["ITLN1",]~design$treatment,col=c("white","grey"),
	las=1,ylab="Log normalized expression",medlwd=2,cex.lab=0.8,cex.axis=0.8,
	outline=F,ylim=c(min(secRna_353_norm_log["ITLN1",]),14),xlab="",names=c("",""))
beeswarm(secRna_353_norm_log["ITLN1",]~design$treatment,method="swarm", corral="random",
	add=TRUE,cex=0.5,corralWidth=3,pch=16)
text(x=1.3:2.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=c("control","IL-13"),srt=0,adj=1,xpd=T,cex=0.8)
dev.off()

# Differential expression analysis
y <- DGEList(counts=secRna_353, group = design$treatment)
y <- calcNormFactors(y)
design.mat<-model.matrix(~subject + treatment, design)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='treatmentIL13')
res_Rna353<-topTags(qlf, n=nrow(secRna_353))
res_Rna353<-res_Rna353[[1]]
res_Rna353["ITLN1","PValue"] #3.63e-19
res_Rna353["ITLN1","FDR"] #5.03e-16








#====================================#
#              Figure 1d             #
#====================================#

#Bring in the ITLN1 network genes
itln1Corrs<-read.xlsx("Data/Brown_module_genes.xlsx")

#Bring in and format Krasnow genes and isolate the cell types
Krasnow_all<-readExcelTabs("Data/Krasnow_lung_cell_markers.xlsx")                                                                                                      
#Cull for only immune cells
names(Krasnow_all)[68]<-"Neutrophil"
Krasnow<-Krasnow_all[-grep("SS$",names(Krasnow_all))]
Krasnow<-Krasnow[c(1:12)]
gene2module_krasnow<-data.frame()
#Convert to a gene2module table
for(i in 1:length(Krasnow)){
	gene2module_krasnow<-rbind(gene2module_krasnow,
		data.frame("Gene"=Krasnow[[i]]$Gene[which(Krasnow[[i]]$p_val_adj < 0.05)],
		"Module"=names(Krasnow)[i],stringsAsFactors=F))
}
gene2module_krasnow$Module<-as.character(gene2module_krasnow$Module)
gene2module_krasnow<-gene2module_krasnow[-which(gene2module_krasnow$Module == "Serous" | 
	gene2module_krasnow$Module == "Mucous" ),]

#Do enrichment
enrichment<-DoEnrichment_alt(genes=itln1Corrs$gene,gene2module=gene2module_krasnow,
	background=nrow(secRna_353_norm))

#Plot
pdf("Fig1d.pdf")
par(mar=c(10,5,2,2))
barplot(-log10(enrichment$padj),ylim=c(0,20),las=1,ylab="Enrichment (-log10[FDR])")
abline(h=-log10(0.05),lty=2)
text(x=seq(from=0.7,to=11.5,length.out=10),y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]),
	labels=enrichment$module,srt=45,adj=1,xpd=T,cex=0.8)
dev.off()










#====================================#
#              Figure 1f             #
#====================================#

######### APICAL

#Bring in raw apical secretome counts
apical<-read.table("Data/secretome_count_matrix.txt")

#Specify design
subject <- sapply(strsplit(colnames(apical), split ="_"),function(x)x[1])
treatment <- factor(sapply(strsplit(colnames(apical), split ="_"),function(x)x[2]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(sapply(strsplit(colnames(apical), split ="_"),function(x)x[4]),levels=c("H","A"))  #Asthmatic versus healthy
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_H","IL13_H","BSA_A","IL13_A"))
design <- data.frame(row.names = colnames(apical), subject = subject, treatment = treatment,
	status = status, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = apical,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
apical_norm<-counts(dds_sizeFactors, normalized=T)

#Make boxplot
pdf("Fig1f.1.pdf",height=3.5,width=2.5)
boxplot(simplify2array(apical_norm["ITLN1_ITLN1",])~design$treatment,names=c("",""),
	las=1,ylab="Normalized ITLN-1",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.8,xlab="",
	outline=F,cex.out=0.5,ylim=c(min(simplify2array(apical_norm["ITLN1_ITLN1",])),14))
beeswarm(apical_norm["ITLN1_ITLN1",]~design$treatment, method="swarm", corral="random",
	ylab="",xlab="",cex=0.6,add=T,pch=16)
text(x=1.3:2.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=c("control","IL-13"),srt=0,adj=1,xpd=T,cex=0.8)
dev.off()

#Get p-value
y <- DGEList(counts=apical, group = design$treatment)
y <- calcNormFactors(y)
design.mat<-model.matrix(~subject + treatment, design)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='treatmentIL13')
res_edgeR<-topTags(qlf, n=nrow(apical))
res_edgeR<-res_edgeR[[1]]
res_edgeR["ITLN1","logFC"] #5.77
res_edgeR["ITLN1","PValue"] #6.01e-21
res_edgeR["ITLN1","FDR"] #5.193e-19






######### MUCUS

#Bring in the apical mucus secretome counts
mucus<-read.table("Data/secretome_count_matrix_newMucus.txt")

#Specify design
subject <- sapply(strsplit(colnames(mucus), split ="_"),function(x)x[2])
treatment <- factor(sapply(strsplit(colnames(mucus), split ="_"),function(x)x[5]),levels=c("BSA","IL13")) #BSA versus IL13
status <- factor(sapply(strsplit(colnames(mucus), split ="_"),function(x)x[4]),levels=c("H","A"))  #Asthmatic versus healthy
treatment_status<-factor(paste(treatment,status,sep="_"),levels=c("BSA_H","IL13_H","BSA_A","IL13_A"))
design <- data.frame(row.names = colnames(mucus), subject = subject, treatment = treatment,
	status = status, treatment_status = treatment_status)

#Normalize data
dds <- DESeqDataSetFromMatrix(
  countData = mucus,
  colData = design,
  design = ~subject + treatment)
dds_sizeFactors<-estimateSizeFactors(dds)
mucus_norm<-counts(dds_sizeFactors, normalized=T)

#Make boxplot
pdf("Fig1f.2.pdf",height=3.5,width=2.5)
boxplot(simplify2array(mucus_norm["ITLN1",])~design$treatment,names=c("",""),
	las=1,ylab="Normalized ITLN-1",cex=0.6,medlwd=2,cex.lab=0.8,cex.axis=0.8,xlab="",
	outline=F,cex.out=0.5,ylim=c(min(simplify2array(mucus_norm["ITLN1",])),10))
beeswarm(mucus_norm["ITLN1",]~design$treatment, method="swarm", corral="random",
	ylab="",xlab="",cex=0.6,add=T,pch=16)
text(x=1.3:2.3,y=par()$usr[3]-0.15*(par()$usr[4]-par()$usr[3]),
	labels=c("control","IL-13"),srt=0,adj=1,xpd=T,cex=0.8)
dev.off()

#Get p-value
y <- DGEList(counts=mucus, group = design$treatment)
y <- calcNormFactors(y)
design.mat<-model.matrix(~subject + treatment, design)
y <- estimateDisp(y, design.mat, robust=TRUE)
fit <- glmQLFit(y, design.mat, robust=TRUE)
qlf <- glmQLFTest(fit, coef='treatmentIL13')
res_edgeR<-topTags(qlf, n=nrow(mucus))
res_edgeR<-res_edgeR[[1]]
res_edgeR["ITLN1","logFC"] #4.25
res_edgeR["ITLN1","PValue"] #0.00022
res_edgeR["ITLN1","FDR"] #0.0167

















