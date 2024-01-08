library(tidyverse)
library(WGCNA)
library(DESeq2)
library(edgeR)
library(beeswarm)
library(openxlsx)
library(enrichR)

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
vstMat <- assay(varianceStabilizingTransformation(dds_Rna))

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
#              Figure 1c             #
#====================================#

############################# Do WGCNA

allData<-t(vstMat)
#Pick power
powers <- 1:20
sft = pickSoftThreshold(allData, powerVector = powers, verbose = 5)
#Plot
pdf(file = "WGCNA_softThreshold.bF.pdf", width = 9, height = 5);
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
        main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
softPower <- 10

#Create Similarity Matrix
pearson <- WGCNA::cor(as.matrix(allData),method="pearson")

#Convert Similarity Matrix to Adjacency Matrix using Adjacency Function
adjacency.p <- adjacency.fromSimilarity(pearson,type = "signed",power=softPower)

#Convert Adjacency to TOM dissimilarity
TOMdissim.p <- 1 - TOMsimilarity(adjacency.p,TOMType = "signed",TOMDenom = "min")

#Perform hierarchical clustering on the dissimilarity matrix
geneTree <- hclust(as.dist(TOMdissim.p), method = "average")

#Cut tree
x <- 0.5
modules <- cutreeDynamic(dendro = geneTree,method='hybrid', distM = TOMdissim.p,pamStage=F, 
	pamRespectsDendro = T, deepSplit=2, cutHeight=0.95, minClusterSize=40)
modcolors <- labels2colors(modules)

#Merge similar modules
modMerged <- mergeCloseModules(allData, modcolors, cutHeight=0.1)
moduleColors <- modMerged$colors 

#Make gene2module table
gene2module <- data.frame(gene=colnames(allData), module=moduleColors)
gene2module <- na.omit(gene2module)

#Calculate module eigengenes 
datME<-moduleEigengenes(allData,moduleColors)$eigengenes

#Get eigengene-based connectivity matrix, based on the eigengenes
datKME <- signedKME(allData, datME)

#Get hub genes
ADJ1 <- abs(cor(allData, use="p"))^softPower
Alldegrees1<-intramodularConnectivity(ADJ1, gene2module$module)
hub_genes<-data.frame()
for(i in 1:length(unique(gene2module$module))){
    Alldegrees1$Module = gene2module$module
    tmp = Alldegrees1[Alldegrees1$Module == unique(gene2module$module)[i], ]
    hub_genes<-rbind(hub_genes, head(tmp[order(tmp$kWithin, decreasing=T),], n=nrow(tmp)))
}

#Create master module table
gene2module<-read.table("Data/WGCNA.gene2module.txt",header=T,strings=F)
gene2module_with_cor <- gene2module
gene2module_with_cor$module<-factor(gene2module_with_cor$module,levels=c(sort(unique(gene2module_with_cor$module))))
gene2module_with_cor$cor <- NA
for(i in unique(gene2module_with_cor$module)) {
    kME_name <- paste0("kME",i)
    idx <- which(gene2module_with_cor$module==i)
    gene.idx <- as.character(gene2module_with_cor[idx,"gene"])
    gene2module_with_cor$cor[idx] <- datKME[gene.idx,kME_name]
}
#Add in Kwithin, etc
gene2module_with_cor<-cbind(gene2module_with_cor,hub_genes[as.character(gene2module_with_cor$gene),])
gene2module_with_cor<-gene2module_with_cor[,-ncol(gene2module_with_cor)]
#Sort by module, then by correlation
gene2module_with_cor<-gene2module_with_cor[with(gene2module_with_cor,order(module,-abs(cor))),]
#Place into lists
gene2module_with_cor_list<-list()
for(i in 1:length(levels(gene2module_with_cor$module))){
	gene2module_with_cor_list[[i]]<-gene2module_with_cor[which(gene2module_with_cor$module == 
	levels(gene2module_with_cor$module)[i]),]
	names(gene2module_with_cor_list)[i]<-levels(gene2module_with_cor$module)[i]
}




############################# Plot ITLN1 network with focus on ITLN1 connections

#Create master module table for the brown module (Table S2)

#Subselect brown module
genes<-gene2module$gene[which(gene2module$module == "brown")]
genes<-genes[-grep("ITLN1",genes)]
#Get spearman correlation of all genes with ITLN1
itln1Corrs<-apply(secRna_353_norm[genes,],1,
	function(x)cor.test(secRna_353_norm["ITLN1",],x,method="spearman"))
#Pull out the lists of correlations and pvalues and put into a dataframe
tempMat<-as.data.frame(matrix(ncol=2,nrow=0))
itln1Corrs<-t(sapply(itln1Corrs,function(x)rbind(tempMat,as.data.frame(matrix(c(x$estimate,x$p.value),nrow=1)))))
colnames(itln1Corrs)<-c("correlation","pvalue")
itln1Corrs<-as.data.frame(itln1Corrs)
itln1Corrs[,1]<-simplify2array(itln1Corrs[,1])
itln1Corrs[,2]<-simplify2array(itln1Corrs[,2])
#Order and rank
itln1Corrs<-itln1Corrs[order(itln1Corrs$correlation,decreasing=T),]
itln1Corrs<-data.frame("gene"=rownames(itln1Corrs),"ITLN1_cor"=itln1Corrs$correlation,
	"pvalue_cor"=itln1Corrs$pvalue,"qvalue_cor"=p.adjust(itln1Corrs$pvalue,method="fdr"),
	"rank"=seq(nrow(itln1Corrs)),row.names=seq(nrow(itln1Corrs)),stringsAsFactors=F)
#Add in ITLN1
itln1Corrs<-rbind(data.frame("gene"="ITLN1","ITLN1_cor"=1,"pvalue_cor"=0,"qvalue_cor"=0,"rank"=0),itln1Corrs)
#Add in cor (KME) and kWithin values too
tempTab<-gene2module_with_cor_list$brown
rownames(tempTab)<-tempTab$gene
itln1Corrs<-cbind(itln1Corrs,KME=tempTab[itln1Corrs$gene,]$cor,kWithin=tempTab[itln1Corrs$gene,]$kWithin)
#Also, read in differential expression results
diffexp<-res_Rna353
itln1Corrs<-cbind(itln1Corrs,"logFC_IL13"=diffexp[itln1Corrs$gene,]$logFC,"FDR_IL13"=diffexp[itln1Corrs$gene,]$FDR)
#Save
write.xlsx(itln1Corrs,file="Brown_module_genes.xlsx")



#Export network nodes and edges for cytoscape

#Subselect the genes that are significantly correlated with 
#ITLN1 with cor > 0.5, but also KME > 0.85
currModule<-"brown"
itln1Corrs_culled<-itln1Corrs[which(itln1Corrs$qvalue_cor < 0.05 & itln1Corrs$ITLN1_cor > 0.5 & itln1Corrs$KME > 0.85),]
itln1Corrs_culled<-rbind(itln1Corrs[1,],itln1Corrs_culled)
rownames(itln1Corrs_culled)<-itln1Corrs_culled$gene
modProbes<-itln1Corrs_culled$gene
#Select the corresponding Topological Overlap
modTOM = TOMdissim.p[which(colnames(t(vstMat)) %in% modProbes),which(colnames(t(vstMat)) %in% modProbes)]
dimnames(modTOM) = list(modProbes, modProbes)
#Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
	edgeFile = "CytoscapeInput-edges-brown.txt",
	nodeFile = "CytoscapeInput-nodes-brown.txt",
	weighted = TRUE,
	threshold = 0,
	nodeNames = modProbes,
	nodeAttr = itln1Corrs_culled)


#Do pathway enrichment on all the genes
dapsa_modules_forEnrichr<-gene2module[-which(gene2module$module == "grey"),]
dapsa_modules_forEnrichr<-dapsa_modules_forEnrichr[order(dapsa_modules_forEnrichr$module),]
colnames(dapsa_modules_forEnrichr)<-c("gene","comparison")
dataset<-"dapsa_modules"
doEnrichOneAtATime(DEG_table=dapsa_modules_forEnrichr,dataset=dataset,Enrichr_dir="./")

#Subset only genes related to select functional pathways and gene sets
ench<-read.xlsx("Data/EnrichmentsToPlot.xlsx")
enchTab<-data.frame()
for(i in 1:nrow(ench)){
	enchTab<-rbind(enchTab,data.frame("Genes"=strsplit(ench[i,]$Genes,", ")[[1]],
		"Term"=ench[i,]$Term,stringsAsFactors=F))
}
#Remove duplicate rows
enchTab<-enchTab[-which(duplicated(enchTab$Gene)),]
#Also add in top IL-13 pathology genes (removed due to redundancy: "PRB1","PRB2","CST1","CLCA1","FETUB","LYZ")
IL13genes<-c("SH2D1B","CDH26","ITLN1","CISH","IL1RN","SPDEF","NOS2","DPP4","FOXA3","ALOX15","POSTN","FCGBP","CAPN14")
enchTab<-rbind(enchTab,data.frame("Genes"=IL13genes,"Term"="IL13",stringsAsFactors=F))
#Bring in the master module table, subset to enchTab, and add in the Terms
modTab<-itln1Corrs
rownames(modTab)<-modTab$gene
modTab_culled<-modTab[enchTab$Genes,]
modTab_culled<-cbind(modTab_culled,"Term"=enchTab$Term)


#Export subsetted genes to cytoscape
#Input
modProbes<-modTab_culled$gene
currModule<-"brown"
#Select the corresponding Topological Overlap
modTOM = TOMdissim.p[which(colnames(t(vstMat)) %in% modProbes),which(colnames(t(vstMat)) %in% modProbes)]
dimnames(modTOM) = list(modProbes, modProbes)
#Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
	edgeFile = "CytoscapeInput-edges-FUNCTIONALGENES-brown.txt",
	nodeFile = "CytoscapeInput-nodes-FUNCTIONALGENES-brown.txt",
	weighted = TRUE,
	threshold = 0,
	nodeNames = modProbes,
	nodeAttr = modTab_culled)


#Finally, cull edges for background (non-ILTL1 connected) genes
#Include only the most-highly weighted edges. 
#More strongly connected nodes get to keep a higher proportion of edges.
#filterDegree gives the degree to which edges should be filtered.
cullConnections<-function(br,filterDegree=NULL){
	#Add in vector of normalized weights
	br<-cbind(br,"weight_norm"=(br$weight - min(br$weight)) / (max(br$weight) - min(br$weight)))
	#Now, for each node
	br_culled<-data.frame()
	for(i in 1:length(unique(c(br$fromNode,br$toNode)))){
		#Make a temporary table with to and from rows involving the current gene
		currTab<-br[which(br$fromNode == unique(c(br$fromNode,br$toNode))[i] |
			br$toNode == unique(c(br$fromNode,br$toNode))[i]),]
		if(unique(c(br$fromNode,br$toNode))[i] == "ITLN1"){
			br_culled<-rbind(br_culled,cbind(currTab,"color"="red"))
		}else{
			#Get average normalized weight
			meanWeight<-mean(currTab$weight_norm)
			#Now only keep the rows for which the weight is above the meanWeight quantile
			currTab_culled<-currTab[order(currTab$weight_norm,decreasing=T),]
			currTab_culled<-currTab_culled[1:ceiling((nrow(currTab_culled) * (meanWeight / filterDegree))),]
			br_culled<-rbind(br_culled,cbind(currTab_culled,"color"="grey"))
		}
	}
	return(unique(br_culled))
}
#Bring in edges
br = read.table("CytoscapeInput-edges-FUNCTIONALGENES-brown.txt",h=T)
filterDegree<-10
br_culled<-cullConnections(br=br,filterDegree=filterDegree)
write.table(br_culled, file="CytoscapeInput-edges-FUNCTIONALGENES-CULLED-brown.txt",sep="\t",
	quote=F, row.names=F)












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

















