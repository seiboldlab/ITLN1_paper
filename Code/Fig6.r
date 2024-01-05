library(openxlsx)
library(ggpubr)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
library(DESeq2)


#====================================#
#              Figure 6a             #
#====================================#

#Load SARP data
v2_t2.df<-readRDS("Data/v2_t2.df.rds")

### Make box plot
pdf("Fig6a.pdf", height=5, width=5)
g <- ggboxplot(v2_t2.df[!is.na(v2_t2.df$rs4656959_GT),], x="rs4656959_GT", y="log2_ITLN1", 
	facet.by="T2_status", outlier.shape=NA, color="T2_status", palette=c("blue", "red")) + geom_beeswarm(aes(color=T2_status)) 
print(g)
dev.off()


#Do differential expression analysis
#Load SARP raw expression counts
expr_raw <- as.matrix(read.table("Data/SARP_raw_counts.txt", header=T, sep='\t', row.names=1, strings=F))
#QC genes
good_genes <- rowSums(expr_raw >= 6) >= (ncol(expr_raw) * .10)
expr_raw_filt <- expr_raw[good_genes,]
expr_raw_filt <- expr_raw_filt[-grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^RMRP|^MRPS|^RPL|^RPS|^ENSG", 
	rownames(expr_raw_filt)),]

#Model genotype and T2 status and their interaction
dds_by_T2 <- DESeqDataSetFromMatrix(
    countData = expr_raw_filt[, rownames(v2_t2.df[!is.na(v2_t2.df$rs4656959),])],
    colData = v2_t2.df[!is.na(v2_t2.df$rs4656959),],
    design = ~rs4656959_GT + T2_status + T2_status:rs4656959_GT )
dds_by_T2 <- DESeq(dds_by_T2)

#Get LFC in GG vs AA, in T2-Low
resultsNames(dds_by_T2)
GG_effect_acrossT2_res <- results(dds_by_T2, contrast=c(0,0,1,0,0,0))
GG_effect_acrossT2_res["ITLN1",]
#       baseMean log2FoldChange     lfcSE      stat    pvalue      padj
#      <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
#ITLN1   105.819      -0.864926  0.420707  -2.05589 0.0397935  0.870301
#1/(2^ -0.864926)
#[1] 1.821246

#Get LFC in GG vs AA, in T2-Hi
resultsNames(dds_by_T2)
GG_effect_acrossT2_res <- results(dds_by_T2, contrast=c(0,0,1,0,0,1))
GG_effect_acrossT2_res["ITLN1",]
#       baseMean log2FoldChange     lfcSE      stat      pvalue      padj
#      <numeric>      <numeric> <numeric> <numeric>   <numeric> <numeric>
#ITLN1   105.819       -1.66438  0.383754  -4.33709 1.44379e-05  0.179306
#1/(2^-1.66438)
#[1] 3.169774








#====================================#
#              Figure 6b-c           #
#====================================#

#Load SARP mucus plugging data
melted_muc_phen<-readRDS("Data/melted_muc_phen.rds")

#Make box plots
pdf("Fig6b.pdf", height=5, width=5)
g <- ggboxplot(melted_muc_phen, x="T2_status", y="ceil_mean_seg_score", color="T2_status", 
	outlier.shape=NA, palette=c("blue", "red")) + geom_beeswarm(mapping=aes(color=T2_status))
print(g)
dev.off()

pdf("Fig6c.pdf", height=5, width=5)
g <- ggboxplot(melted_muc_phen, x="rs4656959_GT", y="ceil_mean_seg_score", facet.by="T2_status",color="T2_status", 
	outlier.shape=NA, palette=c("blue", "red")) + geom_beeswarm(mapping=aes(color=T2_status))
print(g)
dev.off()











