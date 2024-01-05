library(tidyverse)
library(plotrix)
library(beeswarm)
library(openxlsx)
library(lmerTest)
library(ggbeeswarm)
library(cowplot)
library(DESeq2)
library(Seurat)
library(ggplot2)
library(Matrix)
library(future)



#====================================#
#              Figure 4b             #
#====================================#

#Read in the GALA dataset
phen<-readRDS("Data/phen.rds")

#Box plot
pdf("Fig4b.pdf", width=4, height=4)
g <- ggplot(phen, aes(x=type2_status, y=ITLN1_log2_norm, color=type2_status)) + geom_boxplot() + geom_beeswarm()
g <- g + scale_color_manual(values=c("blue", "red")) + theme_cowplot()
print(g)
dev.off()

#Get p-value and fold change
#Load raw GALA expression counts
rawMat <- read.table("Data/GALA_raw_counts.txt", header=T, strings=F)
#Cull genes
good_genes <- rowSums(rawMat >= 6) >= (ncol(rawMat) * .10)
rawMat_filt <- rawMat[good_genes,]
rawMat_filt <- rawMat_filt[-grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^RMRP|^MRPS|^RPL|^RPS|^ENSG", 
	rownames(rawMat_filt)),]
#DE analysis between t2 status
dds <- DESeqDataSetFromMatrix(
    countData = rawMat[rownames(rawMat_filt), rownames(phen)],
    colData = phen,
    design = ~ type2_status)
dds <- DESeq(dds)
resultsNames(dds)
t2_res <- results(dds, name="type2_status_type2_high_vs_type2_low")
t2_res <- lfcShrink(dds, type="ashr", coef="type2_status_type2_high_vs_type2_low")
t2_res <- data.frame(t2_res[order(t2_res$pvalue),])
t2_res["ITLN1",] #L2FC=5.73, FC=53.1, FDR < 4.4e-285







#====================================#
#              Figure 4c             #
#====================================#

#Read in bronchial brushing cells
dat<-readRDS("Data/BBsc_nonsquamous.rds")

#Sample name
sample_name <- "BB_fresh_wo_squamous"

#Input into Seurat
sDat <- list()
for(sample_id in names(dat)) {
    sDat[[sample_id]] <- CreateSeuratObject(counts = dat[[sample_id]], project = sample_name)
} 

#Add metadata
for(sample_id in names(dat)) {
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = sDat[[sample_id]]), value = TRUE)
    percent.mito <- Matrix::colSums(x=GetAssayData(sDat[[sample_id]], 
    	slot='counts')[mito.genes,])/Matrix::colSums(x=GetAssayData(sDat[[sample_id]], slot='counts'))

    ribo.genes <- grep(pattern = "^RPS|^RPL", x = rownames(x = sDat[[sample_id]]), value = TRUE)
    percent.ribo <- Matrix::colSums(x=GetAssayData(sDat[[sample_id]], 
    	slot='counts')[ribo.genes,])/Matrix::colSums(x=GetAssayData(sDat[[sample_id]], slot='counts'))

    sDat[[sample_id]][['percent.mito']] <- percent.mito
    sDat[[sample_id]][['percent.ribo']] <- percent.ribo
}

#Filter out cells and scTransform
for(sample_id in names(dat)) {
    sDat[[sample_id]] <- subset(x = sDat[[sample_id]], 
    	subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mito < 0.3)
    sDat[[sample_id]] <- SCTransform(sDat[[sample_id]], 
    	vars.to.regress = "percent.mito", verbose = FALSE)
}

#Set multicore
plan("multicore", workers=47)
options(future.globals.maxSize = 6000 * 1024^2)
RhpcBLASctl::omp_set_num_threads(1L)
RhpcBLASctl::blas_set_num_threads(1L)
Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1")

#Perform data integration
var.features <- SelectIntegrationFeatures(object.list = sDat, nfeatures = 5000)
sDat <- PrepSCTIntegration(object.list = sDat, anchor.features = var.features, verbose = FALSE)
anchors <- FindIntegrationAnchors(sDat, dims=1:30, normalization.method = "SCT", anchor.features=var.features)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
DefaultAssay(object = integrated) <- "integrated"
integrated@meta.data$sampleID <- gsub(".*_","",rownames(integrated@meta.data))
integrated@meta.data$sampleID <- plyr::mapvalues(integrated@meta.data$sampleID, from=c(1:2), to=names(dat))

#Prep for clustering
integrated <- RunPCA(object = integrated, npcs = 50, verbose = FALSE)
integrated <- NormalizeData(integrated, assay="RNA")
integrated <- ScaleData(integrated, assay="RNA")

#Cluster cells
nPCAs <- 30
integrated <- FindNeighbors(object = integrated, dims = 1:nPCAs, k.param=20)
integrated <- RunUMAP(object = integrated, reduction = "pca", dims = 1:nPCAs, n.neighbors=20)
integrated3 <- FindClusters(integrated, resolution=0.4, algorithm=4)

#Set cluster colors
cell_type_cols <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"))

#Set cell type names
cell_type_names <- c("MUC5B Mucus Sec.", "Early Sec.", "Ciliated", "MUC5AC Mucus Sec.", "Macrophage", "Diff. Basal",
                     "Monocytes DCs", "Low Quality", "Squamous", "Basal", "Early Ciliating", "B Cells", 
                     "Prolif. Basal", "T Cells", "Glandular Cells", "Mast Cells", "Ionocytes")
cell_type_levels <- c("Basal", "Diff. Basal", "Prolif. Basal", "Early Sec.", "MUC5B Mucus Sec.", 
                    "MUC5AC Mucus Sec.", "Early Ciliating", "Ciliated", "Squamous", "Glandular Cells", 
                    "Ionocytes", "B Cells", "Macrophage", "Mast Cells", "Monocytes DCs", "T Cells", "Low Quality")
integrated3@meta.data$cell_type <- factor(plyr::mapvalues(integrated3@meta.data$seurat_clusters, from=c(1:17),
                                    to=cell_type_names), levels=cell_type_levels)

#Plot clusters
tiff(paste0("Fig4c.tiff"), res=300, width=6*300, height=4*300, compression="lzw")
DimPlot(integrated3, reduction='umap', cols=cell_type_cols, group.by="cell_type", label=F, repel=T) 
dev.off()

#Find cluster markers
DefaultAssay(object = integrated3) <- "RNA"
markers.list <- list()
nClusters <- 17
for(i in 1:nClusters) {
    markers.list[[i]] <- FindMarkers(integrated3, ident.1=i, min.pct=0.1, logfc.threshold=0.25,  only.pos=T, assay="RNA")
    markers.list[[i]] <- markers.list[[i]][markers.list[[i]]$p_val_adj < 0.05,]
}

#Print markers
marker_list_rm_MT_RB <- lapply(1:17, function(i) {
        tmp <- markers.list[[i]]
        tmp$gene <- rownames(tmp)
        tmp <- tmp[!grepl("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^RMRP|^MRPS|^RPL|^RPS|^ENSG", rownames(tmp)),]
        tmp$cell_type <- cell_type_names[i]
        tmp
    })
marker_list_rm_MT_RB.df <- do.call(rbind, marker_list_rm_MT_RB)
wb <- createWorkbook()
addWorksheet(wb, sheetName="Markers")
writeData(wb, x=marker_list_rm_MT_RB.df, sheet="Markers", rowNames=T)
saveWorkbook(wb, file="TableS3.xlsx", overwrite=T)










#====================================#
#              Figure 4d             #
#====================================#

#Plot ITLN1 across cell types
tiff("Fig4d.tiff", res=300, width=6*300, height=4*300, compression="lzw")
VlnPlot(integrated3, c("ITLN1"), assay="RNA", cols=cell_type_cols, group.by="cell_type") + theme(legend.position="none")
dev.off()










#====================================#
#              Figure 4e and f       #
#====================================#

#Read in the GALA dataset
phen<-readRDS("Data/phen.rds")

# eQTL plot
pdf("Fig4ef.pdf")
g <- ggplot(phen[!is.na(phen$rs4656959_GT),], aes(x=rs4656959_GT, y=ITLN1_log2_norm)) + geom_boxplot() 
g <- g + geom_beeswarm() + theme_cowplot()
print(g)
g <- ggplot(phen[!is.na(phen$rs4656959_GT),], aes(x=rs4656959_GT, y=ITLN1_log2_norm, color=type2_status)) + geom_boxplot() 
g <- g + facet_wrap(~ type2_status) + geom_beeswarm() + theme_cowplot() + scale_color_manual(values=c("blue", "red"))
print(g)
dev.off()

#Get p-values and fold changes
#Load GALA raw expression counts
rawMat <- read.table("Data/GALA_raw_counts.txt", header=T, strings=F)
#Cull genes
good_genes <- rowSums(rawMat >= 6) >= (ncol(rawMat) * .10)
rawMat_filt <- rawMat[good_genes,]
rawMat_filt <- rawMat_filt[-grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^RMRP|^MRPS|^RPL|^RPS|^ENSG", 
	rownames(rawMat_filt)),]
# DE analysis between t2 status
dds_GT <- DESeqDataSetFromMatrix(
    countData = rawMat[rownames(rawMat_filt), rownames(phen[!is.na(phen$rs4656959_GT),])],
    colData = phen[!is.na(phen$rs4656959_GT),],
    design = ~ rs4656959_GT + type2_status + rs4656959_GT:type2_status)
dds_GT <- DESeq(dds_GT, parallel=T)
resultsNames(dds_GT)
#[1] "Intercept"
#[2] "rs4656959_GT_G.A_vs_A.A"
#[3] "rs4656959_GT_G.G_vs_A.A"
#[4] "type2_status_type2_high_vs_type2_low"
#[5] "rs4656959_GTG.A.type2_statustype2_high"
#[6] "rs4656959_GTG.G.type2_statustype2_high"
#Run contrasts
### GG vs AA (L2FC=-2.74; FC=6.68; p: 2.84e-34)
GG_AA_res <- results(dds_GT, contrast=c(0,0,1,0,0,0.5))
GG_AA_res <- lfcShrink(dds_GT, type="ashr", contrast=c(0,0,1,0,0,0.5))
GG_AA_res["ITLN1",]

### AG vs AA (L2FC=-0.86; FC=1.81; p: 1.1e-9)
AG_AA_res <- results(dds_GT, contrast=c(0,1,0,0,0.5,0))
AG_AA_res <- lfcShrink(dds_GT, type="ashr", contrast=c(0,1,0,0,0.5,0))
AG_AA_res["ITLN1",]

### GG vs AG (L2FC=-1.85; FC=3.6; p: 5.7e-8)
GG_AG_res <- results(dds_GT, contrast=c(0,-1,1,0,0.5,0.5))
GG_AG_res <- lfcShrink(dds_GT, type="ashr", contrast=c(0,-1,1,0,0.5,0.5))
GG_AG_res["ITLN1",]

### GG vs AA T2 low (L2FC=-1.94, FC: 3.84 , p:8.9e-10)
GG_AA_T2low_res <- results(dds_GT, contrast=c(0,0,1,0,0,0))
GG_AA_T2low_res <- lfcShrink(dds_GT, type="ashr", contrast=c(0,0,1,0,0,0))
GG_AA_T2low_res["ITLN1",]

### GG vs AA T2 high (L2FC=-3.44, FC:10.85, p: 5.30e-30)
GG_AA_T2high_res <- results(dds_GT, contrast=c(0,0,1,0,0,1))
GG_AA_T2high_res <- lfcShrink(dds_GT, type="ashr", contrast=c(0,0,1,0,0,1))
GG_AA_T2high_res["ITLN1",]













