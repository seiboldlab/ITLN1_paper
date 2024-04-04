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











