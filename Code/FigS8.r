library(ggbeeswarm)
library(cowplot)
library(ggpubr)
library(openxlsx)
library(lme4)
library(lmerTest)
library(ggbeeswarm)




#====================================#
#              Figure S8a-b          #
#====================================#

#Load clinical information for SARP
phen <- readRDS("Data/v2_t2.df.rds")

#Plot
pdf("FigS8ab.pdf", height=5, width=5)
g <- ggboxplot(phen[!is.na(phen$rs4656959_GT),], x="rs4656959_GT", y="log2_MUC5AC", facet.by="T2_status", outlier.shape=NA, color="T2_status", palette=c("blue", "red")) + geom_beeswarm(aes(color=T2_status)) + stat_compare_means(comparisons=list(c("A/A", "G/A"), c("A/A", "G/G")),method=("wilcox"))
print(g)
g <- ggboxplot(phen[!is.na(phen$rs4656959_GT),], x="rs4656959_GT", y="log2_MUC5B", facet.by="T2_status", outlier.shape=NA, color="T2_status", palette=c("blue", "red")) + geom_beeswarm(aes(color=T2_status)) + stat_compare_means(comparisons=list(c("A/A", "G/A"), c("A/A", "G/G")),method=("wilcox"))
print(g)
dev.off()





#====================================#
#              Figure S8cd           #
#====================================#

#Load clinical information for SARP
phen <- readRDS("Data/v2_t2.df.rds")

#Pearson correlation
cor.test(x=phen$ITLN1_vst,y=phen$MUC5AC_vst)$estimate #0.239
cor.test(x=phen$ITLN1_vst,y=phen$MUC5AC_vst)$p.value #9.68e-05

cor.test(x=phen$ITLN1_vst,y=phen$MUC5B_vst)$estimate #0.1139
cor.test(x=phen$ITLN1_vst,y=phen$MUC5B_vst)$p.value #0.0666

pdf("FigS8c.pdf",width=4,height=4)
par(bty="l")
plot(phen$MUC5AC_vst,phen$ITLN1_vst,las=1, xlab="MUC5AC VST expression",ylab="ITLN1 VST expression",cex=0.6)
abline(lm(phen$ITLN1_vst~phen$MUC5AC_vst))
dev.off()

pdf("FigS8d.pdf",width=4,height=4)
par(bty="l")
plot(phen$MUC5B_vst,phen$ITLN1_vst,las=1, xlab="MUC5B VST expression",ylab="ITLN1 VST expression",cex=0.6)
abline(lm(phen$ITLN1_vst~phen$MUC5B_vst))
dev.off()











