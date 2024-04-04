library(ggbeeswarm)
library(cowplot)
library(ggpubr)
library(openxlsx)
library(lme4)
library(lmerTest)
library(ggbeeswarm)



#====================================#
#              Figure S6a-b          #
#====================================#

#Read in GALA phenotype
phen<-readRDS("Data/phen.rds")

#Pearson correlation
cor.test(x=phen$ITLN1_vst,y=phen$MUC5AC_vst)$estimate #0.416
cor.test(x=phen$ITLN1_vst,y=phen$MUC5AC_vst)$p.value #1.81e-30

cor.test(x=phen$ITLN1_vst,y=phen$MUC5B_vst)$estimate #-0.2158
cor.test(x=phen$ITLN1_vst,y=phen$MUC5B_vst)$p.value #9.10e-09

pdf("FigS6a.pdf",width=4,height=4)
par(bty="l")
plot(phen$MUC5AC_vst,phen$ITLN1_vst,las=1, xlab="MUC5AC VST expression",ylab="ITLN1 VST expression",cex=0.6)
abline(lm(phen$ITLN1_vst~phen$MUC5AC_vst))
dev.off()

pdf("FigS6b.pdf",width=4,height=4)
par(bty="l")
plot(phen$MUC5B_vst,phen$ITLN1_vst,las=1, xlab="MUC5B VST expression",ylab="ITLN1 VST expression",cex=0.6)
abline(lm(phen$ITLN1_vst~phen$MUC5B_vst))
dev.off()









#====================================#
#              Figure S6c-e          #
#====================================#

#Bring in GALA dataset
phen<-readRDS("Data/phen.rds")

#Plot MUC5AC and MUC5B
pdf("FigS6cd.pdf")
g <- ggplot(phen[!is.na(phen$rs4656959_GT),], aes(x=rs4656959_GT, y=MUC5AC_log2_norm, color=type2_status)) + geom_boxplot()
g <- g + facet_wrap(~ type2_status) + geom_beeswarm() + theme_cowplot() + scale_color_manual(values=c("blue", "red")) + stat_compare_means(comparisons=list(c("A/A", "G/A"), c("A/A", "G/G")),method=("wilcox"))
print(g)
g <- ggplot(phen[!is.na(phen$rs4656959_GT),], aes(x=rs4656959_GT, y=MUC5B_log2_norm, color=type2_status)) + geom_boxplot()
g <- g + facet_wrap(~ type2_status) + geom_beeswarm() + theme_cowplot() + scale_color_manual(values=c("blue", "red")) + stat_compare_means(comparisons=list(c("A/A", "G/A"), c("A/A", "G/G")),method=("wilcox"))
print(g)
dev.off()

#Plot the T2 mucus secretory network
pdf("FigS6e.pdf")
g <- ggplot(phen[!is.na(phen$rs4656959_GT),], aes_string(x="rs4656959_GT", y="MEpink", color="type2_status")) + geom_boxplot()
g <- g + facet_wrap(~ type2_status) + geom_beeswarm() + theme_cowplot() + scale_color_manual(values=c("blue", "red")) + stat_compare_means(comparisons=list(c("A/A", "G/A"), c("A/A", "G/G")),method=("wilcox"))
print(g)
dev.off()








