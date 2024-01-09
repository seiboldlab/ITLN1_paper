library(lme4)
library(lmerTest)
library(rms)


#====================================#
#              Table S1              #
#====================================#

#Read in GALA phenotype
phen<-readRDS("Data/phen.rds")

#Multivariate model with genotype
ITLN1_geno_form <- as.formula(paste0("ITLN1_eQTL_norm ~ rs4656959_num + type2_status + rs4656959_num:type2_status + sex + asthma + age"))
ITLN1_multi_cov_geno_fit <- ols(ITLN1_geno_form, data=phen)
ITLN1_multi_cov_geno_fit2 <- lm(ITLN1_geno_form, data=phen)

var_of_interest <- rownames(summary(ITLN1_multi_cov_geno_fit2)$coefficients)[-1]
var_of_interest[c(2,3,4,6)] <- c("type2_status", "sex", "asthma", "rs4656959_num * type2_status")

plt <- plot(anova(ITLN1_multi_cov_geno_fit), what='partial R2', pl=F)
ITLN1_multi_cov_geno.df <- data.frame(predictor=var_of_interest, partial_R2=plt[var_of_interest]*100, summary(ITLN1_multi_cov_geno_fit2)$coef[-1,])
colnames(ITLN1_multi_cov_geno.df) <- c("predictor", "partial_R2", "Estimate", "SE", "t", "p-value")







