library(dplyr)
library(survival)
library(survminer)

projectDir <- "C:/Users/informatics3/Documents/GitHub/LSCA"
survivalFP <- paste0(projectDir, "/survival_and_clinical_data/GSE37642_Survival_data.txt")
survival   <- read.delim2(survivalFP)

# GPL96
subFP         <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_GSE37642_GPL96.txt"
ctFrations_FP <- paste0(projectDir, subFP)
ctFrations    <- read.delim2(ctFrations_FP, check.names=FALSE)
GEO_ID        <- ctFrations$Mixture

LSCA              <- ctFrations[, c("GMP", "CMP", "RApos", "MEP", "MPP")] # 50 DEGs
coefficients_LSCA <- c(-2.1480989, -1.6350980, 0.3733158, 0.4942520, 4.5197050)
LSCA              <- apply(LSCA, 2, as.numeric)
rownames(LSCA)    <- GEO_ID
LSCAScore         <- c()
for(i in 1:nrow(LSCA)){
  tempScore <- sum(LSCA[i, ] * coefficients_LSCA)
  LSCAScore  <- c(LSCAScore, tempScore)
}
LSCAScore           <- cbind(GEO_ID, LSCAScore)
LSCAScore           <- as.data.frame(LSCAScore)
LSCAScore$LSCAScore <- as.numeric(LSCAScore$LSCAScore)

medScore                 <- median(LSCAScore$LSCAScore)
Group_by_Score           <- rep(NA, 422)
Inx_high                 <- which(LSCAScore$LSCAScor > medScore)
Inx_low                  <- which(LSCAScore$LSCAScor <= medScore)
Group_by_Score[Inx_high] <- "High"
Group_by_Score[Inx_low]  <- "Low"
LSCAScore                <- cbind(LSCAScore, Group_by_Score)

# KM plot
gse37642                       <- LSCAScore %>% left_join(survival, by="GEO_ID")
gse37642$overall_survival_days <- as.numeric(gse37642$overall_survival_days)
gse37642$life.status           <- ifelse(gse37642$life.status == "alive", 0, 1)

surv_object <- Surv(time=gse37642$overall_survival_days, event=gse37642$life.status)
fit1        <- survfit(surv_object ~ Group_by_Score, data=gse37642)
pValue      <- surv_pvalue(fit=fit1, data=gse37642)
pValue      <- format(pValue$pval, scientific=T, digits=2)

plotTitle <- "GSE37642 (GPL96), LSCA"
ggsurv    <- ggsurvplot(title=plotTitle, fit1, data=gse37642, xlab="Days", 
                        pval=paste0("Log-rank\np = ", pValue), pval.coord=c(3300,0.6), pval.size=10,
                        legend.title="", legend.labs=c("High", "Low"), palette=c("maroon", "forestgreen"),
                        risk.table=TRUE, risk.table.col="strata", conf.int = TRUE, risk.table.fontsize=7,
                        ggtheme=theme(plot.title=element_text(size=25, hjust=0.5, face="bold"),
                                      axis.text=element_text(size=20, color="black"),
                                      legend.text=element_text(size=25, color="black"),
                                      axis.title=element_text(size=25, color="black"),
                                      panel.background=element_rect(fill="white", colour="grey50"),
                                      panel.grid=element_line(colour="grey")))
resultFP <- paste0(projectDir, "/plot/KM_plot/KM_Plot_GSE37642_LSCA_score_DEGs50.tiff")
tiff(filename=resultFP, width=7, height=7, res=300, units="in")
ggsurv
dev.off()
