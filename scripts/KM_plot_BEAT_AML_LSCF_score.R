library(survival)
library(survminer)
library(dplyr)

projectDir <- "C:/Users/informatics3/Documents/GitHub/LSCA"

clinical_BeatAML_FP <- paste0(projectDir, "/survival_data/NIHMS1504008-supplement-BEAT_AML_Clinical_Data.txt")
clinical_BeatAML    <- read.delim2(clinical_BeatAML_FP)
clinical_BeatAML    <- clinical_BeatAML %>% select(LabId, overallSurvival, vitalStatus)


subFP              <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_beat_aml.txt"
ctFractionsFP      <- paste0(projectDir, subFP)
ctFractions        <- read.delim2(ctFractionsFP, check.names=FALSE)
ctFractions[,2:13] <- apply(ctFractions[,2:13], 2, as.numeric)


# LSCA score calculation
CellType    <- c("GMP", "CMP", "RApos", "MEP", "MPP") # 50 DEGs
Coefficient <- c(-2.1480989, -1.6350980, 0.3733158, 0.4942520, 4.5197050)

LSCAscoreCoef             <- cbind(CellType, Coefficient)
LSCAscoreCoef             <- as.data.frame(LSCAscoreCoef)
LSCAscoreCoef$Coefficient <- as.numeric(LSCAscoreCoef$Coefficient)

fractions <- ctFractions[, CellType]
LSCAscore <- c()
for(i in 1:nrow(fractions)){
  temp      <- sum(fractions[i, ] * LSCAscoreCoef$Coefficient)
  LSCAscore <- c(LSCAscore, temp)
}
Group              <- ifelse(LSCAscore > median(LSCAscore), "High", "Low")
LabId              <- ctFractions$Mixture
beat_aml           <- cbind(LabId, LSCAscore, Group)
beat_aml           <- as.data.frame(beat_aml)
beat_aml$LSCAscore <- as.numeric(beat_aml$LSCAscore)


# make survival input
beat_aml <- beat_aml %>% left_join(clinical_BeatAML, by="LabId")

beat_aml$overallSurvival <- as.numeric(beat_aml$overallSurvival)
beat_aml$vitalStatus     <- ifelse(beat_aml$vitalStatus == "Alive", 0, 1)

surv_object <- Surv(time=beat_aml$overallSurvival, event=beat_aml$vitalStatus)
fit1        <- survfit(surv_object ~ Group, data = beat_aml)
pValue      <- surv_pvalue(fit = fit1, data = beat_aml)
pValue      <- format(pValue$pval, scientific=T, digits=2)


plotTitle <- "Beat AML, LSCA"
ggsurv    <- ggsurvplot(title=plotTitle, fit1, data=beat_aml, xlab="Days", 
                        pval=paste0("Log-rank\np = ", pValue), pval.coord=c(3000,0.6), pval.size=10,
                        legend.title="", legend.labs=c("High", "Low"), palette=c("maroon", "forestgreen"),
                        risk.table=TRUE, risk.table.col="strata", conf.int=TRUE, risk.table.fontsize=7,
                        ggtheme=theme(plot.title=element_text(size=25, hjust=0.5, face="bold"),
                                      axis.text=element_text(size=20, color="black"),
                                      legend.text=element_text(size=25, color="black"),
                                      axis.title=element_text(size=25, color="black"),
                                      panel.background=element_rect(fill="white", colour="grey50"),
                                      panel.grid=element_line(colour="grey")))
resultFP <- paste0(projectDir, "/plot/KM_plot/KM_Plot_BEAT_AML_LSCA_score_DEGs50.tiff")
tiff(filename=resultFP, width=17.78, height=17.78, res=300, units="cm")
ggsurv
dev.off()