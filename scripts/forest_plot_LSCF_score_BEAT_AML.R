library(survival)
library(dplyr)
library(forestmodel)

projectDir <- "C:/Users/informatics3/Documents/GitHub/LSCA"


# Beat AML 
clinical_BeatAML_FP <- paste0(projectDir, "/survival_and_clinical_data/NIHMS1504008-supplement-BEAT_AML_Clinical_Data.txt")
clinical_BeatAML    <- read.delim2(clinical_BeatAML_FP)
clinical_BeatAML    <- clinical_BeatAML %>% select(LabId, ageAtDiagnosis, consensus_sex, ELN2017, overallSurvival, vitalStatus)


# LSCA score
# Coefficients
CellType    <- c("GMP", "CMP", "RApos", "MEP", "MPP") # 50 DEGs
Coefficient <- c(-2.1480989, -1.6350980, 0.3733158, 0.4942520, 4.5197050)

coef             <- cbind(CellType, Coefficient)
coef             <- as.data.frame(coef)
coef$Coefficient <- as.numeric(coef$Coefficient)

# cell type fractions
subFP              <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_beat_aml.txt"
ctFrations_BEAT_FP <- paste0(projectDir, subFP)
ctFrations_BEAT    <- read.delim2(ctFrations_BEAT_FP, check.names=FALSE)

LSCA <- c()
for(i in 1:nrow(ctFrations_BEAT)){
  temp_ctfractions <- as.numeric(ctFrations_BEAT[i, coef$CellType])
  temp_coef        <- as.numeric(coef$Coefficient)
  temp_score       <- sum(temp_ctfractions * temp_coef)
  LSCA             <- c(LSCA, temp_score)
}
LabId <- ctFrations_BEAT$Mixture
LSCA  <- cbind(LabId, LSCA)
LSCA  <- as.data.frame(LSCA)


# Merge
clinicalDataForPlot <- clinical_BeatAML %>% left_join(LSCA, by="LabId")
TFindex             <- !(is.na(clinicalDataForPlot$LSCA))
clinicalDataForPlot <- clinicalDataForPlot[TFindex, ]

clinicalDataForPlot$ageAtDiagnosis <- ifelse(clinicalDataForPlot$ageAtDiagnosis >= 60, "กร 60", "< 60")
clinicalDataForPlot$ageAtDiagnosis <- factor(clinicalDataForPlot$ageAtDiagnosis, levels=c("< 60", "กร 60"))

clinicalDataForPlot$consensus_sex <- factor(clinicalDataForPlot$consensus_sex, levels=c("Male", "Female"))

NAindex                               <- which(clinicalDataForPlot$ELN2017 == "Unknown")
clinicalDataForPlot$ELN2017[NAindex]  <- NA
favorInx                              <- which(clinicalDataForPlot$ELN2017 == "FavorableOrIntermediate")
clinicalDataForPlot$ELN2017[favorInx] <- "Favorable"
interInx                              <- which(clinicalDataForPlot$ELN2017 == "IntermediateOrAdverse")
clinicalDataForPlot$ELN2017[interInx] <- "Intermediate"
interInx                              <- which(clinicalDataForPlot$ELN2017 == "Adverse")
clinicalDataForPlot$ELN2017[interInx] <- "Poor"
clinicalDataForPlot$ELN2017           <- factor(clinicalDataForPlot$ELN2017, levels=c("Favorable", "Intermediate", "Poor"))

clinicalDataForPlot$vitalStatus <- ifelse(clinicalDataForPlot$vitalStatus == "Alive", 0, 1)

clinicalDataForPlot$LSCA <- as.numeric(clinicalDataForPlot$LSCA)

clinicalDataForPlot <- clinicalDataForPlot[, -1]

colnames(clinicalDataForPlot) <- c("Age", "Sex", "Cyto_risk", "overallSurvival", "vitalStatus", "LSCA_score")


# Hazard Ratio
model    <- coxph(Surv(overallSurvival, vitalStatus) ~ ., clinicalDataForPlot)
subFP    <- "/plot/forest_plot/BEAT_AML_LSCA_score_DEGs50.tiff"
resultFP <- paste0(projectDir, subFP)
tiff(filename=resultFP, width=12, height=9.6, res=300, units="in")
forest_model(model, limits = c(log(0.1), log(7)), 
             format_options=forest_model_format_options(text_size=9)) + #7
theme(axis.text.x=element_text(size=15, colour="black")) #15
dev.off()