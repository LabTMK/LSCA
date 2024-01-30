library(survival)
library(dplyr)
library(forestmodel)

projectDir <- "C:/Users/informatics3/Documents/GitHub/LSCA"


# coefficients
CellType    <- c("GMP", "CMP", "RApos", "MEP", "MPP") # 50 DEGs
Coefficient <- c(-2.1480989, -1.6350980, 0.3733158, 0.4942520, 4.5197050)

coef1             <- cbind(CellType, Coefficient)
coef1             <- as.data.frame(coef1)
coef1$Coefficient <- as.numeric(coef1$Coefficient)

# cell fractions
subFP             <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_tcga_laml_RPKM.txt"

ctFrationsFP <- paste0(projectDir, subFP)
ctFrations   <- read.delim2(ctFrationsFP, check.names=FALSE)

newScore1 <- c()
for(i in 1:nrow(ctFrations)){
  temp_ctFrations <- as.numeric(ctFrations[i, coef1$CellType])
  tempScore       <- sum(temp_ctFrations * coef1$Coefficient)
  newScore1       <- c(newScore1, tempScore)
}
bcr_patient_barcode <- ctFrations$Mixture
newScore1           <- cbind(bcr_patient_barcode, newScore1)
newScore1           <- as.data.frame(newScore1)


# clinical data
subFP                <- "/survival_and_clinical_data/TCGA_clinical_patient_laml.tsv"
clinicalData_TCGA_FP <- paste0(projectDir, subFP)
clinicalData_TCGA    <- read.delim2(clinicalData_TCGA_FP)

clinicalData_TCGA <- clinicalData_TCGA %>% select(bcr_patient_barcode,
                                                  age_at_initial_pathologic_diagnosis,
                                                  gender,
                                                  acute_myeloid_leukemia_calgb_cytogenetics_risk_category)

# clinical data 2 (OS.time)
tcgaSurvivalFP <- "/survival_and_clinical_data/TCGA_CDR.txt"
tcgaSurvival   <- read.delim2(paste0(projectDir, tcgaSurvivalFP))
tcgaSurvival   <- tcgaSurvival %>% select(bcr_patient_barcode, OS.time, vital_status)

tcgaSurvival$OS.time      <- as.numeric(tcgaSurvival$OS.time)
tcgaSurvival$vital_status <- ifelse(tcgaSurvival$vital_status == "Alive", 0, 1)

clinicalData_TCGA <- clinicalData_TCGA %>% left_join(tcgaSurvival, by="bcr_patient_barcode")
TFindex           <- clinicalData_TCGA$bcr_patient_barcode %in% newScore1$bcr_patient_barcode
clinicalData_TCGA <- clinicalData_TCGA[TFindex, ]
clinicalData_TCGA <- clinicalData_TCGA %>% left_join(newScore1, by="bcr_patient_barcode")


# trimming clinicalData_TCGA
colnames(clinicalData_TCGA) <- c("bcr_patient_barcode", "Age", "Sex", "Cyto_risk", "OS.time", "Vital_Status", "LSCA_score")
clinicalDataForPlot         <- clinicalData_TCGA[, -1]
clinicalDataForPlot$Sex  <- gsub("^MALE", "Male", clinicalDataForPlot$Sex)
clinicalDataForPlot$Sex  <- gsub("FEMALE", "Female", clinicalDataForPlot$Sex)
clinicalDataForPlot$Sex  <- factor(clinicalDataForPlot$Sex, levels=c("Male", "Female"))
cytogeneticRiskNAInx        <- which(clinicalDataForPlot$Cyto_risk == "[Not Available]")

clinicalDataForPlot$Cyto_risk[cytogeneticRiskNAInx] <- NA
clinicalDataForPlot$Cyto_risk                       <- gsub("Intermediate/Normal", "Intermediate", clinicalDataForPlot$Cyto_risk)
clinicalDataForPlot$Cyto_risk                       <- as.factor(clinicalDataForPlot$Cyto_risk)

clinicalDataForPlot$LSCA_score <- as.numeric(clinicalDataForPlot$LSCA_score)

TFindex                           <- clinicalDataForPlot$Age >= 60
clinicalDataForPlot$Age[TFindex]  <- "กร 60"
clinicalDataForPlot$Age[!TFindex] <- "< 60"
clinicalDataForPlot$Age           <- as.factor(clinicalDataForPlot$Age)


# Hazard Ratio
model    <- coxph(Surv(OS.time, Vital_Status) ~ ., clinicalDataForPlot)
resultFP <- paste0(projectDir, "/plot/forest_plot/TCGA_LAML_LSCA_score_DEGs50.tiff")
tiff(filename=resultFP, width=12, height=9.6, res=300, units="in")
forest_model(model, limits=c(log(0.1), log(8)),
             format_options=forest_model_format_options(text_size=9)) + #7
theme(axis.text.x=element_text(size=15, colour="black")) #15
dev.off()