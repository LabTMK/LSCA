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
subFP <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_GSE37642_GPL96.txt"

ctFrationsFP <- paste0(projectDir, subFP)
ctFrations   <- read.delim2(ctFrationsFP, check.names=FALSE)

LSCA_score <- c()
for(i in 1:nrow(ctFrations)){
  temp_ctFrations <- as.numeric(ctFrations[i, coef1$CellType])
  tempScore       <- sum(temp_ctFrations * coef1$Coefficient)
  LSCA_score      <- c(LSCA_score, tempScore)
}
GEO_Accession         <- ctFrations$Mixture
LSCA_score            <- cbind(GEO_Accession, LSCA_score)
LSCA_score            <- as.data.frame(LSCA_score)
LSCA_score$LSCA_score <- as.numeric(LSCA_score$LSCA_score)


# clinical data
clinicalDataFP <- paste0(projectDir, "/survival_data/GSE37642_GPL96_phenotype_data.tsv")
clinicalData   <- read.delim2(clinicalDataFP, check.names=FALSE)
clinicalData   <- clinicalData %>% select(geo_accession, `fab:ch1`, `age:ch1`, `runx1-runx1t1_fusion:ch1`, 
                                          `runx1_mutation:ch1`, `overall survival (days):ch1`, `life status:ch1`)

# RUNX1 Fusion = RUNX1 Fusion
colnames(clinicalData) <- c("GEO_Accession", "FAB", "Age", "RUNX1 Fusion", "RUNX1 Mutation", "OS.time", "Vital_Status")
clinicalData$FAB       <- ifelse(is.na(clinicalData$FAB), clinicalData$FAB, paste0("M", clinicalData$FAB))

clinicalDataForPlot <- clinicalData %>% left_join(LSCA_score, by="GEO_Accession")
clinicalDataForPlot <- clinicalDataForPlot[, -1]

clinicalDataForPlot$FAB <- factor(clinicalDataForPlot$FAB, levels=c("M3", "M0", "M1", "M2", "M3v", "M4", "M5", "M6", "M7"))

TFindex                           <- clinicalDataForPlot$Age >= 60
clinicalDataForPlot$Age[TFindex]  <- "กร 60"
clinicalDataForPlot$Age[!TFindex] <- "< 60"
clinicalDataForPlot$Age           <- as.factor(clinicalDataForPlot$Age)

clinicalDataForPlot$`RUNX1 Fusion` <- factor(clinicalDataForPlot$`RUNX1 Fusion`, levels=c("No", "Yes"))
clinicalDataForPlot$`RUNX1 Mutation`       <- factor(clinicalDataForPlot$`RUNX1 Mutation`, levels=c("No", "Yes"))

Inx_alive                                   <- which(clinicalDataForPlot$Vital_Status == "alive")
Inx_dead                                    <- which(clinicalDataForPlot$Vital_Status == "dead")
clinicalDataForPlot$Vital_Status[Inx_alive] <- 0
clinicalDataForPlot$Vital_Status[Inx_dead]  <- 1
clinicalDataForPlot$Vital_Status            <- as.numeric(clinicalDataForPlot$Vital_Status)

# Hazard Ratio
model    <- coxph(Surv(OS.time, Vital_Status) ~ Age+`RUNX1 Fusion`+`RUNX1 Mutation`+LSCA_score, clinicalDataForPlot)
subFP    <- "/plot/forest_plot/GSE37642_LSCA_score_DEGs50.tiff"
resultFP <- paste0(projectDir, subFP)
tiff(filename=resultFP, width=12, height=9.6, res=300, units="in")
forest_model(model, limits=c(log(0.1), log(7)),
             format_options=forest_model_format_options(text_size=9)) + 
theme(axis.text.x=element_text(size=15, colour="black"))
dev.off()