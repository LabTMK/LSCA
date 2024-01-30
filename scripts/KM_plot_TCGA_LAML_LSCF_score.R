library(survival)
library(survminer)
library(dplyr)

projectDir <- "C:/Users/informatics3/Documents/GitHub/LSCA"


# cell type fractions - TCGA
subFP             <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_tcga_laml_RPKM.txt"
ctFrations_FP     <- paste0(projectDir, subFP)
ctFrations        <- read.delim2(ctFrations_FP, check.names=FALSE)
ctFrations        <- ctFrations[, 1:10]
ctFrations[,2:10] <- apply(ctFrations[,2:10], 2, as.numeric)


# LSCA score calculation
CellType    <- c("GMP", "CMP", "RApos", "MEP", "MPP") # 50 DEGs
Coefficient <- c(-2.1480989, -1.6350980, 0.3733158, 0.4942520, 4.5197050)

newScoreCoef             <- cbind(CellType, Coefficient)
newScoreCoef             <- as.data.frame(newScoreCoef)
newScoreCoef$Coefficient <- as.numeric(newScoreCoef$Coefficient)

fractions <- ctFrations[, CellType]
newScore  <- c()
for(i in 1:nrow(fractions)){
  temp     <- sum(fractions[i, ] * newScoreCoef$Coefficient)
  newScore <- c(newScore, temp)
}
Group               <- ifelse(newScore > median(newScore), "High", "Low")
bcr_patient_barcode <- ctFrations$Mixture
Table_TCGA_LAML     <- cbind(bcr_patient_barcode, Group)
Table_TCGA_LAML     <- as.data.frame(Table_TCGA_LAML)


# survival data
tcgaSurvivalFP <- paste0(projectDir, "/survival_data/TCGA_CDR.txt")
tcgaSurvival   <- read.delim2(tcgaSurvivalFP)

tcgaLAML              <- tcgaSurvival %>% select(bcr_patient_barcode, OS.time, vital_status)
tcgaLAML              <- Table_TCGA_LAML %>% left_join(tcgaLAML, by="bcr_patient_barcode")
tcgaLAML$OS.time      <- as.numeric(tcgaLAML$OS.time)
tcgaLAML$vital_status <- ifelse(tcgaLAML$vital_status == "Alive", 0, 1)

surv_object <- Surv(time=tcgaLAML$OS.time, event=tcgaLAML$vital_status)
fit         <- survfit(surv_object ~ Group, data=tcgaLAML)
pValue      <- surv_pvalue(fit = fit, data = tcgaLAML)
pValue      <- format(pValue$pval, scientific=T, digits=2)

plotTitle <- "TCGA-LAML, LSCA"
ggsurv    <- ggsurvplot(title=plotTitle, fit, data=tcgaLAML, xlab="Days", 
                        pval=paste0("Log-rank\np = ", pValue), pval.coord=c(1900,0.65), pval.size=10,
                        legend.title="", legend.labs=c("High", "Low"), palette=c("maroon", "forestgreen"),
                        risk.table=TRUE, risk.table.col="strata", conf.int = TRUE, risk.table.fontsize=7,
                        ggtheme=theme(plot.title=element_text(size=25, hjust=0.5, face="bold"),
                                      axis.text=element_text(size=20, color="black"),
                                      legend.text=element_text(size=25, color="black"),
                                      axis.title=element_text(size=25, color="black"),
                                      panel.background=element_rect(fill="white", colour="grey50"),
                                      panel.grid=element_line(colour="grey")))
resultFP <- paste0(projectDir, "/plot/KM_plot/KM_Plot_TCGA_LAML_LSCA_score_DEGs50.tiff")
tiff(filename=resultFP, width=17.78, height=17.78, res=300, units="cm")
ggsurv
dev.off()