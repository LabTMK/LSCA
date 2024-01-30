library(dplyr)
library(glmnet)

projectDir  <- "C:/Users/informatics3/Documents/GitHub/LSCA"


# cell type fractions - GSE37642
subFP          <- "/CIBERSORTx_results/CIBERSORTx_HemLin9_DEGs50_GSE37642_GPL96.txt"
ctFractions_FP <- paste0(projectDir, subFP)
ctFractions    <- read.delim2(ctFractions_FP, check.names=FALSE)

# remove P-value, correlation, RSME
ctFractions <- ctFractions[,-c(11,12,13)]

colnames(ctFractions)[1] <- "GEO_ID"


clinicalData_FP <- paste0(projectDir, "/survival_data/GSE37642_Survival_data.txt")
clinicalData    <- read.delim2(clinicalData_FP)
clinicalData    <- clinicalData %>% select(GEO_ID, overall_survival_days, life.status)
gse37642_gpl96  <- ctFractions %>% left_join(clinicalData, by="GEO_ID")

gse37642_gpl96$life.status <- ifelse(gse37642_gpl96$life.status == "alive", 0, 1)

gse37642_gpl96 <- gse37642_gpl96[, -1]
gse37642_gpl96 <- apply(gse37642_gpl96, 2, as.numeric)
gse37642_gpl96 <- as.data.frame(gse37642_gpl96)

X <- as.matrix(gse37642_gpl96[,1:9])
Y <- as.matrix(gse37642_gpl96[,c(10,11)])

TFindex     <- !(is.na(Y[,1])) & Y[,1] != 0
X           <- X[TFindex, ]
Y           <- Y[TFindex, ]
colnames(Y) <- c("time", "status")


#Penalty type (alpha=1 is lasso and alpha=0 is the ridge)
set.seed(123)
coefficientsMatrix <- matrix(nrow=9)
for(i in 1:100){
  cv.lambda.lasso    <- cv.glmnet(x=X, y=Y, family="cox", alpha=1) 
  l.lasso.min        <- cv.lambda.lasso$lambda.min
  lasso.model        <- glmnet(x=X, y=Y, family="cox", alpha=1, lambda=l.lasso.min)
  tempCoef           <- lasso.model$beta #Coefficients
  coefficientsMatrix <- cbind(coefficientsMatrix, tempCoef)
}
coefficientsMatrix   <- as.matrix(coefficientsMatrix)
coefficientsMatrix   <- coefficientsMatrix[, -1]

# Trimming Coefficients Matrix
significantCTinx <- c()
for(i in 1:nrow(coefficientsMatrix)){
  NumOfNonZero <- sum(coefficientsMatrix[i, ] != 0 )
  if(NumOfNonZero >= 95){
    significantCTinx <- c(significantCTinx, i)
  }
}
coefMatrixTrimmed <- coefficientsMatrix[significantCTinx, ]

coefMean <- apply(coefMatrixTrimmed, 1, mean)
coefMean <- as.data.frame(sort(coefMean))
coefMean <- cbind(rownames(coefMean), coefMean)

colnames(coefMean) <- c("CellType", "Coefficient")
