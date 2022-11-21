library(PharmacoGx)

# ==== helper functions ====

train_models <- function(x, y, mydrug, modeltype="elasticnet"){
  # genes <- rownames(x)[1:1000]
  genes <- NULL
  elnmodel <- NULL
  rfmodel <- NULL
  
  if (modeltype == "rf"){
    rfmodel <- pmodel_train(x, y, mydrug, modeltype="rf", folds=5, hyperparams="", genelist=genes)
  } else if (modeltype == "all"){
    elnmodel <- pmodel_train(x, y, mydrug, modeltype="elasticnet", folds=5, hyperparams="", genelist=genes)
    rfmodel <- pmodel_train(x, y, mydrug, modeltype="rf", folds=5, hyperparams="", genelist=genes)
  } else {
    elnmodel <- pmodel_train(x, y, mydrug, modeltype="elasticnet", folds=5, hyperparams="", genelist=genes)
  }
  
  return(list(elnmodel=elnmodel, rfmodel=rfmodel))
}

# ==== prepare data ====

# gdsc <- downloadPSet("GDSC_2020(v2-8.2)")
gdsc <- readRDS("GDSC2.rds")
gdsc <- updateObject(gdsc)
mysample <- sampleInfo(gdsc)

# Use aac_recomputed by default
# Rows are drugs, columns are cell lines
aacs <- summarizeSensitivityProfiles(gdsc, sensitivity.measure="aac_recomputed")

# Get gene expression data - maybe use microarray instead of RNA-seq
mps <- summarizeMolecularProfiles(gdsc, mDataType="Kallisto_0.46.1.rnaseq")

fdata <- assay(mps)

mygeneinfo <- featureInfo(gdsc, mDataType = "Kallisto_0.46.1.rnaseq")

# Slice gene expression matrix to just protein coding genes
# This is gross, clean up
fdata <- fdata[rownames(fdata) %in% rownames(mygeneinfo)[mygeneinfo$gene_type == "protein_coding"], ]
rownames(fdata) <- mygeneinfo$gene_name[match(rownames(fdata), rownames(mygeneinfo))]

# ==== train model ====

setwd("~/Desktop/MBP_PhD/3rd_rotation/synthLethal/R")
# Plug into model training
source("pmodel_train.R")

model <- train_models(fdata, aacs, "Elephantin", "all")
# alpha     lambda       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 0.325 0.03009987 0.06965862 0.1564548 0.05589782 0.007589288 0.06587157 0.005186235
# mtry       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 19957 0.06914015 0.1711985 0.05529098 0.006877269 0.06823276 0.004974111

model2 <- train_models(fdata, aacs, "Sirolimus", "all")
# alpha     lambda       RMSE Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 0.325 0.04033519 0.09402407 0.186458 0.07509004 0.009314951  0.0357338 0.006781816
# mtry       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 1996 0.09332013 0.2117538 0.07521546 0.009566521   0.053586 0.007032985

model3 <- train_models(fdata, aacs, "Irinotecan", "all")
# alpha     lambda       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 0.325 0.02999476 0.08698538 0.4327573 0.06538253 0.01107834  0.1130487 0.007596179
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 19957 0.08591522 0.4540355 0.06519536 0.01038419  0.1142583 0.006811519

model4 <- train_models(fdata, aacs, "Nelarabine", "all")
# alpha     lambda       RMSE   Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1  0.55 0.01661725 0.04746451 0.06764595 0.02968279 0.01090653 0.07129619 0.001975702
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1    2 0.04748694 0.0430369 0.03107479 0.01081124 0.04994968 0.001312485
Sys.time()

model5 <- train_models(fdata, aacs, "Sorafenib", "all")
# alpha     lambda       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 0.325 0.03721838 0.07243965 0.2636793 0.04390854 0.02865489  0.1027081 0.007780694
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 19957 0.06573257 0.3674034 0.04164824 0.02705055 0.05538569 0.008016982

models_compare <- resamples(list(Elephantin_EN=model$elnmodel, Elephantin_RF=model$rfmodel, Sirolimus_EN=model2$elnmodel, Sirolimus_RF=model2$rfmodel, Irinotecan_EN=model3$elnmodel, IrinotecanRF=model3$rfmodel, NelarabineEN=model4$elnmodel, NelarabineRF=model4$rfmodel, Sorafenib_EN=model5$elnmodel, Sorafenib_RF=model5$rfmodel))

# Summary of the models performances
summary(models_compare)

# ==== Train on each tissue type independently with just expression ====

# here we only use lung cell lines
fdata_lung <- fdata[ , which(mysample$unique.tissueid.fromstudies == "lung")]
aacs_lung <- aacs[ , which(mysample$unique.tissueid.fromstudies == "lung")]

model_lung <- train_models(fdata_lung, aacs_lung, "Cisplatin")
# alpha     lambda       RMSE  Rsquared        MAE     RMSESD RsquaredSD      MAESD
# 1 0.775 0.04669098 0.05931393 0.4399339 0.04685312 0.02494641  0.3244823 0.01108179
model_no_lung <- train_models(fdata, aacs, "Cisplatin")
# alpha     lambda      RMSE  Rsquared      MAE      RMSESD RsquaredSD       MAESD
# 1   0.1 0.03060134 0.0596531 0.2755406 0.043068 0.009680539 0.07878425 0.006280853

models_compare_tissue_specific <- resamples(list(lung_specific=model_lung$elnmodel, all_tissue=model_no_lung$elnmodel))
summary(models_compare_tissue_specific)

