library(PharmacoGx)

# ==== helper functions ====

train_models <- function(x, y, mydrug){
  # genes <- rownames(x)[1:1000]
  genes <- NULL
    
  elnmodel <- pmodel_train(x, y, mydrug, modeltype="elasticnet", folds=5, hyperparams="", genelist=genes)
  
  rfmodel <- pmodel_train(x, y, mydrug, modeltype="rf", folds=5, hyperparams="", genelist=genes)
  
  return(list(elnmodel=elnmodel, rfmodel=rfmodel))
}

# ==== prepare data ====

gdsc <- downloadPSet("GDSC_2020(v2-8.2)")

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

model <- train_models(fdata, aacs, "Elephantin")
# alpha     lambda       RMSE  Rsquared        MAE      RMSESD RsquaredSD
# 1   0.2 0.04266578 0.06962304 0.1584558 0.05613209 0.007625451 0.07079398
# MAESD
# 1 0.00529516
# mtry       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 19957 0.06914015 0.1711985 0.05529098 0.006877269 0.06823276 0.004974111

model2 <- train_models(fdata, aacs, "Sirolimus")
# alpha     lambda       RMSE  Rsquared        MAE      RMSESD RsquaredSD
# 1   0.2 0.05717409 0.09400033 0.1865315 0.07516952 0.008769055 0.03192125
# MAESD
# 1 0.006553272
# mtry       RMSE  Rsquared        MAE      RMSESD RsquaredSD       MAESD
# 1 1996 0.09332013 0.2117538 0.07521546 0.009566521   0.053586 0.007032985

model3 <- train_models(fdata, aacs, "Irinotecan")
# alpha     lambda       RMSE Rsquared        MAE     RMSESD RsquaredSD
# 1   0.1 0.09159966 0.08702908 0.434758 0.06584657 0.01063606  0.1091689
# MAESD
# 1 0.007514599
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 19957 0.08591522 0.4540355 0.06519536 0.01038419  0.1142583 0.006811519

model4 <- train_models(fdata, aacs, "Nelarabine")
# alpha     lambda       RMSE  Rsquared        MAE     RMSESD RsquaredSD
# 1   0.4 0.02355453 0.04738242 0.0668078 0.02964353 0.01095683 0.07134862
# MAESD
# 1 0.001856583
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1    2 0.04748694 0.0430369 0.03107479 0.01081124 0.04994968 0.001312485
Sys.time()

model5 <- train_models(fdata, aacs, "Sorafenib")
# alpha     lambda       RMSE  Rsquared       MAE     RMSESD RsquaredSD
# 1   0.2 0.05275609 0.07188754 0.2771664 0.0435744 0.02810783  0.1048557
# MAESD
# 1 0.007621864
# mtry       RMSE  Rsquared        MAE     RMSESD RsquaredSD       MAESD
# 1 19957 0.06573257 0.3674034 0.04164824 0.02705055 0.05538569 0.008016982

# models_compare <- resamples(list(EN=model$elnmodel, RF=model$rfmodel))

# Summary of the models performances
# summary(models_compare)
