library(caret)
library(pROC)
library(stringr)
#######GETTING THE TCGA PHARMACOGENOMIC DATA
gdsc <- downloadPSet("GDSC_2020(v2-8.2)")

#ctrp_drugs <- drugInfo(CTRPv2)
#ctrp_drugs

#ctrp_aacs <- summarizeSensitivityProfiles(CTRPv2, sensitivity.measure="aac_recomputed")

#??PharmacoSet

# Get drug info
mydrugs <- drugInfo(gdsc)

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


mydrug <- "Sorafenib"
mytarget <- "BRAF"

if (length(grep(mydrug, rownames(aacs), ignore.case=TRUE)) != 1){
  print("Warning: number of drugs matched does not equal 1.")
}

yall <- aacs[grep(mydrug, rownames(aacs), ignore.case=TRUE), ]

# Remove NAs - cell lines for which the drug was not assayed, or cell lines without rna-seq data
missingcells <- colSums(is.na(fdata))

x <- fdata[, (complete.cases(yall) & missingcells == 0)]
y <- yall[(complete.cases(yall) & missingcells == 0)]
x<-t(x)







#######Training and evaluating Drug Targets Regression Model on the GSE109211 Sorafenib dataset
##########CHANGE YOUR WORKING DIRECTORY TO WHEREVER YOUR R FILES ARE

my_working_directory <- "C:/Users/andyl/Documents/EngSci/Year 1/Medical Biophysics/Work/SELECT"
setwd(my_working_directory)
source("get_clinical_data.R")
my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),"GSE119262")
ex <- my_data$ex
#nrow(ex)
new_ex <- t(as.matrix(ex))

res <- my_data$res
meta <- my_data$metadata

get_targets <- function(targets){
  
  #parses strings separated by commas to a vector of drug targets
  
  
  sep <- strsplit(targets, ",")
  vec <- c()
  for (word in sep){
    vec <- append(vec, word)
    
  }
  return(vec)
  
}
targets <- get_targets(meta$target_genes)
length(targets)
targets
intersection <-  intersect(colnames(x), targets)
intersection
x_new <- x[, intersection, drop=FALSE]

#print("MTOR" %in% colnames(x))
#rownames(x_new)
#colnames(x_new)
#head(x_new)
#NCOL(x_new)
#NROW(x_new)


#nrow(x_new)

set.seed(42)
cv_5 <- trainControl(method = "cv", number = 5)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}
drug_target_model <- train(
  x_new, y,
  method = "lm",
  trControl = cv_5,
  tuneLength = 10
)
get_best_result(drug_target_model)
#"MTOR" %in% colnames(new_ex)
new_test <-new_ex[, intersection]
ncol(new_test)

x_hat_pre <- predict(drug_target_model, new_test)
x_hat_pre
print(colnames(x_new) == colnames(new_test))
#NCOL(x_hat_pre)
#x_hat_old <- predict(elnet, new_test)###neror number of genes must match
sample_names <- res$sample
sample_names
gt_mat <- cbind(res$response, res$no_response)
colnames(gt_mat) <- c("response", "no_response")
rownames(gt_mat)<- sample_names
head(gt_mat)

pred <- x_hat_pre
obs <- c()
#pred <- c()
for(i in 1:NROW(x_hat_pre)){
  sample_name <- ROWNAMES(x_hat_pre[i])
  
  if(gt_mat[i, "response"] == 1){
    obs <- append(obs, 1)
    
  }
  else{
    obs <- append(obs, 0)
  }
  
}

obs <- factor(obs)

sorafenib_drug_target_eval <- roc(obs, pred)


print(sorafenib_drug_target_eval)
plot(sorafenib_drug_target_eval)
