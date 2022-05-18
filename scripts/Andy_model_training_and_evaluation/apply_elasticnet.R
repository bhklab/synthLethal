library(PharmacoGx)
library(caret)
library(pROC)
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


##########Training and Applying Elasticnet on the GSE109211 Sorafenib dataset

set.seed(42)
cv_5 <- trainControl(method = "cv", number = 5)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}


#print(get_best_result(elnet))
##########CHANGE YOUR WORKING DIRECTORY TO WHEREVER YOUR R FILES ARE
my_working_directory <- "C:/Users/andyl/Documents/EngSci/Year 1/Medical Biophysics/Work/SELECT"
setwd(my_working_directory)
source("get_clinical_data.R")

my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),"GSE119262")
ex <- my_data$ex
#nrow(ex)
new_ex <- t(as.matrix(ex))
#ncol(new_ex)
colnames(new_ex)
nrow(new_ex)
res <- my_data$res
meta <- my_data$metadata

#colnames(new_ex)
#colnames(x)
#intersection <-  intersect(colnames(x), colnames(new_ex))
#length(intersection)
#print(intersection)
#print('ABCC9' %in% colnames(new_ex))
#ncol(x)
#ncol(new_ex)
#print(new_ex["GSM2935406",])

#names <- colnames(x)[(colnames(x) %in% intersection)]
#x_new <- x[, names]

#ncol(x_new)


test.unique <- function(x_new) {  ## function to test unique columns
  
  length1 <- length(colnames(df))
  length2 <- length(unique(colnames(df)))        
  if (length1 - length2 > 0 ) {
    
    print(paste("There are", length1 - length2, " duplicates", sep=" "))
  } 
  else{
    print("no repeating genes")
  }
}
test.unique(new_ex)
test.unique(x)
#length(intersection)
#names <- colnames(x)[(colnames(x) %in% intersection)]
#length(colnames(x))
i = 0
inter <- c()
for (n in colnames(x)){
  if(n %in% colnames(new_ex)){
    #print(n)
    inter <- append(inter, n)
    i = i + 1
  }
}
#ncol(x)
x_new <- x[, inter]


#print(j)
#length(colnames(x_new))
#test.unique(x_new)

#length(unique(intersection))
ncol(x_new)

elnet_restricted <- train(
  x_new, y,
  method = "glmnet", # I also tried "enet" and it produced identical results
  trControl = cv_5,
  tuneLength = 10
)
#length(intersection)

get_best_result(elnet_restricted)
new_test <-new_ex[, inter]
ncol(new_test)
#NCOL(new_ex["GSM2935406",])
#new_ex["GSM2935406",]

#elnet<- train(
#  x, y,
#  method = "glmnet",
#  trControl = cv_5,
#  tuneLength = 10
#)

x_hat_pre <- predict(elnet_restricted, new_test)
x_hat_pre
print(colnames(x_new) == colnames(new_test))
NCOL(x_hat_pre)
#x_hat_pre
      #x_hat_old <- predict(elnet, new_test)###neror number of genes must match
sample_names <- res$sample
length(sample_names)
gt_mat <- cbind(res$response, res$no_response)
colnames(gt_mat) <- c("response", "no_response")
rownames(gt_mat)<- sample_names
head(gt_mat)

#gt_mat
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

sorafenib_eval <- roc(obs, pred)


print(sorafenib_eval)
plot(sorafenib_eval)
      
