library(PharamcoGx)
library(caret)
get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}
pmodel_train <- function(expr_mat, drug_response, drug, modeltype="elasticnet", drugTarget="", root, target_dataset_name){
  #target_dataset_name = name of the clinical trial dataset you will test the model on
  # Filter out NAs
  
  gdsc <- downloadPSet("GDSC_2020(v2-8.2)")
  if (modeltype == "elasticnet"){
    # Use aac_recomputed by default
    # Rows are drugs, columns are cell lines
    # <- summarizeSensitivityProfiles(gdsc, sensitivity.measure="aac_recomputed")
    
    
    # Get gene expression data - maybe use microarray instead of RNA-seq
    #mps <- summarizeMolecularProfiles(gdsc, mDataType="Kallisto_0.46.1.rnaseq")
    
    #fdata <- assay(mps)
    #mygeneinfo <- featureInfo(gdsc, mDataType = "Kallisto_0.46.1.rnaseq")
    # Slice gene expression matrix to just protein coding genes
    # This is gross, clean up
    #fdata <- fdata[rownames(fdata) %in% rownames(mygeneinfo)[mygeneinfo$gene_type == "protein_coding"], ]
    #rownames(fdata) <- mygeneinfo$gene_name[match(rownames(fdata), rownames(mygeneinfo))]
    aacs <- drug_response
    
    
    mydrug <- drug
    mytarget <- drugTarget
    
    if (length(grep(mydrug, rownames(aacs), ignore.case=TRUE)) != 1){
      print("Warning: number of drugs matched does not equal 1.")
    }
    
    #yall <- aacs[grep(mydrug, rownames(aacs), ignore.case=TRUE), ]
    
    # Remove NAs - cell lines for which the drug was not assayed, or cell lines without rna-seq data
    #missingcells <- colSums(is.na(fdata))
    
    #x <- fdata[, (complete.cases(yall) & missingcells == 0)]
    #y <- yall[(complete.cases(yall) & missingcells == 0)]
    y <- aacs
    x <- t(expr_mat)

    set.seed(42)
    cv_5 <- trainControl(method = "cv", number = 5)
    source("create_pmodelobj.R")
    elnet_model <- create_pmodelobj(root, target_dataset_name, 'elasticnet', x)
    inter <- elnet_model$feature_genes
    
    x_new <- x[, inter]
    elnet <- train(
      x_new, y,
      method = "glmnet",
      trControl = cv_5,
      tuneLength = 10
    )
    print(get_best_result(elnet))
    elnet_model$model <- elnet
    return(elnet_model) #S3 class object for elasticnet
  } else if (modeltype == "univariate"){
    
  } else if (modeltype == "drugtarget"){
#not trainable
  } else if (modeltype == "SELECT") {
    #should be pretrained???
    #SELECT_model <- create_pmodelobj(root, target_dataset_name, 'SELECT', NULL)
    #csl_set <- SELECT_model$feature_genes
    #total_sl_pairs <- length(csl_set)
    #thresholds <- matrix(matrix(ncol =1, nrow = total_sl_pairs))
    #rownames(thresholds) <- csl_set
    #colnames(thresholds) <- c("expression_threshold")
    
  } else if (modeltype == "random"){
   #not trainable 
  }
  
  
  return(mymodel)
}


