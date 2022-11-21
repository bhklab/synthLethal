library(PharmacoGx)
library(caret)

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

pmodel_train <- function(expr_mat, drug_response, drug, modeltype="elasticnet", folds=5, hyperparams="", genelist=NULL){
  
  if (!is.null(genelist)){
    inter <- intersect(rownames(expr_mat), genelist)
    if (length(inter) == 0){
      print("Error: No matched genes.")
      return(1)
    }
    x <- expr_mat[inter,]
  } else {
    x <- expr_mat
  }
  
  if (length(grep(drug, rownames(drug_response), ignore.case=TRUE)) != 1){
    print("Warning: number of drugs matched does not equal 1.")
    return(1)
  }
  
  yall <- drug_response[grep(drug, rownames(drug_response), ignore.case=TRUE), ]
  
  # Remove NAs - cell lines for which the drug was not assayed, or cell lines without rna-seq data
  missingcells <- colSums(is.na(x))
  
  x <- x[, (complete.cases(yall) & missingcells == 0)] 
  y <- yall[(complete.cases(yall) & missingcells == 0)]
    
  x <- t(x)
  
  set.seed(42)
  cv_5 <- trainControl(method = "cv", number = folds, savePredictions = "final")
    
  if (modeltype == "elasticnet"){
    mymodel <- train(
      x, y,
      method = "glmnet",
      trControl = cv_5,
      tuneLength = 10
    )
    
  } else if (modeltype == "rf"){
    mymodel <- train(
      x, y,
      method = "rf",
      trControl = cv_5,
      tuneLength = 5
    )
    
  } else {
    print("Error: modeltype is not valid.")
    return(1)
  }
  
  print(get_best_result(mymodel))
  return(mymodel)
}


