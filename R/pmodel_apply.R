library(caret)
library(pROC)

pmodel_apply <- function(mymodel, clinical_data){
  
  genes <- colnames(mymodel$trainingData)[1:length(mymodel$trainingData) - 1]
  
  x <- t(as.matrix(clinical_data$ex)[genes, ]) 

  x_hat_pre <- predict(mymodel, x)
  
  res <- clinical_data$res
  
  obs <- factor(res$response)
    
  eval <- roc(obs, x_hat_pre)
  print(eval)
  plot(eval)
    
  out <- list()
  out$eval <- eval
  out$preictions <- x_hat_pre
  return(out)
}