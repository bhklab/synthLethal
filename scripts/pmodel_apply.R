library(caret)
library(pROC)

calculate_SL_score <- function(ex, sample, csl_pairs, thresholds){
  tot <- length(csl_pairs)
  inactive <- 0
  for (pair in csl_pairs){
    thresh <- thresholds[pair, 1]
    expr <- ex[which(rownames(ex) == sample), pair]
    if(expr <= thresh){
      inactive <- inactive + 1
    }
  }
  return(inactive / tot)
}


pmodel_apply <- function(expr_mat, response_vector, mymodel, drug="", modeltype=""){
  
  ###expr_mat is the expression matrix of the dataset (target clinical trial)
  
  new_ex <- t(as.matrix(expr_mat)) ###expects "ex" dataframe from get_clinical_data function
  
  if(modeltype == "elasticnet"){
    inter <- mymodel$feature_genes
    new_test <-new_ex[, inter]
    x_hat_pre <- predict(mymodel$model, new_test)
    #print(x_hat_pre)
    res <- response_vector
    sample_names <- res$sample
    #length(sample_names)
    gt_mat <- cbind(res$response, res$no_response)
    colnames(gt_mat) <- c("response", "no_response")
    rownames(gt_mat)<- sample_names
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
    
    eval <- roc(obs, pred)
    print(eval)
    plot(eval)
    
    out <- list()
    out$eval <- eval
    out$preictions <- pred
    return(out)
    
    
  }
  
  
  else if(modeltype == "SELECT"){
    csl_set <- mymodel$feature_genes
    total_sl_pairs <- length(csl_set)
    thresholds <- matrix(matrix(ncol =1, nrow = total_sl_pairs))
    rownames(thresholds) <- csl_set
    colnames(thresholds) <- c("expression_threshold")
    index <-  round(length(rownames(new_ex)) * 1/3, 0)
    
    ###getting the threshold expression (bottom tertile from SELECT)
    j = 0
    for (gene in csl_set){

      val <- sort(new_ex[index,which(colnames(new_ex)==gene)][1])
      
      
      thresholds[which(rownames(thresholds) == gene),1 ] <- val
      
      j = j +  1
    }
    
    pred <- matrix(nrow = length(rownames(new_ex)), ncol=1)
    rownames(pred)<-rownames(new_ex)
    
    colnames(pred) <- c("SL_score")
    
    obs <- matrix(nrow = length(rownames(new_ex)), ncol=1)
    rownames(obs)<-rownames(new_ex)
    
    for (sample in rownames(pred)){
      pred[which(rownames(pred) == sample), 1] <- calculate_SL_score(new_ex, sample, csl_set, thresholds)
      if(res[which(res$sample == sample), "response"] == 1){  
        obs[which(rownames(obs) == sample), 1] <- 1
      } else{
        obs[which(rownames(obs) == sample), 1] <- 0
      }
    }
    
    obs <- factor(obs)
    eval <- roc(obs, pred)
    
    
    print(eval)
    plot(eval)
    out <- list()
    out$eval <- eval
    out$preictions <- pred
    return(out)
    
  }
    
}