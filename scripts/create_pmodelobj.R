library(PharmacoGx)
library(caret)
create_pmodelobj <- function(root, dataset_name, model_type,train_dataset_expr_matrix){
  ####note parameter train_dataset_expr_matrix for elasticnet is the gene expression matrix ofthe training pharmacogenomic dataset
  if (model_type == "elasticnet"){
    if(missing(train_dataset_expr_matrix)){
      print("error: please provide the training gene expression matrix")
      return(1)
    }
    my_working_directory <- root
    setwd(my_working_directory)
    source("get_clinical_data.R")
    
    my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),dataset_name)
    ex <- my_data$ex
    #nrow(ex)
    new_ex <- t(as.matrix(ex))
    #ncol(new_ex)
    colnames(new_ex)
    nrow(new_ex)
    res <- my_data$res
    meta <- my_data$metadata
    x <- train_dataset_expr_matrix
    i = 0
    inter <- c()
    for (n in colnames(x)){
      if(n %in% colnames(new_ex)){
        #print(n)
        inter <- append(inter, n)
        i = i + 1
      }
    }
    mymodel <- list(feature_genes = inter)
    class(mymodel)<-"elasticnet"
    return (mymodel)
  }
  
  else if(model_type == "SELECT"){
    my_working_directory <- root
    setwd(my_working_directory)
    source("get_CSL_network.R")
    csl <- get_CSL_network("CSL_pairs.csv")
    
    source("get_clinical_data.R")
    my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),dataset_name)
    ex <- my_data$ex
    #nrow(ex)
    new_ex <- t(as.matrix(ex))
    
    res <- my_data$res
    meta <- my_data$metadata
    
    #source("Drug_Targets.R")
    
    targets <- get_targets(meta$target_genes)
    #targets
    csl_set <- c()
    for (target in targets){
      if(target %in% csl$Gene){
        pairs <- get_CSL_pairs(csl, target)
        csl_set <- c(csl_set, pairs)
      }
    }
    csl_set <- csl_set[!duplicated(csl_set)]
    csl_set <-  intersect(colnames(new_ex), csl_set)
    #total_sl_pairs <- length(csl_set)
    mymodel <- list(feature_genes = csl_set)
    class(mymodel)<-"SELECT"
    return (mymodel)
    
    
    
    
  }
}

#create_pmodelobj <- function(exprmat, response, metadata){
  
#}
