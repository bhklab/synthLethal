library(caret)
library(pROC)
my_working_directory <- "C:/Users/andyl/Documents/EngSci/Year 1/Medical Biophysics/Work/SELECT"
setwd(my_working_directory)
source("get_CSL_network.R")
csl <- get_CSL_network("CSL_pairs.csv")
drug <- "Sorafenib"

source("get_clinical_data.R")
my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),"GSE119262")
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


test.unique <- function(lis) {  ## function to test unique columns
  
  length1 <- length(lis)
  length2 <- length(unique(lis))        
  if (length1 - length2 > 0 ) {
    
    print(paste("There are", length1 - length2, " duplicates", sep=" "))
  } 
  else{
    print("no repeating genes")
  }
}
length(csl_set)
test.unique(csl_set)
csl_set <- csl_set[!duplicated(csl_set)]
csl_set

csl_set <-  intersect(colnames(new_ex), csl_set)
#length(intersection)
###intersection
#csl_set <- csl_set[sapply(csl_set, function(x) any(words %in% x))]
total_sl_pairs <- length(csl_set)
csl_set

#print("KIT" %in% colnames(new_ex))

#new_ex[,which(colnames(new_ex)=="BRAF")]

thresholds <- matrix(matrix(ncol =1, nrow = total_sl_pairs))
rownames(thresholds) <- csl_set
colnames(thresholds) <- c("expression_threshold")
  #setNames(data.frame(matrix(ncol = 1, nrow = total_sl_pairs )), ,csl_set)
#length(rownames(new_ex))
index <-  round(length(rownames(new_ex)) * 1/3, 0)
#index

#thresholds[which(rownames(thresholds) == "DDOST"),1 ] <- 23
#thresholds
#test.unique(colnames(new_ex))
#a <- sort(new_ex[index,which(colnames(new_ex)=="SRRT")[1]])
#print(a)
#NCOL(new_ex[,which(colnames(new_ex)=="KIT")])
#NCOL(new_ex[,which(colnames(new_ex)=="BRAF")])
#print(sum(is.na(new_ex[,which(colnames(new_ex)=="SRRT")])))
#print(is.na(new_ex[1,which(colnames(new_ex)=="SRRT")]))
#colSums(is.na(df))<nrow(df)

#length(csl_set)
j = 0
for (gene in csl_set){
  #print(gene)
  val <- sort(new_ex[index,which(colnames(new_ex)==gene)][1])
  
  #print(val)
  #print(thresholds[which(rownames(thresholds) == gene),1 ])
  #print(val)
  #print(gene)
  thresholds[which(rownames(thresholds) == gene),1 ] <- val
  #print(thresholds[which(rownames(thresholds) == gene),1 ])
  j = j +  1
  print(j)
  #print(thresholds[which(rownames(thresholds) == gene),1 ])
}

#thresholds["AUP1", 1]
#NROW(thresholds)
#new_ex[which(rownames(new_ex)=="GSM2935412"), "AUP1"]
#rownames(new_ex)

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


pred <- matrix(nrow = length(rownames(new_ex)), ncol=1)
rownames(pred)<-rownames(new_ex)

colnames(pred) <- c("SL_score")

obs <- matrix(nrow = length(rownames(new_ex)), ncol=1)
rownames(obs)<-rownames(new_ex)

#colnames(pred) <- c("SL_score")

for (sample in rownames(pred)){
  pred[which(rownames(pred) == sample), 1] <- calculate_SL_score(new_ex, sample, csl_set, thresholds)
  if(res[which(res$sample == sample), "response"] == 1){  
    obs[which(rownames(obs) == sample), 1] <- 1
  } else{
    obs[which(rownames(obs) == sample), 1] <- 0
  }
}

obs <- factor(obs)
SELECT_eval <- roc(obs, pred)


print(SELECT_eval)
plot(SELECT_eval)

if (FALSE){ 
  setwd("C:/Users/andyl/Documents/EngSci/Year 1/Medical Biophysics/Work/SELECT")

  source("get_TCGA_Drug_CSL_network.R")
  net <- get_TCGA_Drug_CSL_network("TCGA_Drug_response_CSL_pairs.csv")
  length(net$targets$Drug)
  print( net$targets$Drug )


  my_targs <- get_TCGA_Drug_Targets(net$targets,"Celecoxib" )
  my_targs
  tcga_csl <- c()
  for (name in my_targs){
    print(name)
    pairs <- get_CSL_pairs(net$csl, name)
    tcga_csl <- c(tcga_csl, pairs)
  }
  length(tcga_csl)
  tcga_csl
  
  
  net2 <- get_CSL_network("sr.final.v13.symbol2.2326.txt default edge.csv")
  csl_set2 <- c()
  for (target in targets){
    if(target %in% net2$Gene){
      pairs <- get_CSL_pairs(net2, target)
      csl_set2 <- c(csl_set2, pairs)
    }
  }
  length(csl_set2)
  test.unique(csl_set2)
  print(csl_set2)
}
