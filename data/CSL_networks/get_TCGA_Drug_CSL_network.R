
library(stringr)

get_TCGA_Drug_CSL_network <- function(path){
  ### Takes in a path to the TCGA Drugs CSL csv file and generates a dataframe representing the network
  ### The dataframe returned has 2 columns: "Gene" and "CSL_pairs". for each gene name in the "Gene" column, 
  ### the entry in "CSL_pairs" is a string of all the CSL pair genes, with each gene name separated by a space 
  
  ### ALSO returns a dataframe with 2 columns: "Drug" and "Target_genes". Each drug in TCGA has its corresponding target genes in the "Target_genes" column, with each gene name separated by a space (" ").
  
  ### This function returns the 2 dataframes in a list
  
  df <- read.csv(path)
  drugs <- c("Letrozole", "Exemestane", "Anastrozole", "Camptothecin", "Celecoxib", "Trastuzumab", "Etoposide", "Doxorubicin", "Epirubicin", "Gemcitabine", "Cetuximab", "Erlotinib", "Dacarbazine", "Vinorelbine", "Bevacizumab", "Taxane")
  
  
  csl <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Gene", "CSL_pairs"))#(Gene = character(), SL_pairs = character())
  drug_targets <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Drug", "Target_genes"))
  for (row in 1:nrow(df)) {
    x <- df[row, "name"]
    gene1 <- sapply(strsplit(x," "), getElement, 1)
    gene2 <- sapply(strsplit(x," "), getElement, 4)
    if (gene1 %in% drugs){
      if (gene1 %in% drug_targets$Drug){
        targets <- drug_targets[which(drug_targets$Drug == gene1), 2]
        targets <- paste(targets, gene2, sep=" ")
        drug_targets[which(drug_targets$Drug == gene1), 2] <- targets
      }
      else{
        
        df2<-data.frame(gene1,gene2)
        names(df2)<-c("Drug","Target_genes")
        drug_targets <- rbind(drug_targets, df2)
      }
    }
    else{
      if (gene1 %in% csl$Gene){
        pairs <- csl[which(csl$Gene == gene1), 2]
        pairs <- paste(pairs, gene2, sep=" ")
        csl[which(csl$Gene == gene1), 2] <- pairs
      }
      else{
        
        df2<-data.frame(gene1,gene2)
        names(df2)<-c("Gene","CSL_pairs")
        csl <- rbind(csl, df2)
      }
      if (gene2 %in% csl$Gene){
        pairs <- csl[which(csl$Gene == gene2), 2]
        pairs <- paste(pairs, gene1, sep=" ")
        csl[which(csl$Gene == gene2), 2] <- pairs
      }
      else{
        df2<-data.frame(gene2, gene1)
        names(df2)<-c("Gene","CSL_pairs")
        csl <- rbind(csl, df2)
      }
    }
  }
  #return(drug_targets)
  out <- list()
  out$csl <- csl
  out$targets <- drug_targets
  return(out)#list(csl, drug_targets))
  
}





get_TCGA_Drug_Targets <- function(df_network, drug_name){
  
  ### Helper function to process the Drug Target data from the dataframe
  ### takes in the dataframe representing the Drugs and their target genes (generated with get_TCGA_Drug_CSL_network) and a drug name
  ### for the given gene name, will return the target gene names in a vector, with each entry representing the gene name of a target gene
  ### To extract the CSL pairs, the function get_CSL_pairs defined in get_CSL_network.R can be used
  targets <- df_network[which(df_network$Drug == drug_name), 2]
  sep <- strsplit(targets, " ")
  vec <- c()
  for (word in sep){
    vec <- append(vec, word)
    
  }
  return(vec)
  
}

