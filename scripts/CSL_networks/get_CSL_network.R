
library(stringr)

get_CSL_network <- function(path){
  ### Takes in a path to the csv file and generates a dataframe representing the network
  ### The dataframe returned has 2 columns: "Gene" and "CSL_pairs". for each gene name in the "Gene" column, 
  ### the entry in "CSL_pairs" is a string of all the CSL pair genes, with each gene name separated by a space (" ")
  
  df <- read.csv(path)
  
  
  
  csl <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Gene", "CSL_pairs"))#(Gene = character(), SL_pairs = character())
  
  for (row in 1:nrow(df)) {
     x <- df[row, "name"]
     gene1 <- sapply(strsplit(x," "), getElement, 1)
     gene2 <- sapply(strsplit(x," "), getElement, 4)
   
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
  
  return(csl)
  
}





get_CSL_pairs <- function(df_network, gene_name){
  
  ### Helper function to process the CSL pair data from the dataframe
  ### takes in the dataframe representing the CSL network (generated with get_CSL_network) and a gene name
  ### for the given gene name, will return the CSL pairs in a vector, with each entry representing the gene name of a CSL pair
  
  pairs <- df_network[which(df_network$Gene == gene_name), 2]
  sep <- strsplit(pairs, " ")
  vec <- c()
  for (word in sep){
    vec <- append(vec, word)
    
  }
  return(vec)
  
}



