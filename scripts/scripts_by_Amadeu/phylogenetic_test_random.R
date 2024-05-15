#' Title
#'
#' @param SL_pairs 
#'
#' @return
#' 
#'
#' @examples
phylogenetic_test_random <- function(SL_pairs){



library(tidyverse)
library(tictoc)


datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"



phylo_path <- file.path(datapath, 'yuval.phylogenetic.profile.random.RData')
load(phylo_path)
### The phylogenetic profile is downloaded from Yuval Tabach et al. Mol Syst Biol. (2013), Supplementary Table 1
weights_path <- file.path(datapath, 'feature.weight.RData')
load(weights_path)
### the feature weights are determined based on the phylogenetic tree (Ensembl database: http://useast.ensembl.org/index.html)

phylo_coefficient <- vector(mode="numeric", length = nrow(SL_pairs))

p <- progressr::progressor(nrow(SL_pairs) %/% 1000)

tic()
for (i in  1:nrow(SL_pairs)){
  if(i %% 1000 == 0){
    p()
  }
  gene1 <- SL_pairs$gene1[i]
  gene2 <-  SL_pairs$gene2[i]
  if (gene1 %in% phylo$genes && gene2 %in% phylo$genes){
    
    sl.phylo =  cbind(match(gene1, phylo$genes), match(gene2, phylo$genes))
    featureMat = (phylo[sl.phylo[,1],-(1:3)] - phylo[sl.phylo[,2],-(1:3)])^2 #Check how different the two genes are
    phylo_coefficient[i] <- sum(featureMat * feature.weight)
    
  } else {
    phylo_coefficient[i] <- 9999999
  }
}
toc()

SL_pairs <- SL_pairs |> mutate(phylo_coefficient = phylo_coefficient)

return(SL_pairs)

# SL_pairs_binarize_expression_path <-
#   file.path(datapath, 'SL_pairs_binarize_expression.csv')
# SL_pairs <- write_csv(SL_pairs, SL_pairs_binarize_expression_path)

}



