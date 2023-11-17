#' Create a synthetic lethal network
#' 
#' createEssentialNetwork takes as arguments a path to Cancer Dependency Map essentiality data
#' and parameters defining the network properties. It then computes a gene essentiality 
#' network. All the variations produce a test statistic and a p-value.
#' 
#' @param datapath  Path to directory containing Dependency Map data
#' @param outpath   Output path to save .rds object, optionally including filename.
#' @param method    Character string describing which method to use. One of c("ISLE")
#' @param n         The number of genes to compute associations for, useful for debugging and testing.
#' Default is 0, which computes the entire network. 
#' @param tissue    Character string defining which tissue to use. Default is "all"
#' 
#' @returns         A dataframe containing a matrix of effect sizes (effectSizes), a matrix of p-values
#' (pvals), and a list of network parameters. 
#' @export
#' @importFrom stats wilcox.test
#' 

createEssentialNetwork <- function(datapath=".", outpath=".", method="ISLE", n=0, 
                                   tissue="all"){
  
  
  
  
  
}