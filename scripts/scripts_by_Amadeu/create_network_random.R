#' Create a synthetic lethal network
#'
#'
#' Create A tibble of gene pairs with their synthetic lethality coefficients according to tests of CRISPR knockout, depletion, survival, and phylogeny.
#' @param cancer_type string establishing the type of tissue of the tumor of interest: `"pan-cancer` (the default), `"breast"`, `"lung"`, and others.
#' @param crispr_method string. The method to use when establishing a relationship between a gene pair's essentiality and expression: `"binarize_expression"`, `"binarize_essentiality"`, or `"pearson"`
#' @param survival_method string. The model to use for the survival data:`"cox"` or `"parametric"` (Exponential)
#' @param n 
#' @param synthetic_rescue 
#' @param save boolean. Whether to save the resulting network in memory
#' @param filter boolean. Whether to filter the unlikely gene pairs to increase speed.
#'
#' @return A tibble of gene pairs in the first two columns and the tests coefficients of those pairs informing how likely they are to have a synthetic lethal relationship
#' @export
#'
#' @import dplyr
#' 
#' @examples 
#' create_network()
#' create_network(cancer_type = "lung", crispr_method = "pearson", survival_method = "parametric")
create_network_random <- ##Create a network of SL and SR candidates
  function(cancer_type = "pan-cancer",
           crispr_method = "binarize_expression", survival_method = "cox", save = TRUE, filter = TRUE) {
    
    library(progressr)
    
    datapath <-
      "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
    
    print("Undergoing the CRISPR and Depletion test:")
    SL_pairs <- with_progress(obtain_pairs_random(method = crispr_method, cancer_type = cancer_type))
    
    if(save){
      SL_pairs_path <-
        file.path(datapath, str_c("SL_pairs_", cancer_type, "_", crispr_method, "_full_network_random", ".csv"))
      
      write_csv(SL_pairs, SL_pairs_path)
    }
    
    if (filter == TRUE){
    SL_pairs <-  SL_pairs |> slice_min(q_value, n = 5000, with_ties = FALSE) #Filter the pairs so it is faster
    }
    
    
    
    print("Undergoing the Survival test:")
    SL_pairs <- progressr::with_progress(survival_filter_random(SL_pairs = SL_pairs, cancer_type = cancer_type, method = survival_method))
    
    print("Undergoing the Phylogenetic test:")
    SL_pairs <- progressr::with_progress(phylogenetic_test_random(SL_pairs))
    
    
    if(save){
    SL_pairs_path <-
      file.path(datapath, str_c("SL_pairs_", cancer_type, "_", crispr_method, "_", survival_method, "_random2.csv"))
    
    write_csv(SL_pairs, SL_pairs_path)
    }
    
    return(SL_pairs)
    
  }