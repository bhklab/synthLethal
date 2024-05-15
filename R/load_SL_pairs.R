#' Title
#'
#' @return The SL_pairs tibble of the chosen characteristic
#' @export
#' 
#' @examples SL_pairs <- load_SL_pairs()
load_SL_pairs <- function(cancer_type = "pan-cancer",
                          crispr_method = "binarize_expression", survival_method = "cox"){
  
  datapath <- file.path(getwd(), "data")
  
  SL_pairs_path <-
    file.path(datapath, stringr::str_c("SL_pairs_", cancer_type, "_", crispr_method, "_", survival_method, ".csv"))
  

  
  SL_pairs <- readr::read_csv(SL_pairs_path, show_col_types = FALSE)
  
  SL_pairs
  
  
  
}