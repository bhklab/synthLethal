library(progressr)


cancer_type = "pan-cancer"

crispr_method = "binarize_expression"


datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

print("Undergoing the CRISPR and Depletion test:")
SL_pairs <- with_progress(obtain_pairs_inverse(method = crispr_method, n = n, cancer_type = cancer_type))


  SL_pairs_path <-
    file.path(datapath, str_c("SL_pairs_", cancer_type, "_", crispr_method, "_full_network_inverse", ".csv"))
  
  write_csv(SL_pairs, SL_pairs_path)
