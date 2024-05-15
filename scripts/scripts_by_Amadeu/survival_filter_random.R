
#' Title
#'
#' @param SL_pairs 
#' @param cancer_type 
#' @param method 
#'
#' @return
#' 
#'
#' @examples
survival_filter_random <- function(SL_pairs, cancer_type = "pan-cancer", method = "cox"){



library(survival)
library(tidyverse)
library(tictoc)
library(ggsurvfit)

datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"



TCGA_t_expression_percentile_path <-
  file.path(datapath, 'TCGA_t_expression_percentile_random.csv')
t_expression_percentile <- read_csv(TCGA_t_expression_percentile_path, show_col_types = FALSE)

TCGA_survival_path <-
  file.path(datapath, 'TCGA_survival_random.csv')
survival <- read_csv(TCGA_survival_path, show_col_types = FALSE)


if (cancer_type != "pan-cancer"){
  
  survival <- survival |> dplyr::filter(tumor_tissue_site == cancer_type)
  t_expression_percentile <- t_expression_percentile |> dplyr::filter(patientID %in% survival$patientID)
}


coefficients <- vector(mode="numeric", length = nrow(SL_pairs))
survival_p_value <- vector(mode="numeric", length = nrow(SL_pairs))


genes <- colnames(t_expression_percentile)[-1]

p <- progressr::progressor(nrow(SL_pairs) %/% 1000)

#SL_pairs <- SL_pairs |> filter(p_value < 0.01, depletion_p_value < 0.01)

coefficients <- vector(mode="numeric", length = nrow(SL_pairs))
survival_p_value <- vector(mode="numeric", length = nrow(SL_pairs))

if(method == "cox"){

########################################################  
#Binary Cox test  

for (i in  1:nrow(SL_pairs)){
  if(i %% 1000 == 0){
    p()
  }
  gene1 <- SL_pairs$gene1[i]
  gene2 <-  SL_pairs$gene2[i]
  if (gene1 %in% genes && gene2 %in% genes){
    
    
    survival <- survival |> mutate(SL_index = if_else(select(t_expression_percentile,{{gene1}}) < 0.33 & select(t_expression_percentile,{{gene2}}) < 0.33, 1, 0)) #Only one "&" because it is vectorised. If the SL pair is underexpressed, so it should have an effect, you get a 1
    if(cancer_type == "pan-cancer"){
    cox_fit <-  coxph(Surv(days, vital_status) ~ gender + years_to_birth + tumor_tissue_site + pathology_T_stage + pathology_N_stage + pathology_M_stage + SL_index , data = survival)
    } else {
      cox_fit <-  coxph(Surv(days, vital_status) ~ gender + years_to_birth + pathology_T_stage + pathology_N_stage + pathology_M_stage + SL_index , data = survival)
    }
    coefficients[i] <- cox_fit$coefficients["SL_index"]
    survival_p_value[i] <- summary(cox_fit)[['coefficients']][,'Pr(>|z|)']["SL_index"]
    
    
  }
}
}

else if (method == "parametric"){
  for (i in  1:nrow(SL_pairs)){
    if(i %% 1000 == 0){
      p()
    }
    gene1 <- SL_pairs$gene1[i]
    gene2 <-  SL_pairs$gene2[i]
    if (gene1 %in% genes && gene2 %in% genes){
      
      
      survival <- survival |> mutate(SL_index = if_else(select(t_expression_percentile,{{gene1}}) < 0.33 & select(t_expression_percentile,{{gene2}}) < 0.33, 1, 0)) #Only one "&" because it is vectorised. If the SL pair is underexpressed, so it should have an effect, you get a 1
      if(cancer_type == "pan-cancer"){
        par_fit <- survreg(Surv(days, vital_status) ~ gender + years_to_birth + tumor_tissue_site + pathology_T_stage + pathology_N_stage + pathology_M_stage + SL_index, data = survival, dist = 'exponential')
      } else {
        par_fit <- survreg(Surv(days, vital_status) ~ gender + years_to_birth + pathology_T_stage + pathology_N_stage + pathology_M_stage + SL_index, data = survival, dist = 'exponential')
      }
      coefficients[i] <- par_fit$coefficients["SL_index"]*(-1) #the coefficients have opposite sign in the parametric case, but given we only use the magnitude and sign of the coefficients, we flip them.
      survival_p_value[i] <- summary(par_fit)$table["SL_index","p"]
      
      
    }
  }
}
  
  
  
  
  

SL_pairs <- SL_pairs |> mutate(survival_coef = coefficients, survival_p_value = survival_p_value) |> filter(survival_coef != '') |> mutate(survival_p_value = if_else(survival_p_value == 0, 1, survival_p_value)) #We want to keep those t which there is a TCGA expression, those without an effect whould have a p_value of 1


# var_survival_q_value <-
#   p.adjust(SL_pairs$survival_p_value, method = "fdr")
# 
# SL_pairs <- SL_pairs |> mutate(survival_q_value = var_survival_q_value)

# SL_pairs_survival_path <-
#   file.path(datapath, 'SL_pairs_survival.csv')
# SL_pairs <- write_csv(SL_pairs, SL_pairs_survival_path)

return(SL_pairs)

}













