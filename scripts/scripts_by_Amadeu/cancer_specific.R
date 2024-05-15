library(tidyverse)
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

model_path <-
  file.path(datapath, 'Model.csv')
model <- read_csv(model_path)

CCLE_cancer_types <- tolower(unique(model$OncotreeLineage))
CCLE_cancer_types

TCGA_survival_path <-
  file.path(datapath, 'TCGA_survival.csv')
survival <- read_csv(TCGA_survival_path)

TCGA_cancer_types <- unique(survival$tumor_tissue_site)
TCGA_cancer_types

drug_targets_path <-
  file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path)

drug_targets_cancer_types <- tolower(unique(drug_targets$tissue))
drug_targets_cancer_types

common_types <- intersect(TCGA_cancer_types, intersect(CCLE_cancer_types, drug_targets_cancer_types))
common_types