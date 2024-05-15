library(ggpubr)
library(pROC)
library(tidyverse)

############################SLKB
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
SL_pairs_path <-
  file.path(datapath, "SLKB_pairs.csv")
SLKB_pairs <- readr::read_csv(SL_pairs_path, show_col_types = FALSE)
drug_targets_path <- file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path) |>
  mutate(genes = str_split(genes, ","))
targeted_genes <- character()
for (i in 1:nrow(drug_targets)) {
  genes <- drug_targets$genes[i][[1]]
  for (gene in genes) {
    targeted_genes <- append(targeted_genes, gene)
  }
}
targeted_genes <- unique(targeted_genes)
SLKB_pairs_filtered <- SLKB_pairs
for (i in 1:nrow(SLKB_pairs)) {
  if (SLKB_pairs$gene2[i] %in% targeted_genes) {
    if (SLKB_pairs$gene1[i] %in% targeted_genes) {
      SLKB_pairs_filtered <- bind_rows(SLKB_pairs_filtered, slice(SLKB_pairs, i))
    }
    ph <- SLKB_pairs_filtered$gene2[i]
    SLKB_pairs_filtered$gene2[i] <- SLKB_pairs_filtered$gene1[i]
    SLKB_pairs_filtered$gene1[i] <- ph
  }
}
SLKB_pairs_filtered <- SLKB_pairs_filtered |> distinct(gene1, gene2, .keep_all = TRUE)
SLKB_pairs_filtered <- SLKB_pairs_filtered |> filter(gene1 %in% targeted_genes) |> arrange(gene1, gene2)
results_SLKB <- validate_enlight(
  SLKB_pairs,
  synthetic_rescue = FALSE,
  filter = FALSE,
  title = "SLKB"
)
###############################SL pairs
SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox.csv')
SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)
results_SL <- validate_enlight(SL_pairs)
print("intersect between me and SLKB:")
length(intersect(
  str_c(SLKB_pairs_filtered$gene1, SLKB_pairs_filtered$gene2),
  str_c(SL_pairs$gene1, SL_pairs$gene2)
))
#52/377 i funciona millor
SL_pairs_filtered <-  SL_pairs |> dplyr::filter(q_value < 0.1,
                                                depletion_q_value < 0.1,
                                                phylo_coefficient < 10.5,
                                                phylo_coefficient > 0)   #USe the same boundary for the phylogenetic coefficient as ISLE

results_SL <- validate_enlight(SL_pairs_filtered)
print("If I reduce it to n = 87, to match ISLE")
length(intersect(
  str_c(SLKB_pairs_filtered$gene1, SLKB_pairs_filtered$gene2),
  str_c(SL_pairs_filtered$gene1, SL_pairs_filtered$gene2)
))
###########################################ISLE pairs
SL_pairs_path <-
  file.path(datapath, 'ISLE_0.2.csv')
ISLE_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)
ISLE_pairs <- ISLE_pairs |> mutate(gene1 = word(ISLE_pairs$name),
                                   gene2 = word(ISLE_pairs$name, 4))
ISLE_pairs_filtered <- ISLE_pairs
for (i in 1:nrow(ISLE_pairs)) {
  if (ISLE_pairs$gene2[i] %in% targeted_genes) {
    if (ISLE_pairs$gene1[i] %in% targeted_genes) {
      ISLE_pairs_filtered <- bind_rows(ISLE_pairs_filtered, slice(ISLE_pairs, i))
    }
    ph <- ISLE_pairs_filtered$gene2[i]
    ISLE_pairs_filtered$gene2[i] <- ISLE_pairs_filtered$gene1[i]
    ISLE_pairs_filtered$gene1[i] <- ph
  }
}
ISLE_pairs_filtered <- ISLE_pairs_filtered |> distinct(gene1, gene2, .keep_all = TRUE)
ISLE_pairs_filtered <- ISLE_pairs_filtered |> filter(gene1 %in% targeted_genes) |> arrange(gene1, gene2)
results_ISLE <- validate_enlight(
  ISLE_pairs_filtered,
  synthetic_rescue = FALSE,
  filter = FALSE,
  title = "ISLE pairs"
)
length(intersect(
  str_c(SLKB_pairs_filtered$gene1, SLKB_pairs_filtered$gene2),
  str_c(ISLE_pairs_filtered$gene1, ISLE_pairs_filtered$gene2)
))
length(intersect(
  str_c(SL_pairs$gene1, SL_pairs$gene2),
  str_c(ISLE_pairs_filtered$gene1, ISLE_pairs_filtered$gene2)
))
length(intersect(
  str_c(SL_pairs_filtered$gene1, SL_pairs_filtered$gene2),
  str_c(ISLE_pairs_filtered$gene1, ISLE_pairs_filtered$gene2)
))
####################
SL_pairs_path <-
  file.path(datapath, 'ISLE_drug.csv')
ISLE_drug <- read_csv(SL_pairs_path, show_col_types = FALSE)
ISLE_drug <- ISLE_drug |> mutate(gene1 = word(ISLE_drug$name),
                                 gene2 = word(ISLE_drug$name, 4))
ISLE_drug_filtered <- ISLE_drug
for (i in 1:nrow(ISLE_drug)) {
  if (ISLE_drug$gene2[i] %in% targeted_genes) {
    if (ISLE_drug$gene1[i] %in% targeted_genes) {
      ISLE_drug_filtered <- bind_rows(ISLE_drug_filtered, slice(ISLE_drug, i))
    }
    ph <- ISLE_drug_filtered$gene2[i]
    ISLE_drug_filtered$gene2[i] <- ISLE_drug_filtered$gene1[i]
    ISLE_drug_filtered$gene1[i] <- ph
  }
}
ISLE_drug_filtered <- ISLE_drug_filtered |> distinct(gene1, gene2, .keep_all = TRUE)
ISLE_drug_filtered <- ISLE_drug_filtered |> filter(gene1 %in% targeted_genes) |> arrange(gene1, gene2)
results_ISLE <- validate_enlight(
  ISLE_drug_filtered,
  synthetic_rescue = FALSE,
  filter = FALSE,
  title = "ISLE drug"
)
length(intersect(
  str_c(SLKB_pairs_filtered$gene1, SLKB_pairs_filtered$gene2),
  str_c(ISLE_drug_filtered$gene1, ISLE_drug_filtered$gene2)
))
length(intersect(
  str_c(SL_pairs$gene1, SL_pairs$gene2),
  str_c(ISLE_drug_filtered$gene1, ISLE_drug_filtered$gene2)
))
length(intersect(
  str_c(SL_pairs_filtered$gene1, SL_pairs_filtered$gene2),
  str_c(ISLE_drug_filtered$gene1, ISLE_drug_filtered$gene2)
))
length(intersect(
  str_c(SLKB_pairs_filtered$gene1, SLKB_pairs_filtered$gene2),
  str_c(ISLE_drug_filtered$gene1, ISLE_drug_filtered$gene2)
))




#######################Random

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox_random.csv')
SL_pairs_random <- read_csv(SL_pairs_path, show_col_types = FALSE) |> slice_min(q_value, n=1000, with_ties = FALSE)
results_SL_random <- validate_enlight(SL_pairs_random, filter = FALSE, synthetic_rescue = FALSE, title="random SL")

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox_random2.csv')
SL_pairs_random <- read_csv(SL_pairs_path, show_col_types = FALSE) |> slice_min(q_value, n=1000, with_ties = FALSE)
results_SL_random <- validate_enlight(SL_pairs_random, filter = FALSE, synthetic_rescue = FALSE, title="random SL")




##########################CREATE PLOTSSSSS

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_breast_binarize_expression_cox.csv')
SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)
results_SL <- validate_enlight(SL_pairs, tissue = "breast")





