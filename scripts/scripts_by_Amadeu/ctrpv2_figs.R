library(scales)
library(ggpubr)
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

SL_pairs <- load_SL_pairs()

SL_pairs_filtered <-SL_pairs |>  group_by(gene1) |>  slice_min(survival_coef, n =20)

results_SL <- validate_ctrpv2(SL_pairs_filtered, filter = FALSE, synthetic_rescue = FALSE)
score_lapa <- results_SL$Lapatinib





 ggscatter(
  data = score_lapa,
  x = "score",
  y = "aac_recomputed",  fill = "#73a9c9", color = "#73a9c9", size = 2, merge = TRUE
) +
  xlab("SL score") + ylab("Lapatinib AAC")+ stat_cor( label.x = 0, alternative = "greater", r.accuracy = 0.001,method = "spearman", size=4.8)+
   font("xlab", size = 14)+
   font("ylab", size = 14)+
   font("x.text", size = 13)+
   font("y.text", size = 13)
 







###################ISLE

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
results_ISLE <- validate_ctrpv2(
  ISLE_drug_filtered,
  synthetic_rescue = FALSE,
  filter = FALSE,
  title = "ISLE drug"
)

score_lapa <- results_ISLE$Lapatinib

  

   ggscatter(
  data = score_lapa,
  x = "score",
  y = "aac_recomputed",  fill = "#73a9c9", color = "#73a9c9", size = 2, merge = TRUE
) +
  xlab("SL score") + ylab("Lapatinib AAC")+ stat_cor( label.x = 0, method="spearman", alternative = "greater", r.accuracy = 0.001, size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)



######################################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox_random2.csv')
SL_pairs_random <- read_csv(SL_pairs_path, show_col_types = FALSE) |> slice_min(q_value, n=1000)
results_SL_random <- validate_ctrpv2(SL_pairs_random, filter = FALSE, synthetic_rescue = FALSE, title="random SL")


score_lapa <- results_SL_random$Lapatinib





ggscatter(
  data = score_lapa,
  x = "score",
  y = "aac_recomputed",  fill = "#73a9c9", color = "#73a9c9", size = 2, merge = TRUE
) +
  xlab("Lapatinib random SL score") + ylab("Lapatinib AAC")+ stat_cor( label.x = 0, method="spearman", alternative = "greater", r.accuracy = 0.001)

 ggscatter(
  data = score_sora,
  x = "score",
  y = "aac_recomputed",  fill = "#73a9c9", color = "#73a9c9", size = 2, merge = TRUE
) +
  xlab("Sorafenib  SL score") + ylab("Sorafenib AAC")+ stat_cor(label.x = 0, method="spearman", alternative = "greater", r.accuracy = 0.001)

score_sora <- results_SL_random$Sorafenib

