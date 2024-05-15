library(tidyverse)


datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"


CTRP_exps_path <- file.path(datapath, "CTRP-exps.csv")

CTRP_exps <- read_csv(CTRP_exps_path)

CTRP_drugs_path <- file.path(datapath, "CTRP-drugs.csv")

CTRP_drugs <- read_csv(CTRP_drugs_path)


CTRP_cells_path <- file.path(datapath, "CTRP-cells.csv")

CTRP_cells <- read_csv(CTRP_cells_path)

CTRP_exps <- CTRP_exps |> dplyr::select(treatmentid, sampleid, aac_recomputed, ic50_recomputed)
treatment_ids <- unique(CTRP_exps$treatmentid)
sample_ids <- unique(CTRP_exps$sampleid)

CTRP_exps <- CTRP_exps |> dplyr::filter(!is.na(ic50_recomputed), !is.na(aac_recomputed)) 




drug_targets_path <- file.path(datapath, 'drug_targets.csv')
drug_targets <- read_csv(drug_targets_path) |>
  mutate(genes = str_split(genes, ","))


targeted_drugs <- unlist(drug_targets$drug)
targeted_drugs <- tolower(unique(word(targeted_drugs, sep = fixed("_"))))

model_path <-
  file.path(datapath, 'Model.csv')
model <- read_csv(model_path, show_col_types = FALSE)

CTRP_exps <- CTRP_exps |> filter(tolower(treatmentid) %in% targeted_drugs)

model <- model |> rename(sampleid = CellLineName) |> select(ModelID, sampleid)



CTRP_exps <- inner_join(CTRP_exps, model, by = "sampleid")


cell_lines_1 <- CTRP_exps$ModelID
cell_lines_2 <- expression[[1]]
cell_lines <- intersect(cell_lines_1, cell_lines_2)

CTRP_exps <- CTRP_exps |> filter(ModelID %in% expression[[1]])


CTRP_exps_path <- file.path(datapath, "CTRP-exps.csv")
write_csv(CTRP_exps, CTRP_exps_path)


















expression_path <- file.path(datapath, 'CCLE_expression.csv')
expression <-
  read_csv(expression_path, show_col_types = FALSE)

expression <- expression[expression[[1]] %in% unlist(CTRP_exps$ModelID), ]
expression <- expression |> dplyr::rename(DepMap_ID = names(expression)[1])
expression <- expression[order(expression$DepMap_ID), ]
colnames(expression) <- word(colnames(expression))
expression <- as_tibble(expression)








t_expression <- expression |>
  pivot_longer(cols = -DepMap_ID, names_to = 'index') |>
  pivot_wider(names_from = DepMap_ID, values_from = value)

t_expression <- as_tibble(t_expression)


t_expression_path <-
  file.path(datapath, 't_expression.csv')
write_csv(t_expression, t_expression_path)
