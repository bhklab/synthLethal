datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

####################################Essentiality Test


essentiality_path <-
  file.path(datapath, 'CRISPR_gene_effect.csv')
expression_path <- file.path(datapath, 'CCLE_expression.csv')

model_path <-
  file.path(datapath, 'Model.csv')
model <- read_csv(model_path, show_col_types = FALSE)

gene_effect <-
  read_csv(essentiality_path, show_col_types = FALSE)
expression <-
  read_csv(expression_path, show_col_types = FALSE)


gene_effect_randomised <- gene_effect 



for (i in 2:ncol(gene_effect)){
  
  gene_effect_randomised[i] <- sample(unlist(gene_effect[i]))
  
  
}

random_essentiality_path <-
  file.path(datapath, 'CRISPR_gene_effect_random.csv')

write_csv(gene_effect_randomised, random_essentiality_path)



expression_randomised <- expression 



for (i in 2:ncol(expression)){
  
  expression_randomised[i] <- sample(unlist(expression[i]))
  
  
}

random_expression_path <-
  file.path(datapath, 'CCLE_expression_random.csv')

write_csv(expression_randomised, random_expression_path)

##########################################################Survival

TCGA_t_expression_percentile_path <-
  file.path(datapath, 'TCGA_t_expression_percentile.csv')
t_expression_percentile <- read_csv(TCGA_t_expression_percentile_path, show_col_types = FALSE)



t_expression_percentile_randomised <- t_expression_percentile 



for (i in 2:ncol(t_expression_percentile)){
  
  t_expression_percentile_randomised[i] <- sample(unlist(t_expression_percentile[i]))
  
  
}


TCGA_t_expression_percentile_random_path <-
  file.path(datapath, 'TCGA_t_expression_percentile_random.csv')

write_csv(t_expression_percentile_randomised, TCGA_t_expression_percentile_random_path)



TCGA_survival_path <-
  file.path(datapath, 'TCGA_survival.csv')
survival <- read_csv(TCGA_survival_path, show_col_types = FALSE)


survival_randomised <- survival



for (i in 2:ncol(survival)){
  
 survival_randomised[i] <- sample(unlist(survival[i]))
  
  
}


survival_random_path <-
  file.path(datapath, 'TCGA_survival_random.csv')

write_csv(survival_randomised, survival_random_path)



###############################Phylogenetic


phylo_path <- file.path(datapath, 'yuval.phylogenetic.profile.RData')
load(phylo_path)

phylo_randomised <- phylo

for (i in 3:ncol(phylo)){
  
  phylo_randomised[i] <- sample(unlist(phylo[i]))
  
  
}

phylo_random_path <-
  file.path(datapath, 'yuval.phylogenetic.profile.random.RData')

phylo <- phylo_randomised

save(phylo, file = phylo_random_path)




