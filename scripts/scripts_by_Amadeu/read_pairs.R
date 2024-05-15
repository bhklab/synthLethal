library(scales)
library(ggpubr)
datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

#Binarize essentiality and binarize expression
##################################################


SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network.csv')


SL_pairs_binarize_expression <- read_csv(SL_pairs_path, show_col_types = FALSE)



SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_essentiality_full_network.csv')


SL_pairs_binarize_essentiality <- read_csv(SL_pairs_path, show_col_types = FALSE)




#Make them share the same genes so they can be plotted 1 to 1
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> filter(gene1 %in% SL_pairs_binarize_essentiality$gene1 & (gene2 %in% SL_pairs_binarize_essentiality$gene2))
SL_pairs_binarize_essentiality <- SL_pairs_binarize_essentiality |> filter(gene1 %in% SL_pairs_binarize_expression$gene1 & (gene2 %in% SL_pairs_binarize_expression$gene2))

#order them so the values match
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> arrange(gene1, gene2)
SL_pairs_binarize_essentiality <- SL_pairs_binarize_essentiality|> arrange(gene1, gene2)

plot_data <- tibble(gene1 = SL_pairs_binarize_expression$gene1, gene2 = SL_pairs_binarize_expression$gene2, q_value_binarize_expression = SL_pairs_binarize_expression$q_value+1e-30, q_value_binarize_essentiality = SL_pairs_binarize_essentiality$q_value+1e-25)




ggscatter(data = plot_data |> dplyr::filter(gene1 == "EGFR"),
          x = "q_value_binarize_expression", y = "q_value_binarize_essentiality", size = 1, color = "#73a9c9") + scale_x_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) +  scale_y_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) + 
  xlab("Binarize expression q-value") + ylab("Binarize essentiality q-value")+
  stat_cor(label.x = -29, label.y = 0.8,alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)



    
  






#Pearson and expression
###########################################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network.csv')


SL_pairs_binarize_expression <- read_csv(SL_pairs_path, show_col_types = FALSE)



SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_pearson_full_network.csv')


SL_pairs_pearson <- read_csv(SL_pairs_path, show_col_types = FALSE)


#Make them share the same genes so they can be plotted 1 to 1
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> filter(gene1 %in% SL_pairs_pearson$gene1 & (gene2 %in% SL_pairs_pearson$gene2))
SL_pairs_pearson <- SL_pairs_pearson |> filter(gene1 %in% SL_pairs_binarize_expression$gene1 & (gene2 %in% SL_pairs_binarize_expression$gene2))

#order them so the values match
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> arrange(gene1, gene2)
SL_pairs_pearson <- SL_pairs_pearson|> arrange(gene1, gene2)

plot_data <- tibble(gene1 = SL_pairs_binarize_expression$gene1, gene2 = SL_pairs_binarize_expression$gene2, q_value_binarize_expression = SL_pairs_binarize_expression$q_value+1e-30,  q_value_pearson = SL_pairs_pearson$q_value+1e-30)




ggscatter(data = plot_data |> dplyr::filter(gene1 == "EGFR"),
          x = "q_value_binarize_expression", y = "q_value_pearson", size = 1, color = "#73a9c9") + scale_x_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) +  scale_y_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) + 
  xlab("Binarize expression q-value") + ylab("Pearson q-value")+
  stat_cor(label.x = -29, label.y = 0.8,alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)








plot_data <- tibble(gene1 = SL_pairs_binarize_expression$gene1, gene2 = SL_pairs_binarize_expression$gene2, q_value_binarize_expression = SL_pairs_binarize_expression$q_value, q_value_pearson = SL_pairs_pearson$q_value)


ggplot(
  data = plot_data |> dplyr::filter(gene1 == "EGFR"),
  aes(x = q_value_binarize_expression + 1e-25,
      y = q_value_pearson + 1e-25)
) + scale_x_log10(
  breaks = trans_breaks("log10", function(x)
    10 ^ x),
  labels = trans_format("log10", math_format(10 ^
                                               .x))
) +
  scale_y_log10(
    breaks = trans_breaks("log10", function(x)
      10 ^ x),
    labels = trans_format("log10", math_format(10 ^ .x))
  ) + geom_point(size = 1, color = "#73a9c9")  + 
  
  ggtitle(
    str_c(
      "Plot of q-values of the essentiality-expression Wilcoxon q-value of EGFR via binarizing the expression or Pearson"
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # labs(caption = 
  #      str_c("Synthetic lethal pairs obtained via the method ",
  #      "Binarize expression")) +
  xlab("Binarize expression q-value") + ylab("Pearson q-value")+ stat_cor(aes(label = ..r.label..), label.x = -19, alternative = "greater")




#Binarize ess
######################################################



#Survival coefficients
###########################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox.csv')


SL_pairs_cox <- read_csv(SL_pairs_path, show_col_types = FALSE)

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_parametric.csv')


SL_pairs_parametric <- read_csv(SL_pairs_path, show_col_types = FALSE)

plot_data <- tibble(gene1 = SL_pairs_cox$gene1, gene2 = SL_pairs_cox$gene2, cox = SL_pairs_cox$survival_coef, parametric = SL_pairs_parametric$survival_coef)


  
plot <-   ggscatter(
    data = (plot_data |> dplyr::filter(gene1 == "EGFR")),
  x = "cox",
        y = "parametric",  fill = "#73a9c9", color = "#73a9c9", size = 1, merge = TRUE
  ) +
  xlab("Cox proportional hazard coefficient") + ylab("Parametric exponential coefficient")+ stat_cor(label.x = -0.7, alternative = "greater", r.accuracy = 0.0001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)




ggpar(plot, xlim = c(-0.75, 0.25), ylim = c(-0.75, 0.25))



title = "Comparison of the survival coefficients of a Cox proportional hazard model and a parametric exponential model"

#Inverse
####################################

datapath <-
  "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network.csv')


SL_pairs_binarize_expression <- read_csv(SL_pairs_path, show_col_types = FALSE)



SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network_inverse.csv')


SL_pairs_binarize_expression_inverse <- read_csv(SL_pairs_path, show_col_types = FALSE)





#Make them share the same genes so they can be plotted 1 to 1
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> filter(gene1 %in% SL_pairs_binarize_expression_inverse$gene1 & (gene2 %in% SL_pairs_binarize_expression_inverse$gene2))
SL_pairs_binarize_expression_inverse <- SL_pairs_binarize_expression_inverse |> filter(gene1 %in% SL_pairs_binarize_expression$gene1 & (gene2 %in% SL_pairs_binarize_expression$gene2))

#order them so the values match
SL_pairs_binarize_expression <- SL_pairs_binarize_expression |> arrange(gene1, gene2)
SL_pairs_binarize_expression_inverse <- SL_pairs_binarize_expression_inverse|> arrange(gene1, gene2)

plot_data <- tibble(gene1 = SL_pairs_binarize_expression$gene1, gene2 = SL_pairs_binarize_expression$gene2, q_value_binarize_expression = SL_pairs_binarize_expression$q_value, q_value_binarize_expression_inverse = SL_pairs_binarize_expression_inverse$q_value)




ggscatter(
  data = (plot_data |> dplyr::filter(gene1 == "EGFR")),
  x = "q_value_binarize_expression",
  y = "q_value_binarize_expression_inverse",  fill = "#73a9c9", color = "#73a9c9", size = 1, merge = TRUE
) +
  xlab("Essentiality q-value") + ylab("Expression q-value")+ stat_cor(label.x = 0.3, label.y = 1.05, alternative = "greater", r.accuracy = 0.0001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)







ggplot(
  data = plot_data |> dplyr::filter(gene1 == "EGFR"),
  aes(x = q_value_binarize_expression + 1e-15,
      y = q_value_binarize_expression_inverse + 1e-15)
) + geom_point(size = 1, color = "#73a9c9")  + 
  
  ggtitle(
    str_c(
      "Distribution of the q-values of the essentiality-expression Wilcoxon q-value and of the expression-essentiality Wilcoxon q-value of EGFR"
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("EGFR essentiality q-value") + ylab("EGFR expression q-value")




#Survival vs things
#######################################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox.csv')


SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)

pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(q_value),
    method = "pearson",
    use = "complete.obs"
  )
ggscatter(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = q_value + 1e-20, y = survival_coef)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + 
  ggtitle(
    str_c(
      "Plot of q-values of the Essentiality test and coefficients of the Cox model of EGFR"
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("EGFR essentiality q-value") + ylab("EGFR Cox model coefficient")+
  stat_cor( label.x = -10, label.y = 0.8,alternative = "greater", r.accuracy = 0.001)


#essentiality vs phylogenetic
#############################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_cox.csv')


SL_pairs_binarize_expression <- read_csv(SL_pairs_path, show_col_types = FALSE)

pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(q_value),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = q_value + 1e-20, y = phylo_coefficient)) + geom_point(size = 1, color = "#73a9c9") + 
  ggtitle(
    str_c(
      "Plot of q-values of the Essentiality test and coefficients of the phylogenetic test of EGFR"
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("EGFR essentiality q-value") + ylab("EGFR phylogenetic coefficient")+
  stat_cor( label.x = 0.5, label.y = 170,alternative = "greater", r.accuracy = 0.001)



#essentialiy vs depletion
######################################

ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = q_value + 1e-20, y = depletion_q_value + 1e-15)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + scale_y_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       )+
  ggtitle(
    str_c(
      "Plot of q-values of the Essentiality test and the depletion test of EGFR"
    )
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("EGFR essentiality q-value") + ylab("EGFR depletion q-value")+
  stat_cor(aes(label = ..r.label..), label.x = -10, label.y = 0.8,alternative = "greater", r.accuracy = 0.001)



#Parametric
###################################################

SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_parametric.csv')


SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)

pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(q_value),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = depletion_p_value + 1e-15, y = survival_coef)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + labs(
         title = str_c(
           "Plot of p-values of the CRISPR and coefficients of the Cox model of EGFR, Pearson Correlation = ",
           round(pearson, 3)
         ),
         x =
           "Hypergeometric Depletion p-value",
         y = "Cox model coefficient"
       )



##############################################
SL_pairs_path <-
  file.path(datapath, 'SL_pairs_breast_binarize_expression_cox.csv')


SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE) |> mutate(q_value = q_value+1e-30, depletion_q_value = depletion_q_value+1e-25)


ggscatter(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
          x = "q_value", y = "depletion_q_value", size = 1, color = "#73a9c9") + scale_x_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) +  scale_y_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) + 
  xlab("Essentiality q-value") + ylab("Depletion q-value")+
  stat_cor(label.x = -29, label.y = 0.8,alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)



##############################################

ggscatter(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
          x = "q_value", y = "survival_coef", size = 1, color = "#73a9c9") + scale_x_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) + 
  xlab("Essentiality q-value") + ylab("Cox model coefficient")+
  stat_cor(label.x = -28, label.y = 0.6,alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)



###########################################
ggscatter(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
          x = "q_value", y = "phylo_coefficient", size = 1, color = "#73a9c9") + scale_x_log10(
            breaks = trans_breaks("log10", function(x)
              10 ^ x),
            labels = trans_format("log10", math_format(10 ^
                                                         .x))
          ) + 
  xlab("Essentiality q-value") + ylab("Phylogenetic coefficient")+
  stat_cor(label.x = -28, label.y = 170,alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)

###############Others

ggscatter(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
          x = "phylo_coefficient", y = "q_value", size = 1, color = "#73a9c9") + 
  xlab("Essentiality q-value") + ylab("Depletion q-value")+
  stat_cor( r.accuracy = 0.001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)





