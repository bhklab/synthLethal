library(tidyverse)
library(scales)


#Expression vs essentiality


x <- gene_effect |> select("EGFR") |> pull()
y <- expression |> select("PRR3") |> pull()


plot_data <- tibble(Essentiality = x, Expression = y)

# ggplot(data = plot_tibble, aes(x=Essentiality, y=Expression)) + geom_point(color = "#73a9c9") + labs(title="Expression of PRR3 vs Essentiality of EGFR",
#                                                                                    x ="Essentiality of EGFR", y = "Expression of PRR3")


p1 <- ggscatter(
  data = plot_data,
  x = "Essentiality",
  y = "Expression", color = "#73a9c9", size = 1
) +
  
  
  xlab("Essentiality of EGFR") + ylab("Expression of PRR3") + stat_cor(alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)





#Essentiality vs expression

x <- gene_effect |> select("PRR3") |> pull()
y <- expression |> select("EGFR") |> pull()


plot_data <- tibble(Essentiality = x, Expression = y)

# ggplot(data = plot_tibble, aes(x=Essentiality, y=Expression)) + geom_point(color = "#73a9c9") + labs(title="Expression of PRR3 vs Essentiality of EGFR",
#                                                                                    x ="Essentiality of EGFR", y = "Expression of PRR3")


p2 <- ggscatter(
  data = plot_data,
  x = "Essentiality",
  y = "Expression", color = "#73a9c9", size = 1
) +
  
  
  xlab("Essentiality of PRR3") + ylab("Expression of EGFR") + stat_cor(alternative = "greater", r.accuracy = 0.001, method = "spearman", size=4.8)+  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("x.text", size = 13)+
  font("y.text", size = 13)





ggarrange(p1,p2, widths = c(1,1))





########comparison of p-values


ggplot(
  data = SL_pairs |> dplyr::filter(gene1 == "ERBB2"),
  aes(
    x = q_value,
    y = depletion_q_value,
    colour = "cadetblue"
  )
) + geom_point() + labs(title = "Comparison of p-values of EGFR",
                        x =
                          "Wilcoxon p-value", y = "Hypergeometric Depletion p-value")







pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(p_value),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(depletion_p_value),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(
  data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
  aes(x = p_value + 1e-30,
      y = depletion_p_value + 1e-25)
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
  ) + geom_point(size = 1, color = "#73a9c9")  + labs(
    title = str_c(
      "Plot of p-values of the Wilcoxon and Hypergeometric Depletion Test of EGFR, Pearson Correlation = ",
      round(pearson, 3)
    ),
    x = "Wilcoxon p-value",
    y = "Hypergeometric Depletion p-value"
  )


pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(p_value),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    method = "pearson",
    use = "complete.obs"
  )
ggplot(data = SL_pairs |> dplyr::filter(gene1 == "EGFR"),
       aes(x = p_value + 1e-20, y = survival_coef)) + geom_point(size = 1, color = "#73a9c9") + scale_x_log10(
         breaks = trans_breaks("log10", function(x)
           10 ^ x),
         labels = trans_format("log10", math_format(10 ^
                                                      .x))
       ) + labs(
         title = str_c(
           "Plot of p-values of the Wilcoxon Test and coefficients of the Cox model of EGFR, Pearson Correlation = ",
           round(pearson, 3)
         ),
         x =
           "Wilcoxon p-value",
         y = "Cox model coefficient"
       )

pearson <-
  cor(
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(survival_coef),
    SL_pairs |> dplyr::filter(gene1 == "EGFR") |> select(depletion_p_value),
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
           "Plot of p-values of the Hypergeometric Depletion Test and coefficients of the Cox model of EGFR, Pearson Correlation = ",
           round(pearson, 3)
         ),
         x =
           "Hypergeometric Depletion p-value",
         y = "Cox model coefficient"
       )

####Histogram of p-value and q-value

ggplot(data = SL_pairs, aes(x=p_value, y="density")) + geom_histogram(binwidth = 0.005, colour =  "black", fill =  "#73a9c9") + labs(title="Histogram of Wilcoxon p-values",
                                                                                                                       x ="Wilcoxon p-value", y = "Frequency")

ggplot(data = SL_pairs, aes(x=q_value)) + geom_histogram(binwidth = 0.005, colour =  "black", fill =  "#73a9c9") + labs(title="Histogram of Wilcoxon q-values",
                                                                                                                        x ="Wilcoxon q-value", y = "Frequency")



SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network_random.csv')


SL_pairs_random <- read_csv(SL_pairs_path, show_col_types = FALSE)

p1 <- gghistogram(data = SL_pairs_random, x = "p_value", colour =  "#73a9c9", fill =  "#73a9c9", bins=50) + xlab("Essentiality test p-value") 

p2 <- gghistogram(data = SL_pairs_random, x = "depletion_p_value", colour =  "#73a9c9", fill =  "#73a9c9", bins=50)+ xlab("Depletion test p-value") 

ggarrange(p1,p2)


SL_pairs_path <-
  file.path(datapath, 'SL_pairs_pan-cancer_binarize_expression_full_network.csv')


SL_pairs <- read_csv(SL_pairs_path, show_col_types = FALSE)

p1 <- gghistogram(data = SL_pairs, x = "p_value", fill =  "#73a9c9", colour =  "#73a9c9",  bins=50) + xlab("Essentiality test p-value") 

p2 <- gghistogram(data = SL_pairs, x = "depletion_p_value", colour =  "#73a9c9", fill =  "#73a9c9", bins=50)+ xlab("Depletion test p-value") 

ggarrange(p1,p2)








