#' Obtain synthetic lethal pairs using the DepMap CRISPR datasets
#'
#' @param method A string. How should the pairs be computed: : `"binarize_expression"`, `'binarize_essentiality'`, `"pearson"`
#'
#' @param n A number. How many genes to use for the computation
#' @param max_p_value A number between 0 and 1. Below which p_value should the pairs be saved.
#'
#' @return A tibble with each pairs and its corresponding statistic
#' 
#'
#' @examples SL_pairs <- obtain_pairs()
obtain_pairs <-
  function(method = "binarize_expression",
           parallel = FALSE,
           cancer_type = "pan-cancer") {
    library(tidyverse)
    datapath <-
      "C:/Users/aikon/OneDrive/Desktop/TFG/SL/synthLethal/data"
    
    essentiality_path <-
      file.path(datapath, 'CRISPR_gene_effect.csv')
    expression_path <- file.path(datapath, 'CCLE_expression.csv')
    
    model_path <-
      file.path(datapath, 'Model.csv')
    model <- read_csv(model_path, show_col_types = FALSE)
    
    
    if (TRUE) {
      #!exists("gene_effect") or doesnt have the size it has to have
      gene_effect <-
        read_csv(essentiality_path, show_col_types = FALSE)
      expression <-
        read_csv(expression_path, show_col_types = FALSE)
      
      cell_lines_1 <- gene_effect$DepMap_ID
      cell_lines_2 <- expression[[1]]
      cell_lines <- intersect(cell_lines_1, cell_lines_2)
      
      
      
      gene_effect <-
        gene_effect[gene_effect$DepMap_ID %in% cell_lines, ]
      colnames(gene_effect) <- word(colnames(gene_effect))
      gene_effect <- as_tibble(gene_effect)
      gene_effect <- gene_effect[order(gene_effect$DepMap_ID), ]
      expression <- expression[expression[[1]] %in% cell_lines, ]
      expression <- expression |> dplyr::rename(DepMap_ID = names(expression)[1])
      expression <- expression[order(expression$DepMap_ID), ]
      colnames(expression) <- word(colnames(expression))
      expression <- as_tibble(expression)
      
    }
    
    if (cancer_type != "pan-cancer") {
      type_patients <-
        model |> mutate(OncotreeLineage = tolower(OncotreeLineage)) |>  dplyr::filter(OncotreeLineage == cancer_type) |>  select(ModelID) |> unlist()
      expression <-
        expression |> dplyr::filter(DepMap_ID %in% type_patients)
      gene_effect <-
        gene_effect |> dplyr::filter(DepMap_ID %in% type_patients)
    }
    
    expression_genes <-
      colnames(expression)[-1] #Remove the DepMapID
    effect_genes <- colnames(gene_effect)[-1]
    
    n_expression_genes <- ncol(expression)
    
    drug_targets_path <- file.path(datapath, 'drug_targets.csv')
    targeted_genes <-
      read_csv(drug_targets_path, show_col_types = FALSE) |>
      transmute(genes = str_split(genes, ",")) |>
      unlist() |>
      unique()
    
    p <- progressr::progressor(length(targeted_genes) + 3)
    p()
    
    
    total_cell_lines <- nrow(gene_effect)
    cat(total_cell_lines, '\n')
    benchmark <- total_cell_lines %/% 3
    
    # if (n == "all") {
    #   total_gene_As <- ncol(gene_effect)
    #   total_gene_Bs <- ncol(expression)
    #
    # } else {
    #   total_gene_As <- n + 1
    #   total_gene_Bs <- n + 1
    #
    # }
    
    
    SL_pairs <-
      tibble(
        gene1 = character(),
        gene2 = character(),
        p_value = numeric(),
        statistic = numeric(),
        
      )
    
    
    
    
    
    
    
    if (cancer_type == "pan-cancer") {
      thirds_path <-
        file.path(datapath, 'thirds.RData')
      
      load(thirds_path) #Get the bottom_thirds and top_thirds
      
    }
    
    if (!exists("bottom_thirds") || !exists("top_thirds")) {
      top_thirds <- lst()
      bottom_thirds <- lst()
      
      for (gene in expression_genes) {
        #Guardaho que mai canviara
        
        y <- expression |> select(DepMap_ID, {
          {
            gene
          }
        })
        
        bottom_third <-
          y |>  slice_min(order_by = .data[[gene]], prop = 0.33)  #It will get more than the benchmark if the minimum value is the same for more than benchmark cell lines (0) Better than the alternative because otherwise it would be more arbitrary
        top_third <-
          y |>  slice_max(order_by = .data[[gene]], prop = 0.33) |> filter(.data[[gene]] > 0)  #given that sometimes anything above 0 is selected, bimodality could be implemented
        
        
        
        top_thirds[[gene]] <- top_third
        bottom_thirds[[gene]] <- bottom_third
        
      }
    }
    
    p()
    
    if (!parallel) {
      if (method == "binarize_expression") {
        for (gene_A in targeted_genes) {
          p()
          
          
          
          if (gene_A %in% effect_genes) {
            x <- gene_effect |> select(DepMap_ID, {
              {
                gene_A
              }
            })
      
            #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
            for (gene_B in expression_genes) {
              bottom_third <- bottom_thirds[[gene_B]]
              top_third <- top_thirds[[gene_B]]
              
              
              essentiality_bottom <-
                x |>  dplyr::filter(DepMap_ID %in% bottom_third$DepMap_ID) |> select({
                  {
                    gene_A
                  }
                }) |>  unlist() |> as.numeric()
              essentiality_top <-
                x |>  dplyr::filter(DepMap_ID %in% top_third$DepMap_ID) |> select({
                  {
                    gene_A
                  }
                }) |> unlist() |> as.numeric()
              
              #Do the wilcoxon test
              if (length(top_third$DepMap_ID) > 0) {
                wilcoxon_sl <-
                  wilcox.test(essentiality_top,
                              essentiality_bottom,
                              alternative = "greater") #Check whether essentiality_top is greater than essentiality_bottom (when gene y is missing the cell is more likely to die when gene x is removed, low essentiality means the cell dies)
                
              }
              else {
                wilcoxon_sl <- list("p.value" = NA)
              }
              
              
              if (!is.nan(wilcoxon_sl$p.value)) {
                SL_pairs <- SL_pairs |>
                  add_row(
                    gene1 = word(gene_A),
                    gene2 = word(gene_B),
                    p_value = wilcoxon_sl$p.value,
                    statistic = wilcoxon_sl$statistic
                  )
              }
              
              
            }
          }
        }
      } else if (method == 'binarize_essentiality') {
        if (!exists("bottom_essentialities") ||
            !exists("top_essentialities")) {
          top_essentialities <- lst()
          bottom_essentialities <- lst()
          
          for (gene in effect_genes) {
            #Guardaho que mai canviara
            
            y <- gene_effect |> select(DepMap_ID, {
              {
                gene
              }
            })
            
            bottom_essentiality <-
              y |>  slice_min(order_by = .data[[gene]], prop = 0.33)  #It will get more than the benchmark if the minimum value is the same for more than benchmark cell lines (0) Better than the alternative because otherwise it would be more arbitrary
            top_essentiality <-
              y |>  slice_max(order_by = .data[[gene]], prop = 0.33) |> filter(.data[[gene]] > 0)  #given that sometimes anything above 0 is selected, bimodality could be implemented
            
            
            
            top_essentialities[[gene]] <- top_essentiality
            bottom_essentialities[[gene]] <- bottom_essentiality
            
          }
        }
        
        
        
        for (gene_A in expression_genes) {
         
          
          
          
        
            
            x <- expression |> select(DepMap_ID, {
              {
                gene_A
              }
            })
        
            #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
            for (gene_B in targeted_genes) {
              if(gene_B %in%effect_genes){
              bottom_essentiality <- bottom_essentialities[[gene_B]]
              top_essentiality <- top_essentialities[[gene_B]]
              
              
              expression_bottom <-
                x |>  dplyr::filter(DepMap_ID %in% bottom_essentiality$DepMap_ID) |> select({
                  {
                    gene_A
                  }
                }) |>  unlist() |> as.numeric()
              expression_top <-
                x |>  dplyr::filter(DepMap_ID %in% top_essentiality$DepMap_ID) |> select({
                  {
                    gene_A
                  }
                }) |> unlist() |> as.numeric()
              
              #Do the wilcoxon test
              if (length(top_essentiality$DepMap_ID) > 0) {
                wilcoxon_sl <-
                  wilcox.test(expression_top,
                              expression_bottom,
                              alternative = "greater") #Check whether expression_top is greater than expression_bottom (when gene y is missing the cell is more likely to die when gene x is removed, low essentiality means the cell dies)
                
              }
              else {
                wilcoxon_sl <- list("p.value" = NA)
              }
              
              if (!is.nan(wilcoxon_sl$p.value)) {
                SL_pairs <- SL_pairs |>
                  add_row(
                    gene1 = word(gene_B),
                    gene2 = word(gene_A),
                    p_value = wilcoxon_sl$p.value,
                    statistic = wilcoxon_sl$statistic
                  )
              }
              }
              
            }
          
        }
        
        
        #Continuous
      } else if (method == "pearson") {
        if (!exists("bottom_thirds") || !exists("top_thirds")) {
          top_thirds <- lst()
          bottom_thirds <- lst()
          
          
          for (gene in expression_genes) {
            #Guardaho que mai canviara
            
            y <- expression |> select(DepMap_ID, {
              {
                gene
              }
            })
            
            bottom_third <-
              y |>  slice_min(order_by = .data[[gene]], prop = 0.33)  #It will get more than the benchmark if the minimum value is the same for more than benchmark cell lines (0) Better than the alternative because otherwise it would be more arbitrary
            top_third <-
              y |>  slice_max(order_by = .data[[gene]], prop = 0.33) |> filter(.data[[gene]] > 0)  #given that sometimes anything above 0 is selected, bimodality could be implemented
            
            
            
            top_thirds[[gene]] <- top_third
            bottom_thirds[[gene]] <- bottom_third
            
          }
          
        }
        
        for (gene_A in targeted_genes) {
         p()
          
          
          
          if (gene_A %in% effect_genes) {
            x <- gene_effect |> select(DepMap_ID, {
              {
                gene_A
              }
            })
            
            
            #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
            for (gene_B in expression_genes) {
              y <- expression |> select(DepMap_ID, {
                {
                  gene_B
                }
              })
              
              
              cor <-
                cor.test(
                  as.numeric(unlist(x[2])),
                  as.numeric(unlist(y[2])),
                  method = "pearson",
                  alternative = "greater"
                )
              
              
              
              
              
              if (!is.nan(cor$p.value)) {
                SL_pairs <- SL_pairs |>
                  add_row(
                    gene1 = word(gene_A),
                    gene2 = word(gene_B),
                    p_value = cor$p.value,
                    statistic = cor$estimate,
                    
                  )
              }
            }
          }
        }
      }
    }
    
    if (parallel) {
      if (method == "binarize_expression") {
        SL_pairs <-
          foreach (
            i = 1:length(targeted_genes),
            .packages = c("tidyverse", "stats"),
            .combine = rbind
          ) %dopar% {
            p()
            
            SL_pairs_i <-
              tibble(
                gene1 = character(),
                gene2 = character(),
                p_value = numeric(),
                depletion_p_value = numeric()
              )
            
            gene_A <- targeted_genes[i]
            
            if (gene_A %in% effect_genes) {
              x <- gene_effect |> select(DepMap_ID, {
                {
                  gene_A
                }
              })
              
              
              if (gene_A %in% colnames(expression)) {
                bottom_third_A <- bottom_thirds[[gene_A]]
              }
              #Only check gene A because it is the gene that has been "removed", which is what a drug would do too (VEGFA included because it is not in gene_effect for some reason)
              for (gene_B in expression_genes) {
                bottom_third <- bottom_thirds[[gene_B]]
                top_third <- top_thirds[[gene_B]]
                
                
                essentiality_bottom <-
                  x |>  dplyr::filter(DepMap_ID %in% bottom_third$DepMap_ID) |> select({
                    {
                      gene_A
                    }
                  }) |>  unlist() |> as.numeric()
                essentiality_top <-
                  x |>  dplyr::filter(DepMap_ID %in% top_third$DepMap_ID) |> select({
                    {
                      gene_A
                    }
                  }) |> unlist() |> as.numeric()
                
                #Do the wilcoxon test
                if (length(top_third$DepMap_ID) > 0) {
                  wilcoxon_sl <-
                    wilcox.test(essentiality_bottom,
                                essentiality_top,
                                alternative = "greater")$p.value #Check whether essentiality_bottom is greater than essentiality_top (when gene y is missing the cell is more likely to die when gene x is removed)
                  
                }
                else {
                  wilcoxon_sl <- NA
                }
                #Do the hypergeometric depletion test
                
                if (gene_A %in% colnames(expression)) {
                  overlap <-
                    length(intersect(
                      bottom_third$DepMap_ID,
                      bottom_third_A$DepMap_ID
                    ))
                  
                  #phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE) (the depletion test)
                  
                  depletion_p_value <-
                    phyper(
                      overlap,
                      length(bottom_third$DepMap_ID),
                      total_cell_lines - length(bottom_third$DepMap_ID),
                      length(bottom_third_A$DepMap_ID),
                      lower.tail = TRUE
                    )
                  
                }
                else {
                  depletion_p_value <- 2 #fesho millorrrrrrrrrrrrrrrrrrrrrr
                }
                
                if (!is.nan(wilcoxon_sl)) {
                  SL_pairs_i <- SL_pairs_i |>
                    add_row(
                      gene1 = word(gene_A),
                      gene2 = word(gene_B),
                      p_value = wilcoxon_sl,
                      depletion_p_value = depletion_p_value
                    )
                }
                
              }
              SL_pairs_i
            }
            
          }
      }
    }
    
    
    ####Depletion test
    depletions <- numeric(nrow(SL_pairs))
    
    for (i in  1:nrow(SL_pairs)) {
      if (i %% 1000 == 0) {
        #p()
      }
      gene_A <- SL_pairs$gene1[i]
      gene_B <-  SL_pairs$gene2[i]
      
      bottom_third_A <- bottom_thirds[[gene_A]]
      bottom_third <- bottom_thirds[[gene_B]]
      
      if (gene_A %in% colnames(expression) && gene_B %in% colnames(expression)) {
        overlap <-
          length(intersect(bottom_third$DepMap_ID,
                           bottom_third_A$DepMap_ID))
        
        #phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE) (the depletion test)
        
        depletion_p_value <-
          phyper(
            overlap,
            length(bottom_third$DepMap_ID),
            total_cell_lines - length(bottom_third$DepMap_ID),
            length(bottom_third_A$DepMap_ID),
            lower.tail = TRUE
          )
        
      } else {
        depletion_p_value <- 2 #fesho millorrrrrrrrrrrrrrrrrrrrrr
      }
      depletions[i] <- depletion_p_value
      
      
    }
    
    
    SL_pairs <-  SL_pairs |>  mutate(depletion_p_value = depletions)
    
    
    
    
    
    
    
    
    
    SL_pairs <-
      SL_pairs[order(SL_pairs$p_value, decreasing = FALSE),]
    
    q_value <- p.adjust(SL_pairs$p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(q_value = q_value)
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_p_value = 1 - SL_pairs$p_value)   #No estic segur que aixo estigui be
    
    SR_DD_q_value <-
      p.adjust(SL_pairs$SR_DD_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_q_value = SR_DD_q_value)
    
    
    depletion_q_value <-
      p.adjust(SL_pairs$depletion_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(depletion_q_value = depletion_q_value)
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_depletion_p_value = 1 - SL_pairs$depletion_p_value)
    
    SR_DD_depletion_q_value <-
      p.adjust(SL_pairs$SR_DD_depletion_p_value, method = "fdr")
    
    SL_pairs <-
      SL_pairs |> add_column(SR_DD_depletion_q_value = SR_DD_depletion_q_value)
    
    
    return(SL_pairs)
  }