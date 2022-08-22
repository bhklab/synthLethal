library(openxlsx)

achilles_metadata <- read.delim("Achilles_metadata.csv", header=TRUE, sep=",")

ccle_expr <- read.delim("CCLE_expression.csv")

genenames <- colnames(ccle_expr)

filter_gene <- strsplit(genenames, "..", fixed=TRUE)
genenames_clean <- sapply(filter_gene, FUN=function(x) x[[1]])

bimodtable <- read.xlsx("wail_bimodal_genelist_can-21-2395.xlsx",startRow = 4, colNames = TRUE)

bimodgenes <- bimodtable$X2[2:length(bimodtable$X2)]

