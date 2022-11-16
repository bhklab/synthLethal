# Synthetic Lethality project

datadir <- "/Users/iansmith/Work/bhk/data/depmap"

# Load data
expr_ds <- read.delim(file.path(datadir, "CCLE_expression.csv"), header=TRUE, sep=",")
essen_ds <- read.delim(file.path(datadir, "CRISPR_gene_effect.csv"), header=TRUE, sep=",")

meta_ds <- read.delim(file.path(datadir, "Achilles_metadata.csv"), header=TRUE, sep=",")

common_cells <- intersect(essen_ds$DepMap_ID, expr_ds$X)
expr_ds <- expr_ds[match(common_cells, expr_ds$X),]
essen_ds <- essen_ds[match(common_cells, essen_ds$DepMap_ID),]

expr_genes <- sapply(strsplit(colnames(expr_ds), "..", fixed=TRUE), FUN=function(x) x[[1]])
essen_genes <- sapply(strsplit(colnames(essen_ds), "..", fixed=TRUE), FUN=function(x) x[[1]])


plot(density(as.matrix(essen_ds[, 2:10]), bw=0.01), xlab="Essentiality", ylab="Density", main="Essentiality Scores of 100 genes", lwd=3)

plot(density(as.matrix(expr_ds[, 1200:2255]), bw=0.01), xlab="Expression", ylab="Density", main="CCLE Expression", lwd=3, col="red")