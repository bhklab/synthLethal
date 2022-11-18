library(PharmacoGx)

# PSet directory, set to location of PSet RDS files:
psetdir <- "."

ctrp <- readRDS(file.path(psetdir, "CTRPv2_2015.rds"))

# Accessors for metadata about the PSet:

drugInfo(ctrp)
drugNames(ctrp)

cellInfo(ctrp)
cellNames(ctrp)


# List the molecular data types:
mDataNames(ctrp)

# List types of drug sensitivity Measures:
sensitivityMeasures(ctrp)


# Extract metadata for convenience:
ctrp_drugs <- drugInfo(ctrp)

ctrp_aacs <- summarizeSensitivityProfiles(ctrp, sensitivity.measure="aac_recomputed")




# Run pharmaco_train
gdsc <- readRDS("GDSC_2020(v2-8.2).rds")

# Get available sensitivity measures
sensitivityMeasures(gdsc)

# Get drug info
mydrugs <- drugInfo(gdsc)

# Use aac_recomputed by default
# Rows are drugs, columns are cell lines
aacs <- summarizeSensitivityProfiles(gdsc, sensitivity.measure="aac_recomputed")


# List available mDataNames
mDataNames(gdsc)

# Get gene expression data - maybe use microarray instead of RNA-seq
mps <- summarizeMolecularProfiles(gdsc, mDataType="Kallisto_0.46.1.rnaseq")

fdata <- assay(mps)

mygeneinfo <- featureInfo(gdsc, mDataType = "Kallisto_0.46.1.rnaseq")

# Slice gene expression matrix to just protein coding genes
# This is gross, clean up
fdata <- fdata[rownames(fdata) %in% rownames(mygeneinfo)[mygeneinfo$gene_type == "protein_coding"], ]
rownames(fdata) <- mygeneinfo$gene_name[match(rownames(fdata), rownames(mygeneinfo))]


# mydrug <- "Sirolimus"
# mytarget <- "MTOR"
mydrug <- "Sorafenib"
mytarget <- "BRAF"

if (length(grep(mydrug, rownames(aacs), ignore.case=TRUE)) != 1){
  print("Warning: number of drugs matched does not equal 1.")
}

yall <- aacs[grep(mydrug, rownames(aacs), ignore.case=TRUE), ]

# Remove NAs - cell lines for which the drug was not assayed, or cell lines without rna-seq data
missingcells <- colSums(is.na(fdata))

x <- fdata[, (complete.cases(yall) & missingcells == 0)] 
y <- yall[(complete.cases(yall) & missingcells == 0)]

# Plug into model training
source("pmodel_train.R")
my_working_directory <- 'some_path'
md_en <- pmodel_train(x,y, mydrug, "elasticnet", mytarget, my_working_directory, 'GSE119262')
# md_en <- pharmcomodel_train(x,y, mydrug, model, mytarget)

my_data <- get_clinical_data(paste(my_working_directory, "Data", sep="/"),"GSE119262")

source("pmodel_apply.R")
# Apply models
pmodel_apply(my_data$ex, my_data$res, md_en, drug="", modeltype="elasticnet")

#x_clinical <- ReadClinicalData
# Fix clinical data
pharmacomodel_apply(x_clinical, md_select, model="select")




gdsc_x <- get_genex_pset(gdsc)
gdsc_y <- get_response_pset(gdsc, "Sirolimus")

mymodel <- pmodel_train(gdsc_x, gdsc_y, "Sirolimus", modeltype="cows")
gse123ds <- get_clinical_data(datapath, "GSE123")

res_ctrp <- pmodel_apply(ctrp$exprmat, ctrp...)
res <- pmodel_apply(get123ds$exprmat, get123ds$response, "...")