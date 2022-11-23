#BiocManager::install("GEOquery") ###to install GEOquery package
library(stringr)
library(GEOquery)
library("XML")
library("methods")

get_targets <- function(targets){
  #used to get vector of  drug targetsfrom the metadata
  #parses strings separated by commas to a vector of drug targets
  #you would call it like this:
  #my_clinical_data <- get_clinical_data(path, datasetname)
  #metadata <- my_clinical_data$metadata
  #drug_targets <- get_targets(metadata$drug_targets)
  #drug_targets now contains a vector of drug target genes of you drug
  
  sep <- strsplit(targets, ",")
  vec <- c()
  for (word in sep){
    vec <- append(vec, word)
    
  }
  return(vec)
  
}
get_clinical_data <- function(pth, dataset){
  #takes in path to "Data" folder, dataset name
  #returns a list of 3 dataframes: a gene expression matrix,  a response dataframe and a metadata dataframe
  if (dataset == "GSE109211"){
    ###THIS ONE HAS DUPLICATE GENES
    gset <- getGEO("GSE109211", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL13938", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    
    
    
    dataset <- "GSE109211"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE109211_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    end-start-2
    
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    
    rownames(ex) <- df_family$Symbol
    
    data <- xmlParse("./family/GSE109211_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    sample <- c()
    response <- c()
    no_response <- c()
    #print( xmlValue(root[[9]][[6]][[5]], trim=TRUE) == "non-responder")
    for (i in 0:len){
      treat <- xmlValue(root[[9+i]][[6]][[4]], trim=TRUE)
      resp <- xmlValue(root[[9+i]][[6]][[5]], trim=TRUE)
      if(treat == "Sor"){
        sample <- append(sample, colnames(ex)[i])
        if (resp=="responder"){
          response <- append(response, 1)
          no_response <- append(no_response, 0)
        }
        else{
          response <- append(response, 0)
          no_response <- append(no_response, 1)      
        }
      }
    }
    
    ex <- subset( ex, select = sample )
    res <- data.frame(sample, response, no_response)
    
    study_name <- c("GSE109211")
    number_of_samples <- c(67)
    number_of_responders <- c(46)
    number_of_nonresponders <- c(21)
    drug <- c("Sorafenib")
    target_genes<-c("BRAF,RAF1,FLT4,KDR,FLT3,PDGFRB,KIT,FGFR1,RET,FLT1")
    cancer_type <- c("LIHC")
    treatment_type <- c("targeted")
    platform <- c("microarray")
    criteria_for_response <- c("All 70 RFS(recurrence-free-survival) events were recurrences, thus time to recurrence equalled RFS")
    sample <- c("FFPE")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    #head(metadata)
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)
    
  }
  else if (dataset == "GSE68871"){
    gset <- getGEO("GSE68871", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    
    #head(ex)
    dataset <- "GSE68871"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE68871_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    
    
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    
    #print(c(nrow(df_family), ncol(df_family)))
    rownames(ex) <- df_family$Gene.Symbol
    
    ex <- ex[row.names(ex) != "", ,drop = FALSE]
    #print(c(nrow(ex), ncol(ex)))
    
    data <- xmlParse("./family/GSE68871_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    sample <- c()
    response <- c()
    no_response <- c()
    CR <- c()
    nCR <- c()
    VGPR <- c()
    SD <- c()
    PR <- c()
    #root[[7]]
    #xmlValue(root[[7]][[6]][[5]], trim=TRUE)=="nCR"
    #length(colnames(ex))
    responders <- 0
    for (i in 0:len){
      sample <- append(sample, colnames(ex)[i+1])
      resp <- xmlValue(root[[7+i]][[6]][[5]], trim=TRUE)
      if(resp == "CR" | resp == "nCR"){
        
        response <- append(response, 1 )
        no_response <- append(no_response, 0)
        if(resp == "CR"){
          CR <- append(CR, 1)
          nCR <- append(nCR, 0)
          VGPR <- append(VGPR, 0)
          SD <- append(SD, 0)
          PR <- append(PR, 0)
        }
        else if (resp == "nCR"){
          CR <- append(CR, 0)
          nCR <- append(nCR, 1)
          VGPR <- append(VGPR, 0)
          SD <- append(SD, 0)
          PR <- append(PR, 0)
        }
        responders <- responders + 1
      }
      else{
        response <- append(response, 0 )
        no_response <- append(no_response, 1)
        if(resp == "PR"){
          CR <- append(CR, 0)
          nCR <- append(nCR, 0)
          VGPR <- append(VGPR, 0)
          SD <- append(SD, 0)
          PR <- append(PR, 1)
        }
        else if(resp == "SD"){
          CR <- append(CR, 0)
          nCR <- append(nCR, 0)
          VGPR <- append(VGPR, 0)
          SD <- append(SD, 1)
          PR <- append(PR, 0)
        }
        else if(resp == "VGPR"){
          CR <- append(CR, 0)
          nCR <- append(nCR, 0)
          VGPR <- append(VGPR, 1)
          SD <- append(SD, 0)
          PR <- append(PR, 0)
        }
        else{
          print("mistake!")
        }
      }
    }
    
    res <- data.frame(sample, response, no_response, CR, nCR, PR, SD, VGPR)
    
    study_name <- c("GSE68871")
    number_of_samples <- c(118)
    number_of_responders <- c(29)
    number_of_nonresponders <- c(89)
    drug <- c("Bortezomib/Thalidomide/Dexamethasone")
    target_genes<-c("PSMB5,PSMB1,NOS2")
    cancer_type <- c("MM")
    treatment_type <- c("targeted")
    platform <- c("microarray")
    criteria_for_response <- c("CR, nCR (as reported by European Group for Blood and Marrow Transplantation)")
    sample <- c("NA")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)
  }
  else if(dataset == "GSE119262"){
    ###NOTE: expression matrix has 22177 rows, but family has 22185. Had to use for loop to change rownames
    gset <- getGEO("GSE119262", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL6104", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    #head(ex)
    #nrow(ex)
    dataset <- "GSE119262"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE119262_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    #head(df_family)
    #"ILMN_1815924"
    #tail(rownames(ex))
    #nrow(ex)
    #nrow(df_family)
    for(rowname in rownames(ex)){
      rownames(ex)[rownames(ex) == rowname] <- df_family[which(df_family$ID == rowname), "Symbol"]
    }
    
    #tail(ex)
    #length(which(row.names(ex) == "NUP107"))
    
    
    data <- xmlParse("./family/GSE119262_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    sample <- c()
    response <- c()
    no_response <- c()
    delta_ki67 <- c()
    #xmlValue(root[[5]][[6]][[10]], trim=TRUE)
    #xmlValue(root[[5]][[6]][[9]], trim=TRUE)
    #a <- xmlValue(root[[5]][[2]], trim=TRUE)
    #str_detect(a, "Post", negate = FALSE)
    for (i in 0:len){
      status <- xmlValue(root[[5+i]][[2]], trim=TRUE)
      resp <- xmlValue(root[[5+i]][[6]][[10]], trim=TRUE)
      dki67 <- as.double(xmlValue(root[[5+i]][[6]][[9]], trim=TRUE))
      if(str_detect(status, "Pre", negate = FALSE)){
        sample <- append(sample, colnames(ex)[i+1])
        if(resp == "Responder_Ki67"){
          response <- append(response, 1)
          no_response <- append(no_response, 0)
        }
        else if(resp == "Non-Responder_Ki67"){
          response <- append(response, 0)
          no_response <- append(no_response, 1)
        }
        else{
          print("mistake")
        }
        
        delta_ki67 <- append(delta_ki67, dki67)
      }
    }
    #print(length(colnames(ex)))
    #print(length(sample))
    res <- data.frame(sample, response, no_response, delta_ki67)
    ex <- subset( ex, select = sample )
    
    study_name <- c("GSE119262")
    number_of_samples <- c(46)
    number_of_responders <- c(28)
    number_of_nonresponders <- c(18)
    drug <- c("Everolimus")
    target_genes<-c("MTOR")
    cancer_type <- c("BRCA")
    treatment_type <- c("targeted")
    platform <- c("microarray")
    criteria_for_response <- c("delta-ki67<top tertile")
    sample <- c("Frozen")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)
  }
  
  else if ( dataset == "GSE41998"){
    
    gset <- getGEO("GSE41998", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    
    #tail(ex)
    #nrow(ex)
    dataset <- "GSE41998"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE41998_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    #nrow(df_family)
    #nrow(ex)
    #length(which(df_family$Gene.Title == ""))
    #nrow(df_family)
    rownames(ex) <- df_family$Gene.Symbol
    ex <- ex[row.names(ex) != "", ,drop = FALSE]
    #ail(ex)
    #nrow(ex)
    
    
    data <- xmlParse("./family/GSE41998_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    sample <- c()
    response <- c()
    no_response <- c()
    treatment <- c()
    zero <- c()
    complete_response <- c()
    partial_response <- c()
    stable_disease <- c()
    progressive_disease <- c()
    
    #root[[7]]
    #xmlValue(root[[7]][[6]][[7]], trim=TRUE)# response
    #xmlValue(root[[7]][[6]][[6]], trim=TRUE)#treatment
    #root[[7]]
    #responders <- 0
    #nr <- 0
    for (i in 0:len){
      sample <- append(sample, colnames(ex)[i+1])
      resp <- xmlValue(root[[7+i]][[6]][[7]], trim=TRUE)
      treat <- xmlValue(root[[7+i]][[6]][[6]], trim=TRUE)
      treatment <- append(treatment, treat)
      if(resp == "complete response" | resp == "partial response"){
        
        response <- append(response, 1 )
        no_response <- append(no_response, 0)
        if(resp == "complete response"){
          zero <- append(zero, 0)
          complete_response <- append(complete_response, 1)
          partial_response <- append(partial_response, 0)
          stable_disease <- append(stable_disease, 0)
          progressive_disease <- append(progressive_disease, 0)
        }
        else if(resp == "partial response"){
          zero <- append(zero, 0)
          complete_response <- append(complete_response, 0)
          partial_response <- append(partial_response, 1)
          stable_disease <- append(stable_disease, 0)
          progressive_disease <- append(progressive_disease, 0)
        }
        #responders <- responders + 1
      }
      else{
        response <- append(response, 0 )
        no_response <- append(no_response, 1)
        if (resp == "progressive disease"){
          zero <- append(zero, 0)
          complete_response <- append(complete_response, 0)
          partial_response <- append(partial_response, 0)
          stable_disease <- append(stable_disease, 0)
          progressive_disease <- append(progressive_disease, 1)
        }
        else if(resp == "stable disease"){
          zero <- append(zero, 0)
          complete_response <- append(complete_response, 0)
          partial_response <- append(partial_response, 0)
          stable_disease <- append(stable_disease, 1)
          progressive_disease <- append(progressive_disease, 0)
        }
        else{
          
          zero <- append(zero, 1)
          complete_response <- append(complete_response, 0)
          partial_response <- append(partial_response, 0)
          stable_disease <- append(stable_disease, 0)
          progressive_disease <- append(progressive_disease, 0)
          
        }
        
        #nr <- nr+1
        
        
      }
      
    }
    #ncol(ex)
    res <- data.frame(sample, response, no_response, treatment, complete_response, partial_response, stable_disease, progressive_disease, zero)
    #length(which(res$no_response == 0))
    
    study_name <- c("GSE41998")
    number_of_samples <- c(279)
    number_of_responders <- c(201)
    number_of_nonresponders <- c(78)
    ### mutually exclusive drug treatments
    drug <- c("Doxorubicin or Cyclophosphamide or Ixabepilone or Paclitaxel")
    target_genes<-c("TUBB3,TOP2A,TUBB1,BCL2")
    cancer_type <- c("BRCA")
    treatment_type <- c("cytotoxic")
    platform <- c("microarray")
    criteria_for_response <- c("complete response or partial response")
    sample <- c("Frozen")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    
    
    #######NOTE : res has samples where the treatment is "none" also some outcomes are indeterminate (see rows where "zer0" column is 1). From spreadsheet, Lee seems to have classfiied this as no response
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)
  }
  
  else if(dataset == "GSE10255"){
    
    gset <- getGEO("GSE10255", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    #nrow(ex)
    #head(ex)
    
    dataset <- "GSE10255"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE10255_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    
    
    
    
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    #df_family[which(df_family$ID == "1007_s_at"), "Gene.Symbol"]
    
    
    
    rownames(ex) <- df_family$Gene.Symbol
    ex <- subset(ex, row.names(ex) != "")
    
    data <- xmlParse("./family/GSE10255_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    sample <- c()
    response <- c()
    no_response <- c()
    WBCdelta3 <- c()
    #root[[16]]
    #xmlValue(root[[16]][[6]][[5]], trim=TRUE)# WBCdelta3
    
    
    for (i in 0:len){
      sample <- append(sample, colnames(ex)[i+1])
      WBCdelta3 <- append(WBCdelta3 ,as.double(xmlValue(root[[16+i]][[6]][[5]], trim=TRUE)))
    }
    
    med <- median(WBCdelta3)
    
    for (val in WBCdelta3){
      if(val < med){
        response <- append(response, 1)
        no_response <- append(no_response, 0)
      }
      else{
        response <- append(response, 0)
        no_response <- append(no_response, 1)
      }
    }
    
    res <- data.frame(sample, response, no_response, WBCdelta3)
    
    study_name <- c("GSE10255")
    number_of_samples <- c(161)
    number_of_responders <- c(80)
    number_of_nonresponders <- c(81)
    drug <- c("Methotrexate")
    target_genes<-c("DHFR,TYMS")
    cancer_type <- c("ALL (acute lymphoblastic leukemia)")
    treatment_type <- c("cytotoxic")
    platform <- c("microarray")
    criteria_for_response <- c("WBCdelta3<median")
    sample <- c("Frozen")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)
    
  }
  
  
  
  else if(dataset == "GSE3964"){
    # load series and platform data from GEO
    ######THIS TRIAL DOES NOT SEEM TO GIVE LOG2 VALUES IN GENE EX MATRIX
    gset <- getGEO("GSE3964", GSEMatrix =TRUE, getGPL=FALSE)
    if (length(gset) > 1) idx <- grep("GPL3282", attr(gset, "names")) else idx <- 1
    gset <- gset[[idx]]
    
    ex <- exprs(gset)
    # log2 transform
    qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
    LogC <- (qx[5] > 100) ||
      (qx[6]-qx[1] > 50 && qx[2] > 0)
    if (LogC) { ex[which(ex <= 0)] <- NaN
    ex <- log2(ex) }
    #head(ex)
    
    dataset <- "GSE3964"
    cur_path <- file.path(pth, dataset)
    setwd(cur_path)
    family_file <- "GSE3964_family.soft"
    start <- which(grepl("!platform_table_begin", readLines(family_file)))
    end <- which(grepl("!platform_table_end", readLines(family_file)))
    df_family<- read.delim(family_file, header = TRUE, sep = "\t", skip = start, nrows= end-start-2)
    #length(which(df_family$SYMBOL == ""))
    rownames(ex) <- df_family$SYMBOL
    ex <- subset(ex, row.names(ex) != "")
    #nrow(ex)
    #ncol(ex)
    
    data <- xmlParse("./family/GSE3964_family.xml")
    
    root  <- xmlRoot(data)
    
    len <- length(colnames(ex))-1
    
    sample <- c()
    response <- c()
    no_response <- c()
    #root[[16]]
    #str_detect(xmlValue(root[[16]][[2]], trim = TRUE), "Cy3")
    #xmlValue(root[[16]][[6]][[9]], trim = TRUE)
    
    for (i in 0:len){
      if(str_detect(xmlValue(root[[16+i]][[2]], trim = TRUE), "Cy3") | str_detect(xmlValue(root[[16+i]][[2]], trim = TRUE), "CY3")){
        sample <- append(sample, colnames(ex)[i+1])
        if(xmlValue(root[[16+i]][[6]][[9]], trim = TRUE)=="Sensitive (S)"){
          response <- append(response, 1)
          no_response <- append(no_response, 0)
        }
        else{
          response <- append(response, 0)
          no_response <- append(no_response, 1)
          
        }
        
      }
      
      
    }
    ex <- subset( ex, select = sample )
    res <- data.frame(sample, response, no_response)
    #length(which(res$response == 1))
    
    study_name <- c("GSE3964")
    number_of_samples <- c(29)
    number_of_responders <- c(14)
    number_of_nonresponders <- c(15)
    drug <- c("FA/5FU/Irinotecan")
    target_genes<-c("FOLR1,FOLR2,FOLR3,TYMS,TOP2A")
    cancer_type <- c("COADREAD")
    treatment_type <- c("cytotoxic")
    platform <- c("microarray")
    criteria_for_response <- c("Initial (primary) response rates were assessed after each series of two treatment cycles based on WHO response criteria [8], considering complete or partial regression, stabilization or progression of the disease")
    sample <- c("Frozen")
    metadata <- data.frame(study_name,number_of_samples, number_of_responders, number_of_nonresponders,
                           drug, target_genes, cancer_type, treatment_type, platform,criteria_for_response, sample  )
    
    out <- list()
    out$ex <- ex
    out$res <- res
    out$metadata <- metadata
    return(out)    
  }
  else{
    print("ERROR: Not a valid Dataset name")
    return(list())
  }
  
}


#########TO RUN THIS FUNCTIONS
#directory <- ###YOUR DIRECTORY TO THE "Data" FOLDER HERE
#get_clinical_data(directory, ###YOUR DATASET NAME HERE (e.g., "GSE3964"))
#?read.delim


