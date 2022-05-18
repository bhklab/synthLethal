download_data <- function(path){
  ###The "path" parameter is whatever directory you want the "Data" folder to be created inn. The "Data" folder will house subfolders for each of the clinical trials. 
  ###function returns the path to the data folder
  
  datapath <- file.path(path, "Data")
  dir.create(datapath, showWarnings = FALSE)
  ###get data for GSE109211
  out_dir <- file.path(datapath, "GSE109211")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109211/matrix/GSE109211_series_matrix.txt.gz", destfile="GSE109211_series_matrix.txt.gz")
  R.utils::gunzip("GSE109211_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109211/miniml/GSE109211_family.xml.tgz", destfile="GSE109211_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE109211_family.xml.tgz", exdir = family)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109211/soft/GSE109211_family.soft.gz", destfile="GSE109211_family.soft.gz")
  R.utils::gunzip("GSE109211_family.soft.gz", overwrite=TRUE)
  
  
  ###get data for GSE68871
  out_dir <- file.path(datapath, "GSE68871")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68871/matrix/GSE68871_series_matrix.txt.gz", destfile="GSE68871_series_matrix.txt.gz")
  R.utils::gunzip("GSE68871_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68871/miniml/GSE68871_family.xml.tgz", destfile="GSE68871_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE68871_family.xml.tgz", exdir = family)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68871/soft/GSE68871_family.soft.gz", destfile="GSE68871_family.soft.gz")
  R.utils::gunzip("GSE68871_family.soft.gz", overwrite=TRUE)
  
  ###get data for GSE119262
  out_dir <- file.path(datapath, "GSE119262")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119262/matrix/GSE119262_series_matrix.txt.gz", destfile="GSE119262_series_matrix.txt.gz")
  R.utils::gunzip("GSE119262_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119262/miniml/GSE119262_family.xml.tgz", destfile="GSE119262_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE119262_family.xml.tgz", exdir = family)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119262/soft/GSE119262_family.soft.gz", destfile="GSE119262_family.soft.gz")
  R.utils::gunzip("GSE119262_family.soft.gz", overwrite=TRUE)
  
  ###get data for GSE3964
  out_dir <- file.path(datapath, "GSE3964")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3964/matrix/GSE3964_series_matrix.txt.gz", destfile="GSE3964_series_matrix.txt.gz")
  R.utils::gunzip("GSE3964_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3964/miniml/GSE3964_family.xml.tgz", destfile="GSE3964_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE3964_family.xml.tgz", exdir = family)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3964/soft/GSE3964_family.soft.gz", destfile="GSE3964_family.soft.gz")
  R.utils::gunzip("GSE3964_family.soft.gz", overwrite=TRUE)
  
  ###get data for GSE41998
  
  out_dir <- file.path(datapath, "GSE41998")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41998/matrix/GSE41998_series_matrix.txt.gz", destfile="GSE41998_series_matrix.txt.gz")
  R.utils::gunzip("GSE41998_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41998/miniml/GSE41998_family.xml.tgz", destfile="GSE41998_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE41998_family.xml.tgz", exdir = family)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE41nnn/GSE41998/soft/GSE41998_family.soft.gz", destfile="GSE41998_family.soft.gz")
  R.utils::gunzip("GSE41998_family.soft.gz", overwrite=TRUE)
  
  ###get data for GSE10255
  
  out_dir <- file.path(datapath, "GSE10255")
  dir.create(out_dir, showWarnings = FALSE)
  setwd(out_dir)
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10255/matrix/GSE10255_series_matrix.txt.gz", destfile="GSE10255_series_matrix.txt.gz")
  R.utils::gunzip("GSE10255_series_matrix.txt.gz", overwrite=TRUE)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10255/miniml/GSE10255_family.xml.tgz", destfile="GSE10255_family.xml.tgz")
  family <- file.path(out_dir, "family")
  dir.create(family, showWarnings = FALSE)
  untar("GSE10255_family.xml.tgz", exdir = family)
  
  download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE10nnn/GSE10255/soft/GSE10255_family.soft.gz", destfile="GSE10255_family.soft.gz")
  R.utils::gunzip("GSE10255_family.soft.gz", overwrite=TRUE)
  return(datapath) 
}


###RUN BELOW TO DOWNLOAD DATA
#directory <- ##your directory here
#download_data(directory)
