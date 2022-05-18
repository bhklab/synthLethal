load_pdataset <- function(datapath, dataset){
  
}



download.file(url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109211/suppl/GSE109211_Non-normalized_data.txt.gz", destfile="GSE109211_nnd.txt.gz")
gunzip("GSE109211_nnd.txt.gz")
download.file(url="https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109211&format=file", destfile="GSE109211_RAW.tar")
untar(tarfile = "GSE109211_RAW.tar")
gunzip("GPL13938_HumanHT-12_V4_0_R2_15002873_B_WGDASL.txt.gz")
