###load packages
library(dbscan)
library(Seurat)
library(ggplot2)
library(optparse)
library(ggpubr)
library(stringr)
library(dplyr)
library(parallel)
library(MAESTRO)
library(purrr)

###set options
option_list <- list(
  make_option(c("-i","--inputFile"), type = "character", help = "path of meta file include rds path and batch information"),
  make_option(c("-o","--outPath"), type = "character", help = "output path"),
  make_option(c("-t","--threads"), type = "integer", help = "number of threads")
)

opt_parser <- OptionParser(option_list = option_list);
opts <- parse_args(opt_parser);

###set parameters
inputFile <- opts$inputFile
outPath <- opts$outPath
threads <- opts$threads
scriptPath <- opts$script_path

###import functions
source("4_QC/basic_qc_function.R")
source("4_QC/batch_function.R")

###Read in meta data
batchInfo <- read.table(inputFile, sep='\t', stringsAsFactors = F, header = T)
file_path <- batchInfo$rds_path

###batch mod is for evaluating the batch effect in data and removing it
batch_module <- function(file_path, fun, threads) {
  batch_info_list <- lapply(file_path, fun)
  if (length(batch_info_list)==1){
  batch_info <- t(as.data.frame(batch_info_list[[1]]))
} else{
  batch_info <- batch_info_list %>% reduce(rbind)
}
  colnames(batch_info) <- c("Sample","batch_dectect","sample_path")
  write.table(batch_info[,c("Sample","batch_dectect")],file = paste0(outPath,"/batch_detect.txt"),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
  return(batch_info[,"sample_path"])
}

###basic mod is for generating plots indicating the basic quality(nfeature,ncount,ncell,ncluster...) of data and will summarize the statistics into a table 
basicQC_module <- function(file_path, fun, threads) {
  # basic_info_list <- mclapply(file_path, fun, threads)
  basic_info_list <- lapply(file_path, fun)
  if (length(basic_info_list)==1){
    basic_info <- t(as.data.frame(basic_info_list[[1]]))
  } else{
    basic_info <- basic_info_list %>% reduce(rbind)
  }
  colnames(basic_info) <- c("Sname","ncell","ncluster","nfeature_median","nCount_median","up_mean","Ncell","Nfeature","Ncount","Data_quality")
  write.table(basic_info,file = paste0(outPath,"/basic_info.txt"),col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t")
}

###main function
Abatch_path <- batch_module(file_path, batch_main, threads)
basicQC_module(Abatch_path, basic_qc_main, threads)
