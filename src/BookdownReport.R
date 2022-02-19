#!/usr/bin/env Rscript
#BookdownReport.R

args = commandArgs(trailingOnly=T)
if(length(args)<2){
  message("'EAinit' can be used to create a config file for an RNAseq project")
  stop("config yaml file is required!")
}
message("Loading resource ...")
suppressMessages(source(paste0(args[1],"utility.R"),chdir=T))
config <- sapply(yaml::read_yaml(args[3]),unlist)
sysConfig <- yaml::read_yaml(paste0(args[1],"sys.yml"))

#Initiate a bookdown document
message("Rendering book ...")
suppressMessages(invisible(capture.output(
  bookdown::render_book(args[2]), 
  file = NULL)))
